# Copyright 2026 Steven Kearnes
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Score every representation against the teacher's bioactivity similarity.

The headline metric is deliberately *not* overall agreement with the teacher.
Structure and bioactivity are already correlated, so a plain structural
fingerprint scores well on a global average without knowing anything about
biology, and a student that merely relearns ECFP would look successful.

What makes a bioactivity fingerprint worth having is where it *disagrees* with
structural similarity:

- **scaffold-hop regime** (ECFP4 < ``HOP_MAX``): compounds the teacher relates
  but structure does not. Recovering these is the entire point.
- **cliff regime** (ECFP4 > ``CLIFF_MIN``): near-identical compounds whose
  activities diverge. Here a good method must *disagree* with ECFP.

A representation is only interesting if it beats ECFP inside those regimes, so
ECFP is reported as a control on exactly the same restricted pair sets.
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

import fingerprints
import student
import teacher as teacher_mod

HOP_MAX = 0.35
CLIFF_MIN = 0.60

N_WORKERS = 12
CHUNK = 10_000

# Conformer generation costs ~10-50 ms per molecule, so the 3D descriptors are
# scored on a subsample of pairs rather than the whole evaluation set.
N_3D_SUBSAMPLE = 40_000

# Reference-library size for the neighbour control, bounded by memory.
REFERENCE_CAP = 300_000

# Training pairs for the distillation student.
N_DISTILL_PAIRS = 400_000

# A pair counts as bioactively related if its background-corrected teacher score
# clears this many standard deviations above size-matched random pairs.
RELATED_Z = 3.0


def regime_metrics(
    predicted: np.ndarray, teacher_z: np.ndarray, ecfp: np.ndarray
) -> dict[str, float]:
    """Agreement with the teacher, overall and split by structural regime."""
    related = teacher_z >= RELATED_Z

    def block(mask: np.ndarray) -> dict[str, float]:
        if mask.sum() < 100:
            return {"spearman": float("nan"), "auc": float("nan"), "n": int(mask.sum())}
        rho = spearmanr(predicted[mask], teacher_z[mask]).statistic
        labels = related[mask]
        # AUC needs both classes present; a regime can be all-related by
        # construction (every shared-target pair, say).
        auc = (
            roc_auc_score(labels, predicted[mask])
            if 0 < labels.sum() < len(labels)
            else float("nan")
        )
        return {"spearman": float(rho), "auc": float(auc), "n": int(mask.sum())}

    everything = np.ones(len(predicted), dtype=bool)
    return {
        "overall": block(everything),
        "scaffold_hop": block(ecfp < HOP_MAX),
        "cliff": block(ecfp > CLIFF_MIN),
    }


def _bits_chunk(args: tuple[list[str], str]) -> np.ndarray:
    chunk, kind = args
    return fingerprints.bit_fingerprints(chunk, kind)


def _counts_chunk(chunk: list[str]) -> np.ndarray:
    return fingerprints.count_fingerprints(chunk, "ecfp4")


def _usr_chunk(chunk: list[str]) -> np.ndarray:
    return fingerprints.usr_descriptors(chunk, catz=True)


def _descriptor_chunk(chunk: list[str]) -> np.ndarray:
    return fingerprints.descriptors_2d(chunk)


def parallel_stack(fn, smiles: list[str], extra=None) -> np.ndarray:
    blocks = [smiles[i:i + CHUNK] for i in range(0, len(smiles), CHUNK)]
    payload = [(b, extra) for b in blocks] if extra is not None else blocks
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        return np.vstack(list(pool.map(fn, payload)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eval", type=Path, required=True)
    parser.add_argument("--teacher", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--epochs", type=int, default=15)
    parser.add_argument(
        "--distill", action="store_true",
        help="also train the pair-similarity student (slower)",
    )
    args = parser.parse_args()

    pairs_frame = pd.read_parquet(args.eval / "pairs.parquet")
    compounds = pd.read_parquet(args.eval / "compounds.parquet")
    pairs = pairs_frame[["row_1", "row_2"]].to_numpy()
    teacher_z = pairs_frame["teacher_z"].to_numpy()
    ecfp_similarity = pairs_frame["ecfp4"].to_numpy()

    smiles = compounds["canonical_smiles"].tolist()
    is_train = (compounds["split"] == 0).to_numpy()
    print(f"compounds {len(smiles):,d}  pairs {len(pairs):,d}")

    raw_profiles = sparse.load_npz(args.teacher / "profiles_raw.npz")
    molregnos = pd.read_parquet(args.teacher / "profile_molregnos.parquet")["molregno"]
    profile_row = pd.Series(np.arange(len(molregnos)), index=molregnos.to_numpy())
    rows = compounds["molregno"].map(profile_row).to_numpy()
    profiles = raw_profiles[rows]

    results: dict[str, dict] = {}

    # Only compounds appearing in a pair need test-side features, and only
    # training compounds need supervision features. Featurising all 770k for
    # every representation would dominate the runtime for no benefit.
    used = np.unique(pairs)
    local = -np.ones(len(compounds), dtype=np.int64)
    local[used] = np.arange(len(used))
    local_pairs = local[pairs]
    used_smiles = [smiles[i] for i in used]
    print(f"compounds appearing in pairs: {len(used):,d}")

    # ---------------------------------------------------------------- baselines
    for kind in ("ecfp4", "ecfp6", "fcfp4", "rdkit", "atompair"):
        print(f"baseline {kind} ...", flush=True)
        fps = parallel_stack(_bits_chunk, used_smiles, extra=kind)
        predicted = fingerprints.tanimoto_pairs(fps, local_pairs)
        results[kind] = regime_metrics(predicted, teacher_z, ecfp_similarity)

    print("baseline rdkit2d ...", flush=True)
    descriptors = parallel_stack(_descriptor_chunk, used_smiles)
    # Descriptors span wildly different scales, so standardise before cosine or
    # molecular weight alone dominates the similarity.
    descriptors = (descriptors - descriptors.mean(0)) / (descriptors.std(0) + 1e-8)
    results["rdkit2d"] = regime_metrics(
        fingerprints.cosine_pairs(descriptors, local_pairs), teacher_z, ecfp_similarity
    )

    # 3D needs a conformer per molecule, so it is scored on a pair subsample.
    print("baseline usrcat (subsampled) ...", flush=True)
    rng = np.random.default_rng(0)
    subset = rng.choice(len(pairs), min(N_3D_SUBSAMPLE, len(pairs)), replace=False)
    sub_used = np.unique(pairs[subset])
    sub_local = -np.ones(len(compounds), dtype=np.int64)
    sub_local[sub_used] = np.arange(len(sub_used))
    usrcat = parallel_stack(_usr_chunk, [smiles[i] for i in sub_used])
    results["usrcat"] = regime_metrics(
        fingerprints.cosine_pairs(usrcat, sub_local[pairs[subset]]),
        teacher_z[subset],
        ecfp_similarity[subset],
    )

    # ------------------------------------------------------------------ students
    print("student: neighbour (Bioturbo control) ...", flush=True)
    train_rows = np.flatnonzero(is_train)
    # The reference library is capped for memory: as float32 the full 615k
    # training set is ~5 GB before the dense output array is allocated. 300k
    # annotated compounds is still a larger reference set than the original
    # Bioturbo work used, so the control is not weakened.
    if len(train_rows) > REFERENCE_CAP:
        train_rows = rng.choice(train_rows, REFERENCE_CAP, replace=False)
    train_bits = parallel_stack(
        _bits_chunk, [smiles[i] for i in train_rows], extra="ecfp4"
    )
    query_bits = parallel_stack(_bits_chunk, used_smiles, extra="ecfp4")
    inferred = student.neighbour_profiles(
        query_bits, train_bits, profiles[train_rows], n_neighbours=25
    )
    results["neighbour"] = regime_metrics(
        fingerprints.cosine_pairs(inferred, local_pairs), teacher_z, ecfp_similarity
    )

    print("student: multitask ...", flush=True)
    train_counts = parallel_stack(_counts_chunk, [smiles[i] for i in train_rows])
    model = student.train_multitask(
        train_counts, profiles[train_rows], n_epochs=args.epochs
    )
    del train_counts
    # Saved so the transfer benchmark can reuse the encoder without retraining.
    torch.save(
        {
            "state_dict": model.state_dict(),
            "n_input": fingerprints.N_BITS,
            "n_targets": profiles.shape[1],
        },
        args.out / "student_multitask.pt",
    )
    used_counts = parallel_stack(_counts_chunk, used_smiles)
    embedding = student.embed(model, used_counts)
    results["multitask"] = regime_metrics(
        fingerprints.cosine_pairs(embedding, local_pairs), teacher_z, ecfp_similarity
    )

    if args.distill:
        print("student: distill ...", flush=True)
        # Trained on TRAIN-side pairs only, so the test pairs it is scored on
        # were never seen; sampling them from the same distribution as the
        # evaluation would be leakage.
        train_rows_all = np.flatnonzero(is_train)
        left = rng.choice(train_rows_all, N_DISTILL_PAIRS)
        right = rng.choice(train_rows_all, N_DISTILL_PAIRS)
        keep = left != right
        distill_pairs = np.column_stack([left[keep], right[keep]])

        smoothed = sparse.load_npz(args.teacher / "profiles.npz")
        teacher_rows_all = compounds["molregno"].map(profile_row).to_numpy()
        target = teacher_mod.cosine_similarity(
            smoothed, teacher_rows_all[distill_pairs]
        ).astype(np.float32)

        distill_smiles = [smiles[i] for i in train_rows_all]
        distill_counts = parallel_stack(_counts_chunk, distill_smiles)
        compact = -np.ones(len(compounds), dtype=np.int64)
        compact[train_rows_all] = np.arange(len(train_rows_all))

        distill_model = student.train_distill(
            distill_counts, compact[distill_pairs], target, n_epochs=args.epochs
        )
        del distill_counts
        distill_embedding = student.embed(distill_model, used_counts)
        results["distill"] = regime_metrics(
            fingerprints.cosine_pairs(distill_embedding, local_pairs),
            teacher_z, ecfp_similarity,
        )

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "results.json").open("w") as fh:
        json.dump(results, fh, indent=2)

    table = pd.DataFrame({
        name: {
            f"{regime}_{metric}": values[regime][metric]
            for regime in ("overall", "scaffold_hop", "cliff")
            for metric in ("spearman", "auc")
        }
        for name, values in results.items()
    }).T
    print(table.round(3).to_string())
    table.to_csv(args.out / "results.csv")


if __name__ == "__main__":
    main()
