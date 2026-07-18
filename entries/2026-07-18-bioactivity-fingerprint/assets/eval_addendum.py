"""Follow-up scoring: the fair comparison, plus the repaired descriptor baseline.

Round one compared a 512-dimensional learned embedding against a 7,189-dimensional
inherited profile and the inherited profile won. That comparison confounds two
things -- how the representation was obtained, and how wide it is -- so it cannot
support a claim about distillation.

This scores the multi-task model's *profile head* instead: same target space and
same width as the neighbour control, differing only in whether the profile was
predicted by a trained network or borrowed from structural neighbours. That is
the comparison the conclusion actually rests on.
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import torch

import evaluate
import fingerprints
import student

N_WORKERS = 12
CHUNK = 10_000


def _counts_chunk(chunk: list[str]) -> np.ndarray:
    return fingerprints.count_fingerprints(chunk, "ecfp4")


def _descriptor_chunk(chunk: list[str]) -> np.ndarray:
    return fingerprints.descriptors_2d(chunk)


def parallel_stack(fn, smiles: list[str]) -> np.ndarray:
    blocks = [smiles[i:i + CHUNK] for i in range(0, len(smiles), CHUNK)]
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        return np.vstack(list(pool.map(fn, blocks)))


def predict_profiles(model: student.Encoder, features: np.ndarray, batch: int = 2048) -> np.ndarray:
    """The model's full predicted target profile, matching the teacher's width."""
    model.eval()
    chunks = []
    with torch.no_grad():
        for start in range(0, len(features), batch):
            block = torch.from_numpy(features[start:start + batch]).float()
            block = block.to(student.DEVICE)
            chunks.append(torch.tanh(model.predict_profile(block)).cpu().numpy())
    return np.vstack(chunks)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eval", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()

    pairs_frame = pd.read_parquet(args.eval / "pairs.parquet")
    compounds = pd.read_parquet(args.eval / "compounds.parquet")
    pairs = pairs_frame[["row_1", "row_2"]].to_numpy()
    teacher_z = pairs_frame["teacher_z"].to_numpy()
    ecfp_similarity = pairs_frame["ecfp4"].to_numpy()

    smiles = compounds["canonical_smiles"].tolist()
    used = np.unique(pairs)
    local = -np.ones(len(compounds), dtype=np.int64)
    local[used] = np.arange(len(used))
    local_pairs = local[pairs]
    used_smiles = [smiles[i] for i in used]
    print(f"scoring {len(used):,d} compounds over {len(pairs):,d} pairs", flush=True)

    with (args.results / "results.json").open() as fh:
        results = json.load(fh)

    print("rdkit2d (repaired) ...", flush=True)
    descriptors = parallel_stack(_descriptor_chunk, used_smiles)
    descriptors = (descriptors - descriptors.mean(0)) / (descriptors.std(0) + 1e-8)
    results["rdkit2d"] = evaluate.regime_metrics(
        fingerprints.cosine_pairs(descriptors, local_pairs), teacher_z, ecfp_similarity
    )
    del descriptors

    print("multitask profile head ...", flush=True)
    state = torch.load(args.checkpoint, map_location=student.DEVICE, weights_only=False)
    model = student.Encoder(state["n_input"], state["n_targets"]).to(student.DEVICE)
    model.load_state_dict(state["state_dict"])

    counts = parallel_stack(_counts_chunk, used_smiles)
    predicted = predict_profiles(model, counts)
    results["multitask_profile"] = evaluate.regime_metrics(
        fingerprints.cosine_pairs(predicted, local_pairs), teacher_z, ecfp_similarity
    )

    with (args.results / "results.json").open("w") as fh:
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
    table.to_csv(args.results / "results.csv")


if __name__ == "__main__":
    main()
