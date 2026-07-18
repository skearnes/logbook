"""Sample compound pairs for evaluation, with a scaffold-disjoint split.

Two sampling decisions carry the whole evaluation.

**Pairs are stratified, not random.** Random pairs drawn from 750k compounds are
almost all structurally unrelated *and* bioactively unrelated, so a metric
averaged over them is dominated by easy negatives and every method scores well.
The pairs that discriminate are the ones where structure and bioactivity
disagree, and those must be deliberately sampled: compounds sharing a target but
not a chemotype (scaffold hops), and near-identical compounds with divergent
activity (cliffs).

**The split is scaffold-disjoint, not random.** Martin et al. (JCIM 2017, 57,
2077) showed random splits massively overstate accuracy for models built on a
bioactivity matrix -- single-assay models scored r2 0.05 under a realistic split
versus far higher under a random one. Test compounds here share no Bemis-Murcko
scaffold with any training compound.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy import sparse

import fingerprints
import teacher as teacher_mod

RDLogger.DisableLog("rdApp.*")

# Fingerprinting and scaffold perception are per-molecule and independent, and
# there are ~770k molecules; running them on one core dominates the runtime.
N_WORKERS = 12
CHUNK = 20_000

TEST_FRACTION = 0.2
SEED = 20260718

# Pairs are drawn from three sources with different structure/activity
# relationships, so the evaluation spans the whole plane rather than its
# diagonal.
N_SHARED_TARGET = 300_000   # same target, any chemotype -> finds scaffold hops
N_STRUCTURAL = 300_000      # ECFP neighbours -> finds activity cliffs
N_RANDOM = 200_000          # background


def murcko_scaffold(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:
        return ""


def _scaffold_chunk(chunk: list[str]) -> list[str]:
    return [murcko_scaffold(s) for s in chunk]


def _fingerprint_chunk(chunk: list[str]) -> np.ndarray:
    return fingerprints.bit_fingerprints(chunk, "ecfp4")


def parallel_map(fn, smiles: list[str]) -> list:
    blocks = [smiles[i:i + CHUNK] for i in range(0, len(smiles), CHUNK)]
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        return list(pool.map(fn, blocks))


def scaffold_split(
    scaffolds_list: list[str], rng: np.random.Generator
) -> np.ndarray:
    """Assign each molecule to train (0) or test (1) by whole scaffold groups."""
    scaffolds: dict[str, list[int]] = {}
    for i, name in enumerate(scaffolds_list):
        scaffolds.setdefault(name, []).append(i)

    groups = list(scaffolds.values())
    rng.shuffle(groups)

    assignment = np.zeros(len(scaffolds_list), dtype=np.int8)
    target = int(TEST_FRACTION * len(scaffolds_list))
    taken = 0
    for group in groups:
        if taken >= target:
            break
        assignment[group] = 1
        taken += len(group)
    return assignment


def sample_shared_target(
    calls: pd.DataFrame, row_of: pd.Series, eligible: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Pairs of compounds active on a common target."""
    active = calls[calls["value"] > 0].copy()
    active["row"] = active["molregno"].map(row_of)
    active = active[active["row"].notna()]
    active["row"] = active["row"].astype(int)
    active = active[eligible[active["row"].to_numpy()]]

    pairs = []
    for _, group in active.groupby("tid", sort=False):
        rows = group["row"].to_numpy()
        if len(rows) < 2:
            continue
        # Cap per target: a handful of exhaustively-screened targets would
        # otherwise supply most of the evaluation set and bias it toward their
        # chemotypes.
        n_draw = min(len(rows) * 2, 400)
        left = rng.choice(rows, n_draw)
        right = rng.choice(rows, n_draw)
        pairs.append(np.column_stack([left, right]))

    stacked = np.vstack(pairs)
    stacked = stacked[stacked[:, 0] != stacked[:, 1]]
    if len(stacked) > N_SHARED_TARGET:
        stacked = stacked[rng.choice(len(stacked), N_SHARED_TARGET, replace=False)]
    return stacked


def sample_structural(
    fps: np.ndarray, eligible: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Pairs that are close in ECFP space, found by blocked brute force.

    The reference matrix is converted to float32 once and its bit counts cached;
    doing either inside the loop dominates the runtime, since the reference side
    is hundreds of thousands of rows wide.
    """
    rows = np.flatnonzero(eligible)
    anchors = rng.choice(rows, min(len(rows), 20_000), replace=False)

    reference = np.ascontiguousarray(fps[rows], dtype=np.float32)
    reference_counts = reference.sum(axis=1)[None, :]

    pairs = []
    block = 512
    for start in range(0, len(anchors), block):
        chunk = anchors[start:start + block]
        query = np.ascontiguousarray(fps[chunk], dtype=np.float32)
        intersection = query @ reference.T
        union = query.sum(axis=1)[:, None] + reference_counts - intersection
        similarity = np.divide(
            intersection, union, out=np.zeros_like(intersection), where=union > 0
        )
        # Self-matches sit at 1.0 and would otherwise dominate the top-k.
        similarity[similarity > 0.999] = 0.0
        top = np.argpartition(-similarity, 8, axis=1)[:, :8]
        pairs.append(np.column_stack([
            np.repeat(chunk, 8), rows[top.ravel()]
        ]))

    stacked = np.vstack(pairs)
    if len(stacked) > N_STRUCTURAL:
        stacked = stacked[rng.choice(len(stacked), N_STRUCTURAL, replace=False)]
    return stacked


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--teacher", type=Path, required=True)
    parser.add_argument("--curated", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    rng = np.random.default_rng(SEED)

    molregnos = pd.read_parquet(args.teacher / "profile_molregnos.parquet")["molregno"]
    calls = pd.read_parquet(args.teacher / "activity_calls.parquet")
    molecules = pd.read_parquet(args.curated / "molecules.parquet")

    frame = pd.DataFrame({"molregno": molregnos}).merge(molecules, on="molregno", how="left")
    frame = frame[frame["canonical_smiles"].notna()].reset_index(drop=True)
    print(f"profiled compounds with structures: {len(frame):,d}")

    smiles = frame["canonical_smiles"].tolist()
    print("computing ECFP4 ...", flush=True)
    fps = np.vstack(parallel_map(_fingerprint_chunk, smiles))

    print("scaffold split ...", flush=True)
    scaffolds = [s for block in parallel_map(_scaffold_chunk, smiles) for s in block]
    split = scaffold_split(scaffolds, rng)
    print(f"train {int((split == 0).sum()):,d}  test {int((split == 1).sum()):,d}", flush=True)

    row_of = pd.Series(np.arange(len(frame)), index=frame["molregno"])
    eligible = split == 1

    print("sampling pairs ...", flush=True)
    shared = sample_shared_target(calls, row_of, eligible, rng)
    structural = sample_structural(fps, eligible, rng)
    rows = np.flatnonzero(eligible)
    random_pairs = np.column_stack([
        rng.choice(rows, N_RANDOM), rng.choice(rows, N_RANDOM)
    ])
    random_pairs = random_pairs[random_pairs[:, 0] != random_pairs[:, 1]]

    pairs = np.vstack([shared, structural, random_pairs])
    source = np.concatenate([
        np.full(len(shared), "shared_target"),
        np.full(len(structural), "structural"),
        np.full(len(random_pairs), "random"),
    ])
    # The same pair can arise from two sources; keep one copy so it is not
    # double-counted in the metrics.
    order = np.lexsort((pairs[:, 1], pairs[:, 0]))
    pairs, source = pairs[order], source[order]
    _, unique = np.unique(pairs, axis=0, return_index=True)
    pairs, source = pairs[unique], source[unique]

    print(f"pairs: {len(pairs):,d}")

    profiles = sparse.load_npz(args.teacher / "profiles.npz")
    profile_row = pd.Series(np.arange(len(molregnos)), index=molregnos.to_numpy())
    teacher_rows = frame["molregno"].map(profile_row).to_numpy()

    teacher_pairs = np.column_stack([
        teacher_rows[pairs[:, 0]], teacher_rows[pairs[:, 1]]
    ])
    raw = teacher_mod.cosine_similarity(profiles, teacher_pairs)
    sizes = np.diff(profiles.indptr)
    corrected = teacher_mod.background_correct(
        raw, teacher_pairs, sizes, profiles, rng
    )
    ecfp = fingerprints.tanimoto_pairs(fps, pairs)

    out = pd.DataFrame({
        "row_1": pairs[:, 0],
        "row_2": pairs[:, 1],
        "source": source,
        "teacher_raw": raw,
        "teacher_z": corrected,
        "ecfp4": ecfp,
    })

    args.out.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out / "pairs.parquet", index=False)
    frame.assign(split=split).to_parquet(args.out / "compounds.parquet", index=False)
    np.save(args.out / "ecfp4.npy", fps)

    print(out.groupby("source")[["teacher_raw", "ecfp4"]].mean())
    print(f"wrote -> {args.out}")


if __name__ == "__main__":
    main()
