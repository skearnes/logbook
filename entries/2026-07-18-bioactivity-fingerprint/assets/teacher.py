"""Build the teacher: a bioactivity similarity between compounds from assay data.

The teacher is the thing a ligand-only fingerprint is being asked to imitate. It
is defined only for compounds ChEMBL has actually measured, which is exactly the
limitation that motivates distilling it into a structure-only student.

Construction, in order:

1. **Aggregate to comparability units.** Replicates are pooled only within a
   (compound, target, assay-conditions, potency-type) group; pairs whose
   measurements span more than ``MAX_SPREAD`` log units are dropped as
   irreconcilable. Values are never averaged across assay-condition groups --
   the recommendation of Landrum & Riniker (JCIM 2024, 64, 1560) -- so the
   best-evidenced group is kept per compound-target pair.

2. **Call activity on a graded scale.** A hard cut at pChEMBL 6 would be
   arbitrary given inter-assay MAE of 0.27-0.50 log units, so potency is mapped
   through a soft ramp and censored records enter as negatives at their bound.

3. **Smooth over related targets.** Roughly 70% of ChEMBL compounds carry a
   single target annotation, so raw target-set overlap is degenerate: two
   compounds either share their one target or score zero. Spreading each
   activity onto homologous and same-family targets makes those pairs
   comparable. The relatedness comes from protein sequence and curated protein
   family only -- never from ligand chemistry, which would inject structural
   similarity into a ground truth used to score a structure-only model.

4. **Correct for profile size.** Promiscuous and heavily-assayed compounds
   dominate any raw profile overlap. Following the logic of SEA (Keiser et al.,
   Nat. Biotechnol. 2007, 25, 197), a raw score is converted to a z-score
   against a background of random pairs matched on profile size, so "similar"
   means "more alike than two arbitrary compounds with these profiles".
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

# Measurements of one compound-target pair spanning more than this many log
# units cannot both be right; the group is discarded rather than averaged.
MAX_SPREAD = 1.0

# Potency below ACTIVE_P is unambiguously active, above INACTIVE_P
# unambiguously not; between them the call ramps linearly rather than snapping,
# because a hard threshold sits well inside the inter-assay noise band.
ACTIVE_P = 6.0
INACTIVE_P = 4.5

# Number of random pairs used to estimate the size-matched background.
BACKGROUND_PAIRS = 200_000

# Profile-size bins for background matching; a compound with 40 annotations has
# a very different score distribution from one with 2.
SIZE_BINS = (1, 2, 3, 5, 8, 13, 21, 34, 10**6)


def aggregate(actives: pd.DataFrame, inactives: pd.DataFrame) -> pd.DataFrame:
    """Collapse replicates to one graded activity call per compound-target pair."""
    group_keys = ["molregno", "tid", "standard_type", "ach"]

    stats = actives.groupby(group_keys, sort=False)["pchembl_value"].agg(
        ["median", "min", "max", "size"]
    )
    stats = stats[(stats["max"] - stats["min"]) <= MAX_SPREAD]
    stats = stats.reset_index()

    # Keep the best-evidenced assay-conditions group per compound-target pair
    # rather than averaging across groups that measured different things.
    stats = stats.sort_values("size").drop_duplicates(
        ["molregno", "tid"], keep="last"
    )
    positive = stats[["molregno", "tid"]].copy()
    positive["value"] = _grade(stats["median"].to_numpy())

    # Censored records state only a lower bound on the dissociation constant, so
    # they can establish inactivity but never activity. A bound weaker than
    # INACTIVE_P is uninformative -- ">100 nM" excludes nothing interesting.
    bounds = inactives.groupby(["molregno", "tid"], sort=False)["p_bound"].max()
    bounds = bounds[bounds <= INACTIVE_P].reset_index()
    negative = bounds[["molregno", "tid"]].copy()
    negative["value"] = -1.0

    combined = pd.concat([positive, negative], ignore_index=True)
    # An exact measurement outranks a censored bound for the same pair.
    return combined.drop_duplicates(["molregno", "tid"], keep="first")


def _grade(pchembl: np.ndarray) -> np.ndarray:
    """Map potency to [-1, 1] with a soft ramp across the ambiguous band."""
    ramp = (pchembl - INACTIVE_P) / (ACTIVE_P - INACTIVE_P)
    return (2.0 * np.clip(ramp, 0.0, 1.0) - 1.0).astype(np.float32)


def build_profiles(
    calls: pd.DataFrame, tids: np.ndarray, relatedness: sparse.csr_matrix
) -> tuple[sparse.csr_matrix, sparse.csr_matrix, np.ndarray]:
    """Return raw and relatedness-smoothed profile matrices, plus the compound index.

    Both are needed and they play different roles. The smoothed matrix defines
    the teacher's *similarity*, since relating targets is what rescues the ~70%
    of compounds carrying a single annotation. The raw matrix is the *supervision*
    for a student: it holds only what was actually measured, and training a model
    to reproduce smoothed values would be asking it to fit an interpolation
    artefact rather than data.
    """
    molregnos = np.sort(calls["molregno"].unique())
    mol_index = pd.Series(np.arange(len(molregnos)), index=molregnos)
    tid_index = pd.Series(np.arange(len(tids)), index=tids)

    rows = calls["molregno"].map(mol_index).to_numpy()
    cols = calls["tid"].map(tid_index).to_numpy()
    keep = ~pd.isna(cols)
    raw = sparse.csr_matrix(
        (calls["value"].to_numpy()[keep], (rows[keep], cols[keep].astype(int))),
        shape=(len(molregnos), len(tids)),
        dtype=np.float32,
    )

    # Targets are weighted by inverse document frequency. Sharing CYP3A4 or a
    # heavily-screened kinase says little, because half the library has been run
    # against it; sharing a rarely-assayed target is strong evidence of a common
    # mechanism. Without this, similarity mostly reports assay popularity.
    hits = np.asarray((raw != 0).sum(axis=0)).ravel()
    idf = np.log(raw.shape[0] / np.maximum(hits, 1)).astype(np.float32)
    weighted = raw @ sparse.diags(idf)

    # Each activity spreads onto related targets, so compounds annotated against
    # different members of one family are no longer orthogonal.
    smoothed = weighted @ relatedness
    return raw.tocsr(), smoothed.tocsr(), molregnos


def cosine_similarity(
    profiles: sparse.csr_matrix, pairs: np.ndarray, chunk: int = 50_000
) -> np.ndarray:
    """Cosine similarity for specific row pairs, without densifying the matrix.

    Smoothing leaves each profile with a few hundred nonzeros, so gathering all
    pair rows at once would materialise hundreds of millions of entries. Pairs
    are processed in chunks to keep that bounded.
    """
    norms = np.sqrt(profiles.multiply(profiles).sum(axis=1)).A.ravel()
    norms[norms == 0] = 1.0

    out = np.empty(len(pairs), dtype=np.float64)
    for start in range(0, len(pairs), chunk):
        block = pairs[start:start + chunk]
        left, right = profiles[block[:, 0]], profiles[block[:, 1]]
        dots = left.multiply(right).sum(axis=1).A.ravel()
        out[start:start + chunk] = dots / (norms[block[:, 0]] * norms[block[:, 1]])
    return out


def background_correct(
    raw: np.ndarray, pairs: np.ndarray, sizes: np.ndarray, profiles: sparse.csr_matrix,
    rng: np.random.Generator,
) -> np.ndarray:
    """Convert raw similarity to a z-score against size-matched random pairs."""
    bins = np.digitize(sizes, SIZE_BINS)
    sample = np.column_stack([
        rng.integers(0, profiles.shape[0], BACKGROUND_PAIRS),
        rng.integers(0, profiles.shape[0], BACKGROUND_PAIRS),
    ])
    sample = sample[sample[:, 0] != sample[:, 1]]
    sample_raw = cosine_similarity(profiles, sample)
    sample_key = bins[sample[:, 0]] * len(SIZE_BINS) + bins[sample[:, 1]]

    stats: dict[int, tuple[float, float]] = {}
    for key in np.unique(sample_key):
        values = sample_raw[sample_key == key]
        if len(values) >= 50:
            stats[int(key)] = (float(values.mean()), float(values.std()) or 1.0)

    overall = (float(sample_raw.mean()), float(sample_raw.std()) or 1.0)
    pair_key = bins[pairs[:, 0]] * len(SIZE_BINS) + bins[pairs[:, 1]]
    mean = np.array([stats.get(int(k), overall)[0] for k in pair_key])
    std = np.array([stats.get(int(k), overall)[1] for k in pair_key])
    return (raw - mean) / std


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curated", type=Path, required=True)
    parser.add_argument("--target-sim", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--min-targets", type=int, default=1)
    args = parser.parse_args()

    actives = pd.read_parquet(args.curated / "actives.parquet")
    inactives = pd.read_parquet(args.curated / "inactives.parquet")
    print(f"actives {len(actives):,d}  inactives {len(inactives):,d}")

    calls = aggregate(actives, inactives)
    print(f"activity calls: {len(calls):,d} over {calls.molregno.nunique():,d} compounds")

    depth = calls.groupby("molregno").size()
    keep = set(depth[depth >= args.min_targets].index)
    calls = calls[calls["molregno"].isin(keep)]
    print(f"after min-targets={args.min_targets}: {calls.molregno.nunique():,d} compounds")

    relatedness = sparse.load_npz(args.target_sim)
    tids = pd.read_parquet(args.target_sim.with_suffix(".tids.parquet"))["tid"].to_numpy()

    raw, profiles, molregnos = build_profiles(calls, tids, relatedness)
    print(f"raw profiles:      {raw.shape}, nnz={raw.nnz:,d}")
    print(f"smoothed profiles: {profiles.shape}, nnz={profiles.nnz:,d}")

    args.out.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(args.out / "profiles_raw.npz", raw)
    sparse.save_npz(args.out / "profiles.npz", profiles)
    pd.DataFrame({"molregno": molregnos}).to_parquet(
        args.out / "profile_molregnos.parquet", index=False
    )
    calls.to_parquet(args.out / "activity_calls.parquet", index=False)
    print(f"wrote teacher profiles -> {args.out}")


if __name__ == "__main__":
    main()
