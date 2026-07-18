"""Transfer test: does the distilled fingerprint help on unrelated endpoints?

The intrinsic evaluation asks whether the student reproduces the teacher, but
both live inside ChEMBL, so a good score there partly reflects reproducing
ChEMBL's assay structure. This script is the external check: freeze the student
trained on ChEMBL, use its embedding as features on the ADMET/potency endpoints
from the 2026-07-11 entry, and compare against the baselines recorded there
under the same non-random splits and the same downstream model.

Those endpoints are unrelated to the ChEMBL targets the student was trained on,
so any gain is transfer rather than memorisation. Reference numbers from that
entry (mean Spearman across 18 endpoints): RDKit2D 0.577, Tanimoto-GP 0.591,
ECFP4 0.560, POTS 0.455.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor

import fingerprints
import student

POTS_ENTRY = Path(__file__).resolve().parents[2] / "2026-07-11-implicit-geometry-similarity"
POTS_DATA = POTS_ENTRY / "assets" / "data"

# Endpoints measured on a log scale in the source benchmark.
LOG_TARGETS = {
    "HLM", "MLM", "KSOL", "HLM CLint", "MLM CLint", "MDR1-MDCKII",
    "Caco-2 Permeability Papp A>B", "Caco-2 Permeability Efflux",
}

RF_KWARGS = dict(n_estimators=500, min_samples_leaf=2, n_jobs=-1, random_state=0)


def load_endpoint(csv: Path, split_json: Path) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Load one endpoint plus the benchmark's own train/test row indices."""
    frame = pd.read_csv(csv)
    with split_json.open() as fh:
        split = json.load(fh)
    train = np.asarray(split["train"], dtype=int)
    test = np.asarray(split["test"], dtype=int)
    return frame, train, test


def discover() -> list[tuple[str, Path, Path]]:
    """Pair each endpoint CSV with its split file."""
    found = []
    for split_path in sorted(POTS_DATA.glob("*time_split*.json")):
        stem = split_path.name.replace("_time_split", "")
        csv = POTS_DATA / stem.replace(".json", ".csv")
        if csv.exists():
            found.append((csv.stem, csv, split_path))
    return found


def target_column(frame: pd.DataFrame) -> str | None:
    """Pick the assay-readout column, ignoring identifiers and structures."""
    ignore = {"smiles", "canonical_smiles", "SMILES", "CXSMILES", "id", "ID"}
    numeric = [
        c for c in frame.columns
        if c not in ignore and pd.api.types.is_numeric_dtype(frame[c])
    ]
    return numeric[0] if numeric else None


def smiles_column(frame: pd.DataFrame) -> str | None:
    for name in ("smiles", "SMILES", "canonical_smiles", "CXSMILES"):
        if name in frame.columns:
            return name
    return None


def score(features: np.ndarray, values: np.ndarray, train: np.ndarray, test: np.ndarray) -> float:
    model = RandomForestRegressor(**RF_KWARGS)
    model.fit(features[train], values[train])
    predicted = model.predict(features[test])
    return float(spearmanr(predicted, values[test]).statistic)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    state = torch.load(args.checkpoint, map_location=student.DEVICE, weights_only=False)
    model = student.Encoder(state["n_input"], state["n_targets"]).to(student.DEVICE)
    model.load_state_dict(state["state_dict"])

    rows = []
    for name, csv, split_json in discover():
        frame, train, test = load_endpoint(csv, split_json)
        smiles_col = smiles_column(frame)
        target_col = target_column(frame)
        if smiles_col is None or target_col is None:
            print(f"skip {name}: no smiles/target column")
            continue

        smiles = frame[smiles_col].astype(str).tolist()
        values = frame[target_col].to_numpy(dtype=float)
        keep = np.isfinite(values)
        if target_col in LOG_TARGETS or name.split("_")[-1] in LOG_TARGETS:
            positive = keep & (values > 0)
            values = np.where(positive, np.log10(np.maximum(values, 1e-12)), values)

        train = train[keep[train]]
        test = test[keep[test]]
        if len(train) < 30 or len(test) < 20:
            print(f"skip {name}: too few rows ({len(train)}/{len(test)})")
            continue

        counts = fingerprints.count_fingerprints(smiles, "ecfp4")
        bits = fingerprints.bit_fingerprints(smiles, "ecfp4").astype(np.float32)
        embedding = student.embed(model, counts)

        result = {
            "endpoint": name,
            "n_train": len(train),
            "n_test": len(test),
            "ecfp4": score(bits, values, train, test),
            "student": score(embedding, values, train, test),
            # Concatenation is the fair test of whether the learned fingerprint
            # adds anything ECFP does not already carry.
            "ecfp4+student": score(
                np.hstack([bits, embedding]), values, train, test
            ),
        }
        rows.append(result)
        print(
            f"{name:42s} ecfp4 {result['ecfp4']:.3f}  "
            f"student {result['student']:.3f}  "
            f"both {result['ecfp4+student']:.3f}",
            flush=True,
        )

    table = pd.DataFrame(rows)
    args.out.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out / "transfer.csv", index=False)
    print("\nmean Spearman:")
    print(table[["ecfp4", "student", "ecfp4+student"]].mean().round(3).to_string())


if __name__ == "__main__":
    main()
