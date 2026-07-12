"""Place POTS on the ExpansionRx blind-challenge FINAL leaderboard.

Reproduces the challenge's exact scoring (from the leaderboard Space's
``evaluate.py`` / ``utils.py``): per-endpoint RAE = MAE / mean|y-mean(y)| on
log10(clip(y,0)+1)-transformed values (LogD left raw), bootstrapped 1000×
(seed 0), then macro-averaged across the 9 endpoints (MA-RAE). Models are
trained on the official challenge train split and predicted on the official
(unblinded) test split — the same 2282 compounds the final leaderboard used.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

import bench

DATA = os.path.join(os.path.dirname(__file__), "data")
N_WORKERS = max(1, (os.cpu_count() or 4) - 1)
ENDPOINTS = ["LogD", "KSOL", "HLM CLint", "MLM CLint",
             "Caco-2 Permeability Papp A>B", "Caco-2 Permeability Efflux",
             "MPPB", "MBPB", "MGMB"]


def clip_and_log(y):
    return np.log10(np.clip(y, 0, None) + 1)


def rae_bootstrap(pred, true, n=1000, seed=0):
    """Bootstrap-mean RAE = MAE / mean|true-mean(true)|, matching the challenge."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(true), size=(n, len(true)), replace=True)
    vals = []
    for ix in idx:
        t, p = true[ix], pred[ix]
        vals.append(np.mean(np.abs(t - p)) / np.mean(np.abs(t - np.mean(t))))
    return float(np.mean(vals))


def build_features(train_smiles, test_smiles):
    """Shared features for all endpoints: ECFP, RDKit2D, and one POTS landmark
    embedding (+global) using landmarks from the full training set."""
    smiles = list(train_smiles) + list(test_smiles)
    tr_idx = np.arange(len(train_smiles))
    ecfp = bench.ecfp4(smiles)
    rdk = bench.rdkit2d(smiles)
    pots_lm = bench.pots_landmark_embedding(smiles, tr_idx, n_landmarks=64,
                                            n_workers=N_WORKERS)
    pots_glob = bench.pots_global_descriptors(smiles)
    n = len(train_smiles)
    feats = {
        "ECFP4": ecfp,
        "RDKit2D": rdk,
        "MoLFormer": bench.molformer(smiles),
        "POTS": np.concatenate([pots_lm, pots_glob], 1),
        "ECFP+POTS": np.concatenate([ecfp, pots_lm, pots_glob], 1),
    }
    return feats, n


def main():
    tr = pd.read_csv(os.path.join(DATA, "expansionrx_official_train.csv"))
    te = pd.read_csv(os.path.join(DATA, "expansionrx_official_test.csv"))
    feats, n_tr = build_features(tr["SMILES"].tolist(), te["SMILES"].tolist())

    rae = {m: {} for m in feats}
    for ept in ENDPOINTS:
        raw = ["logd"]  # endpoints left untransformed
        transform = ept.lower() not in raw
        ytr = tr[ept].to_numpy(float)
        yte = te[ept].to_numpy(float)
        tr_mask = np.isfinite(ytr)
        te_mask = np.isfinite(yte)
        ytr_t = clip_and_log(ytr) if transform else ytr
        yte_t = clip_and_log(yte) if transform else yte
        for m, X in feats.items():
            Xtr, Xte = X[:n_tr], X[n_tr:]
            rf = RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=0)
            rf.fit(Xtr[tr_mask], ytr_t[tr_mask])
            pred = rf.predict(Xte[te_mask])
            rae[m][ept] = rae_bootstrap(pred, yte_t[te_mask])
        print(f"  {ept:32s} " + "  ".join(f"{m}={rae[m][ept]:.3f}" for m in feats),
              flush=True)

    print("\n=== MA-RAE (macro-average over 9 endpoints) ===")
    ma = {m: float(np.mean([rae[m][e] for e in ENDPOINTS])) for m in feats}
    for m in feats:
        print(f"  {m:12s} MA-RAE = {ma[m]:.4f}")

    # Place in the final leaderboard.
    lb = pd.read_csv(os.path.join(DATA, "leaderboard",
                                  "exprx_leaderboard_cld_results.csv"))
    lb["marae"] = lb["RAE_display"].astype(str).str.extract(r"([0-9.]+)").astype(float)
    scores = lb["marae"].to_numpy()
    n_fin = len(lb)
    print(f"\n=== Placement among {n_fin} finalists (winner {scores.min():.4f}, "
          f"median {np.median(scores):.4f}, last {scores.max():.4f}) ===")
    for m in feats:
        rank = int((scores < ma[m]).sum()) + 1  # 1 + #finalists strictly better
        pct = 100 * (rank - 1) / n_fin
        print(f"  {m:12s} MA-RAE {ma[m]:.4f} -> rank ~{rank}/{n_fin} (top {pct:.0f}%)")

    pd.DataFrame(rae).to_csv(os.path.join(DATA, "leaderboard",
                                          "pots_expansionrx_rae.csv"))
    import json
    json.dump({"ma_rae": ma, "n_finalists": int(n_fin),
               "winner": float(scores.min()), "median": float(np.median(scores))},
              open(os.path.join(DATA, "leaderboard", "pots_placement.json"), "w"),
              indent=2)


if __name__ == "__main__":
    main()
