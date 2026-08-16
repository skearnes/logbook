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

"""Place POTS on the OpenADMET PXR-induction challenge (activity track).

Primary metric: RAE (Relative Absolute Error) on pEC50, bootstrapped 1000×
(same OpenADMET methodology as ExpansionRx; pEC50 is left untransformed, being
already a log quantity). Trains on the challenge train split and predicts on the
full unblinded test (phase 1 + phase 2 = the final-leaderboard set).

The real final leaderboard is served from a private S3 bucket, so it cannot be
downloaded; the only public anchor is that RAE ≈ 0.586 placed ~40th of 211
(top ~19%). POTS's RAE is reported against that anchor.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

import bench

DATA = os.path.join(os.path.dirname(__file__), "data")
N_WORKERS = max(1, (os.cpu_count() or 4) - 1)


def rae_bootstrap(pred, true, n=1000, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(true), size=(n, len(true)), replace=True)
    vals = [np.mean(np.abs(true[ix] - pred[ix])) /
            np.mean(np.abs(true[ix] - np.mean(true[ix]))) for ix in idx]
    return float(np.mean(vals))


def main():
    tr = pd.read_csv(os.path.join(DATA, "pxr_train.csv"))[["SMILES", "pEC50"]].dropna()
    te = pd.concat([
        pd.read_csv(os.path.join(DATA, "pxr_test1.csv"))[["SMILES", "pEC50"]],
        pd.read_csv(os.path.join(DATA, "pxr_test2.csv"))[["SMILES", "pEC50"]],
    ]).dropna().drop_duplicates("SMILES")
    smiles = tr["SMILES"].tolist() + te["SMILES"].tolist()
    n_tr = len(tr)
    ytr = tr["pEC50"].to_numpy(float)
    yte = te["pEC50"].to_numpy(float)
    print(f"PXR: train={n_tr} test={len(te)} (phase1+phase2)", flush=True)

    ecfp = bench.ecfp4(smiles)
    rdk = bench.rdkit2d(smiles)
    pots_lm = bench.pots_landmark_embedding(smiles, np.arange(n_tr),
                                            n_landmarks=64, n_workers=N_WORKERS)
    glob = bench.pots_global_descriptors(smiles)
    feats = {
        "ECFP4": ecfp, "RDKit2D": rdk, "MoLFormer": bench.molformer(smiles),
        "POTS": np.concatenate([pots_lm, glob], 1),
        "ECFP+POTS": np.concatenate([ecfp, pots_lm, glob], 1),
    }
    print("\n=== PXR RAE (lower is better; anchor: RAE 0.586 ~ rank 40/211) ===")
    out = {}
    for m, X in feats.items():
        rf = RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=0)
        rf.fit(X[:n_tr], ytr)
        pred = rf.predict(X[n_tr:])
        out[m] = rae_bootstrap(pred, yte)
        print(f"  {m:12s} RAE = {out[m]:.4f}")
    import json
    json.dump(out, open(os.path.join(DATA, "leaderboard", "pots_pxr_rae.json"), "w"),
              indent=2)


if __name__ == "__main__":
    main()
