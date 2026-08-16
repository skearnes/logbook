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

"""Re-run only the POTS representations with a single-ETKDG-conformer geometry.

The ablation (ablate.py) suggested a real conformer beats the bounds-matrix
geometry. Here we swap the POTS ground metric to ``geometry="conformer"`` across
all 18 endpoints and compare to the committed bounds results, reusing the exact
same splits. Only POTS is recomputed; ECFP/RDKit/FM/GP-Tanimoto numbers are read
from the committed results.
"""

from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor

import bench
from run_bench import (ASAP, EXPANSION, DATA, RESULTS, LOG_TARGETS, GP_MAX_TRAIN,
                       N_WORKERS, scaffold_split, _metrics)

GEOM = "conformer"


def _pots_reps(smiles, y, tr, te, log_target):
    y = np.asarray(y, float)
    if log_target:
        y = np.log10(np.clip(y, 1e-6, None))
    lm = bench.pots_landmark_embedding(
        smiles, tr, n_landmarks=64, geometry=GEOM, n_workers=N_WORKERS)
    glob = bench.pots_global_descriptors(smiles, geometry=GEOM)
    ecfp = bench.ecfp4(smiles)
    out = {}
    for name, X in [("pots", lm), ("pots_aug", np.concatenate([lm, glob], 1)),
                    ("ecfp_pots", np.concatenate([ecfp, lm, glob], 1))]:
        rf = RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=0)
        rf.fit(X[tr], y[tr])
        out[f"{name}+rf"] = _metrics(y[te], rf.predict(X[te]))
    if len(tr) <= GP_MAX_TRAIN:
        Dtr = bench.pots_distance_matrix(smiles, tr, tr, geometry=GEOM, n_workers=N_WORKERS)
        Dte = bench.pots_distance_matrix(smiles, te, tr, geometry=GEOM, n_workers=N_WORKERS)
        out["pots+gp"] = _metrics(y[te], bench.gp_regression(Dtr, Dte, y[tr]))
    return out


def _split_csv(csv, split_json, tcol):
    df = pd.read_csv(os.path.join(DATA, csv))
    tcol = tcol if tcol in df.columns else df.columns[-1]
    sp = json.load(open(os.path.join(DATA, split_json)))
    smiles = df["SMILES"].tolist(); y = df[tcol].to_numpy()
    yn = pd.to_numeric(df[tcol], errors="coerce").to_numpy()
    tr = np.array(sp["train"]); te = np.array(sp["test"])
    return smiles, y, tr[np.isfinite(yn[tr])], te[np.isfinite(yn[te])], tcol


def main():
    out = {}
    for name, csv, split, tcol in ASAP:
        s, y, tr, te, tc = _split_csv(csv, split, tcol)
        print(f"[ASAP {name}] n={len(s)} tr={len(tr)} te={len(te)}", flush=True)
        out[f"ASAP/{name}"] = _pots_reps(s, y, tr, te, tc in LOG_TARGETS)
        json.dump(out, open(os.path.join(RESULTS, "conformer.json"), "w"), indent=2)
    for name, stem, tcol in EXPANSION:
        s, y, tr, te, tc = _split_csv(f"exprx_{stem}_{stem}.csv",
                                      f"exprx_{stem}_time_split_{stem}.json", tcol)
        print(f"[ExpansionRx {name}] n={len(s)} tr={len(tr)} te={len(te)}", flush=True)
        out[f"ExpansionRx/{name}"] = _pots_reps(s, y, tr, te, tc in LOG_TARGETS)
        json.dump(out, open(os.path.join(RESULTS, "conformer.json"), "w"), indent=2)
    # PXR
    trd = pd.read_csv(os.path.join(DATA, "pxr_train.csv"))[["SMILES", "pEC50"]].dropna()
    ted = pd.read_csv(os.path.join(DATA, "pxr_test1.csv"))[["SMILES", "pEC50"]].dropna()
    s = trd["SMILES"].tolist() + ted["SMILES"].tolist()
    y = np.concatenate([trd["pEC50"].to_numpy(), ted["pEC50"].to_numpy()])
    tr = np.arange(len(trd)); te = np.arange(len(trd), len(trd) + len(ted))
    print(f"[PXR] n={len(s)} tr={len(tr)} te={len(te)}", flush=True)
    out["PXR/pEC50"] = _pots_reps(s, y, tr, te, False)
    json.dump(out, open(os.path.join(RESULTS, "conformer.json"), "w"), indent=2)
    # OpenBind
    df = pd.read_csv(os.path.join(DATA, "openbind_ev71.csv"))
    s = df["SMILES"].tolist(); y = df["pKD"].to_numpy()
    tr, te = scaffold_split(s, 0.7, 0)
    print(f"[OpenBind] n={len(s)} tr={len(tr)} te={len(te)}", flush=True)
    out["OpenBind/pKD"] = _pots_reps(s, y, tr, te, False)
    json.dump(out, open(os.path.join(RESULTS, "conformer.json"), "w"), indent=2)
    print("CONFORMER ABLATION DONE", flush=True)


if __name__ == "__main__":
    main()
