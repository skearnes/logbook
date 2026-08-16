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

"""Ablations for POTS design choices, on small ASAP endpoints (POTS-GP).

- geometry: bounds (implicit 3D) vs topology (bond counts) vs conformer (ETKDG).
  Tests whether conformer-free distance-geometry actually helps over pure graph
  topology, and how it compares to a real single conformer.
- alpha: color-only (0) -> shape-only (1). Tests the fused trade-off.

Each cell reports test Spearman under the provided time split, using the FGW
kernel GP (the cleanest readout of the metric itself).
"""

from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

import bench
from run_bench import _metrics

DATA = os.path.join(os.path.dirname(__file__), "data")
RESULTS = os.path.join(os.path.dirname(__file__), "results")
N_WORKERS = max(1, (os.cpu_count() or 4) - 1)

# Small ASAP endpoints where the full-rank FGW-GP is affordable.
TASKS = [
    ("HLM", "asap_HLM_hlm.csv", "asap_HLM_time_split_hlm.json", "HLM", True),
    ("LogD", "asap_LogD_logd.csv", "asap_LogD_time_split_logd.json", "LogD", False),
    ("pIC50_SARS", "asap_pIC50_SARS-CoV-2_pic50_sars_cov_2.csv",
     "asap_pIC50_SARS-CoV-2_time_split_pic50_sars_cov_2.json", "pIC50", False),
]


def ablate_task(name, csv, split, tcol, log_target):
    df = pd.read_csv(os.path.join(DATA, csv))
    tcol = tcol if tcol in df.columns else df.columns[-1]
    sp = json.load(open(os.path.join(DATA, split)))
    smiles = df["SMILES"].tolist()
    y = np.asarray(df[tcol], dtype=float)
    if log_target:
        y = np.log10(np.clip(y, 1e-6, None))
    tr, te = np.array(sp["train"]), np.array(sp["test"])
    out = {"geometry": {}, "alpha": {}}
    # geometry sweep at alpha=0.5
    for geom in ["topology", "bounds", "conformer"]:
        Dtr = bench.pots_distance_matrix(smiles, tr, tr, geometry=geom,
                                         n_workers=N_WORKERS)
        Dte = bench.pots_distance_matrix(smiles, te, tr, geometry=geom,
                                         n_workers=N_WORKERS)
        pred = bench.gp_regression(Dtr, Dte, y[tr])
        out["geometry"][geom] = round(_metrics(y[te], pred)["spearman"], 3)
        print(f"  [{name}] geometry={geom:9s} rho={out['geometry'][geom]}")
    # alpha sweep at geometry=bounds (recompute distances per alpha via kwargs)
    for alpha in [0.0, 0.25, 0.5, 0.75, 1.0]:
        Dtr = bench.pots_distance_matrix(smiles, tr, tr, geometry="bounds",
                                         n_workers=N_WORKERS, alpha=alpha)
        Dte = bench.pots_distance_matrix(smiles, te, tr, geometry="bounds",
                                         n_workers=N_WORKERS, alpha=alpha)
        pred = bench.gp_regression(Dtr, Dte, y[tr])
        out["alpha"][str(alpha)] = round(_metrics(y[te], pred)["spearman"], 3)
        print(f"  [{name}] alpha={alpha:4.2f}     rho={out['alpha'][str(alpha)]}")
    return out


if __name__ == "__main__":
    results = {}
    for name, csv, split, tcol, logt in TASKS:
        print(f"\n=== {name} ===")
        results[name] = ablate_task(name, csv, split, tcol, logt)
        json.dump(results, open(os.path.join(RESULTS, "ablation.json"), "w"),
                  indent=2)
