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

"""Place POTS/baselines in the ExpansionRx FINAL leaderboards (RAE).

Reads the per-endpoint RAE that score_leaderboard.py computed on the official
final test and reports, for each model, its rank in (a) the aggregate MA-RAE
final leaderboard and (b) each of the 9 per-endpoint final leaderboards. All 103
finalists were scored on the full blinded test, so the ranks are directly
comparable.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

DATA = os.path.join(os.path.dirname(__file__), "data")
LB = os.path.join(DATA, "leaderboard")

# ENDPOINT name (as in the test CSV / our RAE table) -> per-endpoint LB filename.
EP_FILE = {
    "LogD": "exprx_ep_LogD.csv",
    "KSOL": "exprx_ep_KSOL.csv",
    "HLM CLint": "exprx_ep_HLM_CLint.csv",
    "MLM CLint": "exprx_ep_MLM_CLint.csv",
    "Caco-2 Permeability Papp A>B": "exprx_ep_Caco-2_Permeability_Papp_AB.csv",
    "Caco-2 Permeability Efflux": "exprx_ep_Caco-2_Permeability_Efflux.csv",
    "MPPB": "exprx_ep_MPPB.csv",
    "MBPB": "exprx_ep_MBPB.csv",
    "MGMB": "exprx_ep_MGMB.csv",
}


def _nums(series):
    return series.astype(str).str.extract(r"([0-9.]+)").astype(float)[0].to_numpy()


def rank_of(value, board_scores):
    """1-based rank if `value` were inserted into `board_scores` (lower=better)."""
    return int((board_scores < value).sum()) + 1


def main():
    rae = pd.read_csv(os.path.join(LB, "pots_expansionrx_rae.csv"), index_col=0)
    # rae: index = endpoints, columns = models
    models = list(rae.columns)
    n = 103

    print("=== Per-endpoint FINAL leaderboard placement (rank / 103, by RAE) ===\n")
    header = f"{'endpoint':30s} " + " ".join(f"{m:>11s}" for m in models)
    print(header)
    per_ep_rank = {m: [] for m in models}
    for ept, fn in EP_FILE.items():
        board = _nums(pd.read_csv(os.path.join(LB, fn))["RAE_display"])
        cells = []
        for m in models:
            v = rae.loc[ept, m]
            r = rank_of(v, board)
            per_ep_rank[m].append(r)
            cells.append(f"{v:.3f}#{r}")
        print(f"{ept:30s} " + " ".join(f"{c:>11s}" for c in cells))
    print(f"\n(cell = RAE#rank; winner/median/worst per endpoint vary — "
          f"e.g. LogD winner {(_nums(pd.read_csv(os.path.join(LB, EP_FILE['LogD']))['RAE_display'])).min():.3f})")

    # Aggregate MA-RAE
    agg = _nums(pd.read_csv(os.path.join(LB, "exprx_leaderboard_cld_results.csv"))["RAE_display"])
    print(f"\n=== Aggregate MA-RAE FINAL leaderboard (winner {agg.min():.3f}, "
          f"median {np.median(agg):.3f}, worst {agg.max():.3f}, n={n}) ===\n")
    ma = {m: float(rae[m].mean()) for m in models}
    for m in models:
        r = rank_of(ma[m], agg)
        best_ep = min(per_ep_rank[m]); worst_ep = max(per_ep_rank[m])
        print(f"  {m:12s} MA-RAE {ma[m]:.4f} -> rank ~{r}/{n} (top {100*(r-1)/n:.0f}%)"
              f"   [per-endpoint ranks {best_ep}-{worst_ep}]")

    # Save a tidy table
    rows = []
    for m in models:
        for i, ept in enumerate(EP_FILE):
            rows.append({"model": m, "endpoint": ept, "RAE": rae.loc[ept, m],
                         "rank": per_ep_rank[m][i], "n_finalists": n})
        rows.append({"model": m, "endpoint": "MA-RAE (aggregate)", "RAE": ma[m],
                     "rank": rank_of(ma[m], agg), "n_finalists": n})
    pd.DataFrame(rows).to_csv(os.path.join(LB, "pots_placement_full.csv"), index=False)
    print("\nwrote pots_placement_full.csv")


if __name__ == "__main__":
    main()
