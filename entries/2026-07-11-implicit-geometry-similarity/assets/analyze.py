"""Summarize benchmark results into a markdown table and a summary figure."""

from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(__file__)
RES = os.path.join(HERE, "results")

REPS = ["mw+rf", "ecfp4+rf", "rdkit2d+rf", "chemberta+rf", "molformer+rf",
        "pots+rf", "pots_aug+rf", "ecfp_pots+rf", "ecfp_tanimoto+gp", "pots+gp"]
NICE = {"mw+rf": "MolWt", "ecfp4+rf": "ECFP4", "rdkit2d+rf": "RDKit2D",
        "chemberta+rf": "ChemBERTa", "molformer+rf": "MoLFormer", "pots+rf": "POTS",
        "pots_aug+rf": "POTS+glob", "ecfp_pots+rf": "ECFP+POTS",
        "ecfp_tanimoto+gp": "Tani-GP", "pots+gp": "POTS-GP"}


def load(name):
    p = os.path.join(RES, f"{name}.json")
    return json.load(open(p)) if os.path.exists(p) else {}


def table(metric="spearman"):
    lines = []
    all_rows = []
    for suite in ["asap", "expansionrx", "pxr", "openbind"]:
        d = load(suite)
        if not d:
            continue
        lines.append(f"\n## {suite}  ({metric})\n")
        header = "| endpoint | " + " | ".join(NICE[r] for r in REPS) + " |"
        sep = "|" + "---|" * (len(REPS) + 1)
        lines += [header, sep]
        for ep, res in d.items():
            vals = [res.get(r, {}).get(metric, float("nan")) for r in REPS]
            all_rows.append(vals)
            # bold the best feature-RF model (all *+rf reps)
            rf_idx = [i for i, r in enumerate(REPS) if r.endswith("+rf")]
            best = max((v for i, v in enumerate(vals) if i in rf_idx and np.isfinite(v)),
                       default=float("nan"))
            cells = []
            for i, v in enumerate(vals):
                s = "—" if not np.isfinite(v) else f"{v:.3f}"
                if np.isfinite(v) and v == best and i in rf_idx:
                    s = f"**{s}**"
                cells.append(s)
            lines.append(f"| {ep} | " + " | ".join(cells) + " |")
    # aggregate means
    arr = np.array(all_rows, dtype=float)
    if arr.size:
        lines.append(f"\n## Mean across all endpoints ({metric})\n")
        lines.append("| " + " | ".join(NICE[r] for r in REPS) + " |")
        lines.append("|" + "---|" * len(REPS))
        means = np.nanmean(arr, axis=0)
        lines.append("| " + " | ".join(f"{m:.3f}" for m in means) + " |")
    return "\n".join(lines)


def win_summary(metric="spearman"):
    """How often each POTS variant beats ECFP4 across all endpoints."""
    lines = [f"\n## POTS vs ECFP4 head-to-head ({metric})\n"]
    rows = []
    for suite in ["asap", "expansionrx", "pxr", "openbind"]:
        for ep, res in load(suite).items():
            base = res.get("ecfp4+rf", {}).get(metric, float("nan"))
            rows.append((ep, base, {r: res.get(r, {}).get(metric, float("nan"))
                                    for r in REPS}))
    for rep in ["pots+rf", "pots_aug+rf", "ecfp_pots+rf", "pots+gp"]:
        wins = sum(1 for _, b, d in rows
                   if np.isfinite(d[rep]) and np.isfinite(b) and d[rep] > b)
        tot = sum(1 for _, b, d in rows
                  if np.isfinite(d[rep]) and np.isfinite(b))
        delta = np.nanmean([d[rep] - b for _, b, d in rows
                            if np.isfinite(d[rep]) and np.isfinite(b)])
        lines.append(f"- {NICE[rep]}: beats ECFP4 on {wins}/{tot} endpoints "
                     f"(mean Δ = {delta:+.3f})")
    return "\n".join(lines)


if __name__ == "__main__":
    parts = ["# POTS benchmark results"]
    for metric in ["spearman", "pearson", "mae", "roc_auc", "prec10", "rec10"]:
        parts.append(table(metric))
    parts.append(win_summary("spearman"))
    parts.append(win_summary("roc_auc"))
    md = "\n".join(parts).rstrip() + "\n"
    print(md)
    with open(os.path.join(RES, "table.md"), "w") as fh:
        fh.write(md)
