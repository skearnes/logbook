"""Figures for the bioactivity-fingerprint entry.

Two panels, both built around the same point: a global score hides whether a
representation has learned anything beyond structure, so the interesting
comparison is inside the regimes where structure and bioactivity disagree.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (backend must be set first)

from evaluate import CLIFF_MIN, HOP_MAX, RELATED_Z  # noqa: E402

# Colour-blind-safe; the control and the trained student are deliberately the
# two most distinguishable hues because their comparison is the headline.
COLOURS = {
    "ecfp4": "#4C72B0",
    "ecfp6": "#7A9CC6",
    "fcfp4": "#A8C0DC",
    "rdkit": "#C4C4C4",
    "atompair": "#B0B0B0",
    "rdkit2d": "#9C9C9C",
    "usrcat": "#8C8C8C",
    "neighbour": "#DD8452",
    "multitask": "#55A868",
    "multitask_profile": "#8FBF9F",
    "distill": "#2E7D5B",
}

REGIMES = [
    ("scaffold_hop", "Scaffold-hop regime\n(ECFP4 < 0.35)"),
    ("cliff", "Cliff regime\n(ECFP4 > 0.60)"),
    ("overall", "All pairs"),
]


def load(path: Path) -> dict:
    with (path / "results.json").open() as fh:
        return json.load(fh)


def panel(ax, results: dict, regime: str, metric: str, title: str) -> None:
    names = [n for n in COLOURS if n in results]
    values = [results[n][regime][metric] for n in names]
    order = np.argsort(values)
    names = [names[i] for i in order]
    values = [values[i] for i in order]

    bars = ax.barh(
        names, values, color=[COLOURS[n] for n in names], edgecolor="none"
    )
    # ECFP4 is the reference every other method has to beat, so mark it.
    if "ecfp4" in names:
        ax.axvline(
            results["ecfp4"][regime][metric], color="#4C72B0",
            linestyle="--", linewidth=1, zorder=0,
        )
    if metric == "auc":
        ax.axvline(0.5, color="#999999", linestyle=":", linewidth=1, zorder=0)

    for bar, value in zip(bars, values):
        if np.isfinite(value):
            ax.text(
                value + 0.005, bar.get_y() + bar.get_height() / 2,
                f"{value:.3f}", va="center", fontsize=8,
            )

    ax.set_title(title, fontsize=10)
    ax.set_xlabel(metric.upper() if metric == "auc" else "Spearman")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.25, linewidth=0.5)
    ax.set_axisbelow(True)


def plane(pairs: pd.DataFrame, out: Path) -> None:
    """The structure/bioactivity plane the evaluation is built around.

    The diagonal is where the two agree and where a structural fingerprint is
    already sufficient. The off-diagonal corners are the whole reason to want a
    bioactivity fingerprint, and they are where the methods are compared.
    """
    fig, ax = plt.subplots(figsize=(6.4, 5.2))

    sample = pairs.sample(min(60_000, len(pairs)), random_state=0)
    ax.scatter(
        sample["ecfp4"], sample["teacher_z"],
        s=2, alpha=0.06, color="#4C72B0", edgecolors="none", rasterized=True,
    )

    ax.axvline(HOP_MAX, color="#DD8452", linestyle="--", linewidth=1)
    ax.axvline(CLIFF_MIN, color="#C44E52", linestyle="--", linewidth=1)
    ax.axhline(RELATED_Z, color="#666666", linestyle=":", linewidth=1)

    ax.annotate(
        "scaffold hops\nteacher: similar\nECFP: distant",
        xy=(0.06, 0.80), xycoords="axes fraction", fontsize=9, color="#DD8452",
    )
    ax.annotate(
        "activity cliffs\nECFP: near-identical\nteacher: unrelated",
        xy=(0.62, 0.06), xycoords="axes fraction", fontsize=9, color="#C44E52",
    )

    ax.set_xlabel("ECFP4 Tanimoto (structural similarity)")
    ax.set_ylabel("Teacher z-score (bioactivity similarity)")
    ax.set_title("Where structure and bioactivity disagree", fontsize=11)
    ax.set_ylim(-5, 20)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "plane.png", dpi=200, bbox_inches="tight")
    print(f"wrote {out / 'plane.png'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, default=None)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    results = load(args.results)

    if args.pairs is not None:
        plane(pd.read_parquet(args.pairs / "pairs.parquet"), args.out)

    for metric in ("auc", "spearman"):
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
        for ax, (regime, title) in zip(axes, REGIMES):
            panel(ax, results, regime, metric, title)
        fig.suptitle(
            "Agreement with the assay-derived teacher, by structural regime",
            fontsize=11,
        )
        fig.tight_layout()
        fig.savefig(args.out / f"regimes_{metric}.png", dpi=200, bbox_inches="tight")
        print(f"wrote {args.out / f'regimes_{metric}.png'}")


if __name__ == "__main__":
    main()
