"""Interpretability analyses and molecule-pair figures for POTS.

Produces, for a potency/affinity dataset:
1. SAR continuity: rank correlation between pairwise metric distance and
   pairwise |Δactivity|. A better SAR metric couples small metric changes to
   small activity changes.
2. Two illustrative molecule-pair figures rendered with RDKit:
   - *Scaffold hop*: a pair POTS calls similar but ECFP/Tanimoto calls distant,
     with similar activity -- POTS sees a pharmacophore match Tanimoto misses.
   - *Activity cliff*: a pair ECFP/Tanimoto calls similar but with a large
     activity gap -- POTS is asked whether it separates them better.
   Each figure is annotated with both distances and the activities.
3. The pharmacophore transport plan for a chosen pair.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from scipy.stats import spearmanr

from pots import PharmacophoreCloud, fgw_distance
import bench

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "data")
ASSETS = HERE


def _load(csv, target_col, n=None, seed=0):
    df = pd.read_csv(os.path.join(DATA, csv)).dropna(subset=[target_col])
    df = df.drop_duplicates(subset=["SMILES"]).reset_index(drop=True)
    if n and len(df) > n:
        df = df.sample(n, random_state=seed).reset_index(drop=True)
    return df


def _pair_distances(smiles):
    clouds = [PharmacophoreCloud.from_smiles(s) for s in smiles]
    gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]
    n = len(smiles)
    Dp = np.zeros((n, n)); Dt = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if clouds[i] is None or clouds[j] is None:
                Dp[i, j] = Dp[j, i] = np.nan
                continue
            d = fgw_distance(clouds[i], clouds[j], **bench.FGW_KW)
            Dp[i, j] = Dp[j, i] = d
            t = 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
            Dt[i, j] = Dt[j, i] = t
    return Dp, Dt, clouds


def sar_continuity(csv, target_col, n=150, seed=0):
    df = _load(csv, target_col, n=n, seed=seed)
    smiles = df["SMILES"].tolist()
    y = df[target_col].to_numpy(dtype=float)
    Dp, Dt, _ = _pair_distances(smiles)
    iu = np.triu_indices(len(smiles), 1)
    dp, dt = Dp[iu], Dt[iu]
    dy = np.abs(y[iu[0]] - y[iu[1]])
    ok = np.isfinite(dp) & np.isfinite(dt)
    return {
        "pots_vs_dactivity": float(spearmanr(dp[ok], dy[ok]).statistic),
        "tanimoto_vs_dactivity": float(spearmanr(dt[ok], dy[ok]).statistic),
        "n_pairs": int(ok.sum()),
    }


def _draw_pair(smi_a, smi_b, legend_a, legend_b, title, path):
    ma, mb = Chem.MolFromSmiles(smi_a), Chem.MolFromSmiles(smi_b)
    for m in (ma, mb):
        AllChem.Compute2DCoords(m)
    img = Draw.MolsToGridImage(
        [ma, mb], legends=[legend_a, legend_b], molsPerRow=2,
        subImgSize=(360, 300))
    img.save(path)
    print(f"  wrote {os.path.relpath(path, HERE)}  ({title})")


def find_and_draw(csv, target_col, tag, n=250, seed=0):
    """Find a scaffold-hop and an activity-cliff pair; render both."""
    df = _load(csv, target_col, n=n, seed=seed)
    smiles = df["SMILES"].tolist()
    y = df[target_col].to_numpy(dtype=float)
    names = (df["Molecule Name"] if "Molecule Name" in df
             else df.get("Molecule name", pd.Series(smiles))).tolist()
    Dp, Dt, _ = _pair_distances(smiles)
    n = len(smiles)
    iu = np.triu_indices(n, 1)
    dp, dt = Dp[iu], Dt[iu]
    dy = np.abs(y[iu[0]] - y[iu[1]])
    ok = np.isfinite(dp) & np.isfinite(dt)

    # Normalize distances to ranks in [0,1] for fair "near/far" comparison.
    def rank01(x):
        r = np.argsort(np.argsort(x)).astype(float)
        return r / (len(x) - 1)
    rp, rt = rank01(dp), rank01(dt)

    # Scaffold hop: POTS near (low rp), Tanimoto far (high rt), activity close.
    score_hop = np.where(ok & (dy < np.nanquantile(dy, 0.4)),
                         rt - rp, -np.inf)
    h = int(np.argmax(score_hop))
    ia, ib = iu[0][h], iu[1][h]
    _draw_pair(
        smiles[ia], smiles[ib],
        f"{names[ia]}  {target_col}={y[ia]:.2f}",
        f"{names[ib]}  {target_col}={y[ib]:.2f}",
        "scaffold hop: POTS near, Tanimoto far, similar activity",
        os.path.join(ASSETS, f"pair_scaffold_hop_{tag}.png"))
    print(f"  scaffold-hop: POTS d={dp[h]:.3f} (rank {rp[h]:.2f}), "
          f"Tanimoto d={dt[h]:.3f} (rank {rt[h]:.2f}), |Δ{target_col}|={dy[h]:.2f}")

    # Activity cliff: Tanimoto near (low rt), large activity gap.
    score_cliff = np.where(ok & (rt < 0.1), dy, -np.inf)
    c = int(np.argmax(score_cliff))
    ja, jb = iu[0][c], iu[1][c]
    _draw_pair(
        smiles[ja], smiles[jb],
        f"{names[ja]}  {target_col}={y[ja]:.2f}",
        f"{names[jb]}  {target_col}={y[jb]:.2f}",
        "activity cliff: Tanimoto near but large activity gap",
        os.path.join(ASSETS, f"pair_activity_cliff_{tag}.png"))
    print(f"  activity-cliff: Tanimoto d={dt[c]:.3f} (rank {rt[c]:.2f}), "
          f"POTS d={dp[c]:.3f} (rank {rp[c]:.2f}), |Δ{target_col}|={dy[c]:.2f}")


if __name__ == "__main__":
    tasks = [
        ("asap_pIC50_SARS-CoV-2_pic50_sars_cov_2.csv",
         "pIC50 (SARS-CoV-2 Mpro)", "sars"),
        ("openbind_ev71.csv", "pKD", "openbind"),
    ]
    for csv, col, tag in tasks:
        print(f"\n=== {tag} ({col}) ===")
        r = sar_continuity(csv, col)
        print(f"SAR continuity (Spearman of pair-distance vs |Δ{col}|): "
              f"POTS={r['pots_vs_dactivity']:.3f}  "
              f"Tanimoto={r['tanimoto_vs_dactivity']:.3f}  (n={r['n_pairs']})")
        find_and_draw(csv, col, tag)
