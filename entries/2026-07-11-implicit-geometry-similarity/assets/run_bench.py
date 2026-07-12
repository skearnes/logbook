"""Driver: evaluate representations on ADMET/potency tasks under time splits.

Usage:
    python run_bench.py asap          # 7 ASAP endpoints
    python run_bench.py expansionrx   # 9 ExpansionRx endpoints (scaffold split)
    python run_bench.py pxr           # PXR induction (train -> unblinded test)
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr, pearsonr
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import roc_auc_score

import bench

DATA = os.path.join(os.path.dirname(__file__), "data")
RESULTS = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(RESULTS, exist_ok=True)
LOG_TARGETS = {"HLM", "MLM", "KSOL", "HLM CLint", "MLM CLint", "MDR1-MDCKII",
               "Caco-2 Permeability Papp A>B", "Caco-2 Permeability Efflux"}
N_WORKERS = max(1, (os.cpu_count() or 4) - 1)
# Max training size for the exact FGW-kernel GP (bounded by the O(n^2) FGW cost,
# not GP scaling). Covers the ADMET + potency + affinity endpoints; larger sets
# rely on the landmark-RF features, which scale linearly.
GP_MAX_TRAIN = 1200


def scaffold_split(smiles: list[str], frac_train: float = 0.7, seed: int = 0):
    """Bemis-Murcko scaffold split: whole scaffolds go to train or test."""
    scaffolds: dict[str, list[int]] = {}
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        core = ""
        if m is not None:
            try:
                core = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
            except Exception:
                core = ""
        scaffolds.setdefault(core, []).append(i)
    groups = sorted(scaffolds.values(), key=len, reverse=True)
    n_train_target = int(frac_train * len(smiles))
    train, test = [], []
    for g in groups:
        if len(train) + len(g) <= n_train_target:
            train += g
        else:
            test += g
    return np.array(sorted(train)), np.array(sorted(test))


def _prec_recall_at_k(active, score, frac=0.1):
    """Precision and recall at the top ``frac`` of ranked predictions."""
    n = len(score)
    k = max(1, int(round(frac * n)))
    top = active[np.argsort(-score)[:k]]  # high score = predicted active
    n_active = active.sum()
    prec = float(top.sum() / k)
    rec = float(top.sum() / n_active) if n_active else float("nan")
    return prec, rec


def _metrics(y_true, y_pred):
    """Regression (Spearman, Pearson, MAE, RMSE) + thresholded ranking
    (ROC-AUC, precision@10%, recall@10%). Actives = top quartile of truth."""
    ok = np.isfinite(y_true) & np.isfinite(y_pred)
    yt, yp = y_true[ok], y_pred[ok]
    out = dict(spearman=float("nan"), pearson=float("nan"), mae=float("nan"),
               rmse=float("nan"), roc_auc=float("nan"), prec10=float("nan"),
               rec10=float("nan"), n=int(len(yt)))
    if len(yt) < 5 or np.std(yp) == 0:
        return out
    out["spearman"] = float(spearmanr(yt, yp).statistic)
    out["pearson"] = float(pearsonr(yt, yp)[0])
    out["mae"] = float(np.mean(np.abs(yt - yp)))
    out["rmse"] = float(np.sqrt(np.mean((yt - yp) ** 2)))
    # Threshold at the top quartile to obtain a discrimination task.
    thr = np.quantile(yt, 0.75)
    active = (yt >= thr).astype(int)
    if 0 < active.sum() < len(active):
        out["roc_auc"] = float(roc_auc_score(active, yp))
        out["prec10"], out["rec10"] = _prec_recall_at_k(active, yp, 0.1)
    return out


def evaluate_task(name, smiles, y, train_idx, test_idx, log_target=False,
                  n_landmarks=64, pred_key=None, names=None):
    """Evaluate all representations on one regression task.

    If ``pred_key`` is given, per-rep test-set predictions are saved to a wide
    parquet ``results/predictions/<pred_key>.parquet`` -- one row per test
    molecule (name, SMILES, y_true) with one column per method -- for separate
    analysis.
    """
    y = np.asarray(y, dtype=float)
    if log_target:
        y = np.log10(np.clip(y, 1e-6, None))
    res = {}
    tr, te = train_idx, test_idx
    preds = {}  # method name -> test predictions

    # Compute the (expensive) POTS landmark embedding once; derive variants.
    pots_lm = bench.pots_landmark_embedding(
        smiles, tr, n_landmarks=n_landmarks, n_workers=N_WORKERS)
    pots_glob = bench.pots_global_descriptors(smiles)

    # --- Feature-based reps, all through the same RandomForest ------------- #
    feats = {
        "mw": lambda: bench.molweight(smiles),
        "ecfp4": lambda: bench.ecfp4(smiles),
        "rdkit2d": lambda: bench.rdkit2d(smiles),
        "chemberta": lambda: bench.chemberta(smiles),
        "molformer": lambda: bench.molformer(smiles),
        "pots": lambda: pots_lm,
        "pots_aug": lambda: np.concatenate([pots_lm, pots_glob], axis=1),
        "ecfp_pots": lambda: np.concatenate(
            [bench.ecfp4(smiles), pots_lm, pots_glob], axis=1),
    }
    for rep, fn in feats.items():
        t = time.time()
        X = fn()
        rf = RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=0)
        rf.fit(X[tr], y[tr])
        pred = rf.predict(X[te])
        m = _metrics(y[te], pred)
        m["seconds"] = round(time.time() - t, 1)
        res[f"{rep}+rf"] = m
        preds[f"{rep}+rf"] = pred
        print(f"    {rep+'+rf':16s} rho={m['spearman']:.3f} "
              f"r={m['pearson']:.3f} rmse={m['rmse']:.3f} ({m['seconds']}s)")

    # --- Similarity-as-kernel comparison via Gaussian process -------------- #
    # The metric IS the GP covariance k(A,B)=exp(-gamma*D): Tanimoto (PSD) vs
    # FGW (PSD-repaired). Uses the exact full-rank GP; a Nyström approximation
    # was tested and found to distort the indefinite FGW kernel (see ablations),
    # so we prefer exact kernels wherever the O(n^2) FGW matrix is affordable.
    tr_smi = [smiles[i] for i in tr]
    te_smi = [smiles[i] for i in te]
    Dtt = bench.tanimoto_distance_matrix(tr_smi, tr_smi)
    Dts = bench.tanimoto_distance_matrix(te_smi, tr_smi)
    pred = bench.gp_regression(Dtt, Dts, y[tr])
    res["ecfp_tanimoto+gp"] = _metrics(y[te], pred)
    preds["ecfp_tanimoto+gp"] = pred
    print(f"    {'ecfp-tani+gp':16s} rho={res['ecfp_tanimoto+gp']['spearman']:.3f}")

    if len(tr) <= GP_MAX_TRAIN:  # exact FGW-GP where the O(n^2) matrix is affordable
        Dtr = bench.pots_distance_matrix(smiles, tr, tr, n_workers=N_WORKERS)
        Dte = bench.pots_distance_matrix(smiles, te, tr, n_workers=N_WORKERS)
        pred = bench.gp_regression(Dtr, Dte, y[tr])
        res["pots+gp"] = _metrics(y[te], pred)
        preds["pots+gp"] = pred
        print(f"    {'pots+gp':16s} rho={res['pots+gp']['spearman']:.3f}")

    if pred_key is not None:
        pdir = os.path.join(RESULTS, "predictions")
        os.makedirs(pdir, exist_ok=True)
        df = pd.DataFrame({
            "name": [names[i] for i in te] if names is not None
                    else [smiles[i] for i in te],
            "smiles": [smiles[i] for i in te],
            "y_true": y[te],
            "log_target": log_target,
        })
        for method, p in preds.items():
            df[method] = p
        df.to_parquet(os.path.join(pdir, f"{pred_key}.parquet"), index=False)
    return res


# --------------------------------------------------------------------------- #
# Dataset loaders
# --------------------------------------------------------------------------- #
ASAP = [
    ("HLM", "asap_HLM_hlm.csv", "asap_HLM_time_split_hlm.json", "HLM"),
    ("KSOL", "asap_KSOL_ksol.csv", "asap_KSOL_time_split_ksol.json", "KSOL"),
    ("LogD", "asap_LogD_logd.csv", "asap_LogD_time_split_logd.json", "LogD"),
    ("MDR1-MDCKII", "asap_MDR1-MDCKII_mdr1_mdckii.csv",
     "asap_MDR1-MDCKII_time_split_mdr1_mdckii.json", "MDR1-MDCKII"),
    ("MLM", "asap_MLM_mlm.csv", "asap_MLM_time_split_mlm.json", "MLM"),
    ("pIC50_SARS", "asap_pIC50_SARS-CoV-2_pic50_sars_cov_2.csv",
     "asap_pIC50_SARS-CoV-2_time_split_pic50_sars_cov_2.json", "pIC50"),
    ("pIC50_MERS", "asap_pIC50_MERS-CoV_pic50_mers_cov.csv",
     "asap_pIC50_MERS-CoV_time_split_pic50_mers_cov.json", "pIC50"),
]

# ExpansionRx: scikit-fingerprints per-endpoint mirrors, each with a curated
# time split. (name, csv_stem, target_col).
EXPANSION = [
    ("LogD", "logd", "LogD"),
    ("KSOL", "ksol", "KSOL"),
    ("HLM_CLint", "hlm_clint", "HLM CLint"),
    ("MLM_CLint", "mlm_clint", "MLM CLint"),
    ("Caco2_Efflux", "caco_2_permeability_efflux", "Caco-2 Permeability Efflux"),
    ("Caco2_Papp_AB", "caco2_papp", "Caco-2 Permeability Papp A>B"),
    ("MPPB", "mppb", "MPPB"),
    ("MBPB", "mbpb", "MBPB"),
    ("MGMB", "mgmb", "MGMB"),
]


def _load_results(results_file):
    """Load existing results for resume (skip already-computed endpoints)."""
    path = os.path.join(RESULTS, results_file)
    return json.load(open(path)) if os.path.exists(path) else {}


def _run_split_csv(suite, name, csv, split_json, target_col, out, results_file,
                   log_target):
    """Evaluate one CSV that ships an index-based train/test split JSON."""
    if name in out and "ecfp4+rf" in out[name]:
        print(f"[{suite} {name}] already done, skipping")
        return
    df = pd.read_csv(os.path.join(DATA, csv))
    tcol = target_col if target_col in df.columns else df.columns[-1]
    sp = json.load(open(os.path.join(DATA, split_json)))
    smiles = df["SMILES"].tolist()
    y = df[tcol].to_numpy()
    name_col = next((c for c in ("Molecule name", "Molecule Name") if c in df.columns),
                    None)
    names = df[name_col].tolist() if name_col else None
    tr, te = np.array(sp["train"]), np.array(sp["test"])
    # Drop split indices whose target is missing.
    tr = tr[np.isfinite(pd.to_numeric(df[tcol], errors="coerce").to_numpy()[tr])]
    te = te[np.isfinite(pd.to_numeric(df[tcol], errors="coerce").to_numpy()[te])]
    print(f"\n[{suite} {name}] n={len(smiles)} train={len(tr)} test={len(te)}")
    out[name] = evaluate_task(name, smiles, y, tr, te, log_target=log_target,
                              pred_key=f"{suite}_{name}", names=names)
    json.dump(out, open(os.path.join(RESULTS, results_file), "w"), indent=2)


def run_asap():
    out = _load_results("asap.json")
    for name, csv, split, target_col in ASAP:
        _run_split_csv("ASAP", name, csv, split, target_col, out, "asap.json",
                       target_col in LOG_TARGETS)
    return out


def run_expansionrx():
    out = _load_results("expansionrx.json")
    for name, stem, target_col in EXPANSION:
        _run_split_csv("ExpansionRx", name, f"exprx_{stem}_{stem}.csv",
                       f"exprx_{stem}_time_split_{stem}.json", target_col, out,
                       "expansionrx.json", target_col in LOG_TARGETS)
    return out


def run_openbind():
    df = pd.read_csv(os.path.join(DATA, "openbind_ev71.csv"))
    smiles = df["SMILES"].tolist()
    y = df["pKD"].to_numpy()
    tr, te = scaffold_split(smiles, frac_train=0.7, seed=0)
    print(f"\n[OpenBind EV-A71 pKD] n={len(smiles)} train={len(tr)} test={len(te)}")
    out = {"pKD_EVA71": evaluate_task("OpenBind", smiles, y, tr, te,
                                      pred_key="OpenBind_pKD_EVA71")}
    json.dump(out, open(os.path.join(RESULTS, "openbind.json"), "w"), indent=2)
    return out


def run_pxr():
    tr_df = pd.read_csv(os.path.join(DATA, "pxr_train.csv"))
    te_df = pd.read_csv(os.path.join(DATA, "pxr_test1.csv"))
    tr_df = tr_df[["SMILES", "pEC50"]].dropna()
    te_df = te_df[["SMILES", "pEC50"]].dropna()
    smiles = tr_df["SMILES"].tolist() + te_df["SMILES"].tolist()
    y = np.concatenate([tr_df["pEC50"].to_numpy(), te_df["pEC50"].to_numpy()])
    tr = np.arange(len(tr_df))
    te = np.arange(len(tr_df), len(tr_df) + len(te_df))
    names = tr_df.get("Molecule Name", pd.Series(tr_df["SMILES"])).tolist() + \
        te_df.get("Molecule Name", pd.Series(te_df["SMILES"])).tolist()
    print(f"\n[PXR pEC50] n={len(smiles)} train={len(tr)} test={len(te)}")
    out = {"pEC50": evaluate_task("PXR", smiles, y, tr, te, log_target=False,
                                  pred_key="PXR_pEC50", names=names)}
    json.dump(out, open(os.path.join(RESULTS, "pxr.json"), "w"), indent=2)
    return out


RUNNERS = {"asap": run_asap, "expansionrx": run_expansionrx,
           "pxr": run_pxr, "openbind": run_openbind}

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "asap"
    order = list(RUNNERS) if which == "all" else [which]
    for suite in order:
        print(f"\n{'='*60}\n== {suite.upper()}\n{'='*60}", flush=True)
        RUNNERS[suite]()
