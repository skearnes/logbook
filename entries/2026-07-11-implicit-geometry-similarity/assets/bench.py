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

"""Benchmark harness: POTS vs ECFP / descriptors / foundation models.

Compares molecular representations on real ADMET and potency regression tasks
(ASAP antiviral, OpenADMET ExpansionRx, OpenADMET PXR) under the datasets'
provided non-random (time / scaffold) splits. The scoring metric is Spearman
rank correlation, matching the Polaris/OpenADMET challenge protocol.

Representations
---------------
- ``ecfp4``   : Morgan radius-2, 2048-bit count fingerprint (the gold standard).
- ``rdkit2d`` : ~200 RDKit physicochemical descriptors.
- ``pots``    : POTS landmark embedding -- each molecule is represented by its
                Fused Gromov-Wasserstein distance to K diverse landmark
                molecules chosen from the training set.
- ``ptsim``   : POTS used directly as a similarity metric (kNN / kernel ridge).
- ``ecfp_knn``: Tanimoto-kNN, the metric-quality baseline for ``ptsim``.

All feature-based representations feed the *same* downstream model so the
comparison isolates the representation.
"""

from __future__ import annotations

import hashlib
import os
import pickle
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

from pots import PharmacophoreCloud, fgw_distance

CACHE_DIR = os.path.join(os.path.dirname(__file__), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)

# FGW solver settings (tuned: ~0.995 corr to high-iter reference at ~2x speed).
FGW_KW = dict(alpha=0.5, epsilon=0.05, n_iter=60, sinkhorn_iter=20)


# --------------------------------------------------------------------------- #
# Cloud cache
# --------------------------------------------------------------------------- #
def _clouds_for(smiles: list[str], geometry: str = "bounds") -> list:
    """Build (and disk-cache) pharmacophore clouds for a list of SMILES."""
    key = hashlib.md5(("||".join(smiles) + geometry).encode()).hexdigest()[:16]
    path = os.path.join(CACHE_DIR, f"clouds_{geometry}_{key}.pkl")
    if os.path.exists(path):
        with open(path, "rb") as fh:
            return pickle.load(fh)
    clouds = [PharmacophoreCloud.from_smiles(s, geometry=geometry) for s in smiles]
    with open(path, "wb") as fh:
        pickle.dump(clouds, fh)
    return clouds


# --------------------------------------------------------------------------- #
# Representations
# --------------------------------------------------------------------------- #
def ecfp4(smiles: list[str], n_bits: int = 2048, radius: int = 2) -> np.ndarray:
    """Morgan count fingerprints as a dense float array."""
    gen = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    out = np.zeros((len(smiles), n_bits), dtype=np.float32)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        fp = gen.GetCountFingerprint(m)
        for idx, v in fp.GetNonzeroElements().items():
            out[i, idx] = v
    return out


def molweight(smiles: list[str]) -> np.ndarray:
    """Single-feature molecular-weight baseline (surprisingly hard to beat)."""
    out = np.zeros((len(smiles), 1), dtype=np.float32)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            out[i, 0] = Descriptors.MolWt(m)
    return out


_FM_CACHE: dict[str, dict] = {}


def fm_embedding(smiles: list[str], tag: str) -> np.ndarray:
    """Look up precomputed foundation-model embeddings (see fm_features.py).

    ``tag`` is e.g. ``"chemberta"`` or ``"molformer"``. Missing SMILES get a
    zero row.
    """
    if tag not in _FM_CACHE:
        sm = np.load(os.path.join(CACHE_DIR, f"fm_{tag}_smiles.npy"),
                     allow_pickle=True)
        _FM_CACHE[tag] = {
            "lookup": {s: i for i, s in enumerate(sm)},
            "matrix": np.load(os.path.join(CACHE_DIR, f"fm_{tag}_matrix.npy")),
        }
    lookup, mat = _FM_CACHE[tag]["lookup"], _FM_CACHE[tag]["matrix"]
    out = np.zeros((len(smiles), mat.shape[1]), dtype=np.float32)
    for i, s in enumerate(smiles):
        j = lookup.get(s)
        if j is not None:
            out[i] = mat[j]
    return out


def chemberta(smiles: list[str]) -> np.ndarray:
    return fm_embedding(smiles, "chemberta")


def molformer(smiles: list[str]) -> np.ndarray:
    return fm_embedding(smiles, "molformer")


_DESC_FNS = [f for name, f in Descriptors.descList]


def rdkit2d(smiles: list[str]) -> np.ndarray:
    """RDKit 2D physicochemical descriptors (NaN/inf sanitized)."""
    out = np.zeros((len(smiles), len(_DESC_FNS)), dtype=np.float32)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        for j, fn in enumerate(_DESC_FNS):
            try:
                out[i, j] = fn(m)
            except Exception:
                out[i, j] = 0.0
    out[~np.isfinite(out)] = 0.0
    return out


def ecfp_tanimoto_matrix(smiles_a: list[str], smiles_b: list[str]) -> np.ndarray:
    """Dense 1 - Tanimoto distance matrix between two SMILES lists (ECFP4)."""
    gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)

    def fps(sm):
        out = []
        for s in sm:
            m = Chem.MolFromSmiles(s)
            out.append(gen.GetFingerprint(m) if m else None)
        return out

    fa, fb = fps(smiles_a), fps(smiles_b)
    D = np.ones((len(fa), len(fb)), dtype=np.float32)
    for i, a in enumerate(fa):
        if a is None:
            continue
        sims = DataStructs.BulkTanimotoSimilarity(a, [b for b in fb if b is not None])
        k = 0
        for j, b in enumerate(fb):
            if b is None:
                continue
            D[i, j] = 1.0 - sims[k]
            k += 1
    return D


# --- POTS landmark embedding via parallel FGW ------------------------------ #
def _fgw_row(args):
    """Worker: FGW distance from one query cloud to all landmark clouds."""
    qi, query, landmarks, kw = args
    if query is None:
        return qi, np.full(len(landmarks), np.nan, dtype=np.float32)
    row = np.empty(len(landmarks), dtype=np.float32)
    for j, lm in enumerate(landmarks):
        row[j] = fgw_distance(query, lm, **kw) if lm is not None else np.nan
    return qi, row


def pots_landmark_embedding(
    smiles: list[str],
    train_idx: np.ndarray,
    n_landmarks: int = 64,
    geometry: str = "bounds",
    seed: int = 0,
    n_workers: int | None = None,
    return_landmarks: bool = False,
):
    """Embed molecules by FGW distance to K diverse training-set landmarks.

    Landmarks are picked with the MaxMin (sphere-exclusion) algorithm on ECFP4
    so they span chemical space. The resulting ``(n_mols, K)`` matrix is a
    fixed-length, model-ready feature set derived purely from the POTS metric.
    If ``return_landmarks``, also returns the landmarks' global indices (reused
    as inducing points for the Nyström GP).
    """
    clouds = _clouds_for(smiles, geometry)
    # Diverse landmark selection over the training subset.
    train_smiles = [smiles[i] for i in train_idx]
    gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    fps, valid_local = [], []
    for li, s in enumerate(train_smiles):
        m = Chem.MolFromSmiles(s)
        if m is not None and clouds[train_idx[li]] is not None:
            fps.append(gen.GetFingerprint(m))
            valid_local.append(li)
    k = min(n_landmarks, len(fps))
    picker = MaxMinPicker()
    picks = picker.LazyBitVectorPick(fps, len(fps), k, seed=seed)
    landmark_global = [int(train_idx[valid_local[p]]) for p in picks]
    landmarks = [clouds[g] for g in landmark_global]

    emb = np.full((len(smiles), k), np.nan, dtype=np.float32)
    tasks = [(i, clouds[i], landmarks, FGW_KW) for i in range(len(smiles))]
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for qi, row in ex.map(_fgw_row, tasks, chunksize=8):
            emb[qi] = row
    # Impute any failures with per-column medians.
    col_med = np.nanmedian(emb, axis=0)
    inds = np.where(np.isnan(emb))
    emb[inds] = np.take(col_med, inds[1])
    if return_landmarks:
        return emb, np.array(landmark_global)
    return emb


def pots_global_descriptors(smiles: list[str], geometry: str = "bounds") -> np.ndarray:
    """Cheap global descriptors of each pharmacophore cloud.

    Six pharmacophore-type counts plus summary statistics of the intra-molecular
    distance matrix (size and geometric spread). These give a downstream model
    the composition/size context that landmark FGW distances alone can wash out.
    """
    clouds = _clouds_for(smiles, geometry)
    n_types = 6
    out = np.zeros((len(smiles), n_types + 4), dtype=np.float32)
    for i, c in enumerate(clouds):
        if c is None:
            continue
        for t in c.types:
            out[i, t] += 1
        d = c.dist[np.triu_indices(len(c), 1)]
        if d.size:
            out[i, n_types:] = [len(c), d.mean(), d.std(), d.max()]
        else:
            out[i, n_types] = len(c)
    return out


def pots_augmented_embedding(smiles, train_idx, n_landmarks=64,
                             geometry="bounds", n_workers=None):
    """Landmark-FGW embedding concatenated with global cloud descriptors."""
    lm = pots_landmark_embedding(smiles, train_idx, n_landmarks=n_landmarks,
                                 geometry=geometry, n_workers=n_workers)
    gd = pots_global_descriptors(smiles, geometry=geometry)
    return np.concatenate([lm, gd], axis=1)


def _psd_clip(K: np.ndarray, floor: float = 1e-8) -> np.ndarray:
    """Project a symmetric matrix to the nearest PSD matrix (eigenvalue clip)."""
    K = 0.5 * (K + K.T)
    w, V = np.linalg.eigh(K)
    w = np.clip(w, floor, None)
    return (V * w) @ V.T


def gp_nystrom(D_tr_m, D_te_m, D_mm, y_tr,
               gammas=(0.5, 1.0, 2.0, 4.0, 8.0, 16.0),
               noises=(1e-3, 1e-2, 1e-1, 3e-1)):
    """Scalable Gaussian-process regression with a Nyström-approximated kernel.

    Uses ``M`` inducing points (landmarks) so the model scales to arbitrarily
    many training points (cost O(n·M^2 + M^3), not O(n^3)). The covariance is
    the similarity kernel ``k(A,B) = exp(-gamma·D(A,B))``; ``D_*_m`` are
    distances from points to the M landmarks and ``D_mm`` is the
    landmark-landmark distance matrix. This is the subset-of-regressors /
    Nyström sparse GP: prior weights in the Nyström feature space, exact
    Bayesian linear regression, with ``gamma`` and noise chosen by the (low-rank)
    log marginal likelihood.

    Args:
        D_tr_m: ``(n_tr, M)`` train-to-landmark distances.
        D_te_m: ``(n_te, M)`` test-to-landmark distances.
        D_mm: ``(M, M)`` landmark-to-landmark distances.
        y_tr: Training targets.

    Returns:
        Posterior-mean predictions for the test points.
    """
    mu, sd = float(np.mean(y_tr)), float(np.std(y_tr)) + 1e-9
    yz = (y_tr - mu) / sd
    n, M = D_tr_m.shape
    best = None
    for g in gammas:
        Kmm = _psd_clip(np.exp(-g * D_mm)) + 1e-6 * np.eye(M)
        L = np.linalg.cholesky(Kmm)  # Kmm = L L^T
        # Nyström features Z = K_nm @ L^{-T}, so that Z Z^T ~ K.
        Ztr = np.linalg.solve(L, np.exp(-g * D_tr_m).T).T  # (n, M)
        Zte = np.linalg.solve(L, np.exp(-g * D_te_m).T).T  # (n_te, M)
        ZtZ = Ztr.T @ Ztr
        Zty = Ztr.T @ yz
        for noise in noises:
            A = ZtZ + noise * np.eye(M)
            try:
                La = np.linalg.cholesky(A)
            except np.linalg.LinAlgError:
                continue
            w = np.linalg.solve(La.T, np.linalg.solve(La, Zty))  # A^{-1} Z^T y
            # Low-rank log marginal likelihood.
            data_fit = (yz @ yz - Zty @ w) / noise
            logdet = 2 * np.log(np.diag(La)).sum() + (n - M) * np.log(noise)
            lml = -0.5 * (data_fit + logdet + n * np.log(2 * np.pi))
            if best is None or lml > best[0]:
                best = (lml, Zte @ w)
    if best is None:
        return np.full(D_te_m.shape[0], mu)
    return best[1] * sd + mu


def gp_regression(D_tr, D_te_tr, y_tr, kind="rbf",
                  gammas=(0.5, 1.0, 2.0, 4.0, 8.0, 16.0),
                  noises=(1e-3, 1e-2, 1e-1, 3e-1)):
    """Gaussian-process regression with a distance/similarity-derived kernel.

    Uses the POTS (or Tanimoto) matrix directly as the GP covariance:
    ``k(A,B) = exp(-gamma * D(A,B))``. Because FGW is not positive
    semi-definite, the training Gram matrix is eigenvalue-clipped to the nearest
    PSD matrix. ``gamma`` and the noise variance are chosen by maximizing the
    exact GP log marginal likelihood on the training set.

    Args:
        D_tr: ``(n_tr, n_tr)`` train-train distance matrix (0 = identical).
        D_te_tr: ``(n_te, n_tr)`` test-train distance matrix.
        y_tr: Training targets.
        kind: Unused hook for alternative kernels; kept for clarity.

    Returns:
        Posterior-mean predictions for the test points.
    """
    mu, sd = float(np.mean(y_tr)), float(np.std(y_tr)) + 1e-9
    yz = (y_tr - mu) / sd
    n = len(y_tr)
    best = None
    for g in gammas:
        K = _psd_clip(np.exp(-g * D_tr))
        Kst = np.exp(-g * D_te_tr)
        for noise in noises:
            A = K + noise * np.eye(n)
            try:
                L = np.linalg.cholesky(A)
            except np.linalg.LinAlgError:
                continue
            alpha = np.linalg.solve(L.T, np.linalg.solve(L, yz))
            # log marginal likelihood
            lml = (-0.5 * yz @ alpha - np.log(np.diag(L)).sum()
                   - 0.5 * n * np.log(2 * np.pi))
            if best is None or lml > best[0]:
                best = (lml, Kst @ alpha)
    if best is None:
        return np.full(D_te_tr.shape[0], mu)
    return best[1] * sd + mu


def tanimoto_distance_matrix(smiles_a, smiles_b):
    """1 - Tanimoto (ECFP4) as a distance matrix; symmetric-PSD kernel source."""
    return ecfp_tanimoto_matrix(smiles_a, smiles_b)


def pots_distance_matrix(
    smiles: list[str],
    row_idx: np.ndarray,
    col_idx: np.ndarray,
    geometry: str = "bounds",
    n_workers: int | None = None,
    **fgw_overrides,
) -> np.ndarray:
    """Full FGW distance matrix between two index subsets (for kNN / kernels).

    ``fgw_overrides`` (e.g. ``alpha=0.25``) override the default FGW settings,
    used by the ablations.
    """
    clouds = _clouds_for(smiles, geometry)
    kw = {**FGW_KW, **fgw_overrides}
    col_clouds = [clouds[j] for j in col_idx]
    tasks = [(r, clouds[row_idx[r]], col_clouds, kw) for r in range(len(row_idx))]
    D = np.full((len(row_idx), len(col_idx)), np.nan, dtype=np.float32)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for r, row in ex.map(_fgw_row, tasks, chunksize=8):
            D[r] = row
    D[~np.isfinite(D)] = np.nanmax(D[np.isfinite(D)]) if np.isfinite(D).any() else 1.0
    return D
