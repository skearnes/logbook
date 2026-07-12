"""POTS: Pharmacophore Optimal-Transport Similarity.

A conformer-free, interpretable molecular similarity built on Fused
Gromov-Wasserstein (FGW) optimal transport between *pharmacophore point clouds*.

Motivation
----------
A Morgan/ECFP fingerprint is a bag of hashed atom-centered subgraphs compared
with Tanimoto. That similarity is (a) not size-robust -- small molecules set
fewer bits and look systematically less similar, (b) a single global scalar that
cannot say *which* parts of two molecules correspond, and (c) blind to shape and
electrostatics, which is what actually governs protein-ligand binding.

POTS reframes a molecule as a small cloud of pharmacophore features (H-bond
donors/acceptors, aromatics, hydrophobes, cations, anions) embedded in an
*implicit* geometry given by a cheap distance-geometry bounds matrix -- no
conformer generation, no docking. Two molecules are compared by the Fused
Gromov-Wasserstein distance between their clouds, which simultaneously matches
feature *types* (the "color" term) and preserves *internal distances* (the
"shape" term). The optimal transport plan is a soft maximum-common-pharmacophore
alignment: it says which feature of A maps to which feature of B, and how much
mass is left unmatched -- a directly interpretable, size-aware readout.

Because Gromov-Wasserstein compares molecules through their internal distance
matrices, it never needs a shared coordinate frame: the 3D-ness is implicit.

Public API
----------
- ``PharmacophoreCloud.from_mol`` / ``from_smiles``: build a cloud.
- ``fgw_distance``: FGW distance between two clouds (plus the transport plan).
- ``pots_distance`` / ``pots_similarity``: convenience wrappers on SMILES.
"""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures, rdDistGeom

# The six classical pharmacophore families, plus a graded type-mismatch cost.
# Families not in this map (e.g. ZnBinder, LumpedHydrophobe) are folded onto a
# representative below.
_FAMILY_CANON = {
    "Donor": "Donor",
    "Acceptor": "Acceptor",
    "Aromatic": "Aromatic",
    "Hydrophobe": "Hydrophobe",
    "LumpedHydrophobe": "Hydrophobe",
    "PosIonizable": "Cation",
    "NegIonizable": "Anion",
    "ZnBinder": "Anion",
}
_TYPES = ("Donor", "Acceptor", "Aromatic", "Hydrophobe", "Cation", "Anion")
_TYPE_IDX = {t: i for i, t in enumerate(_TYPES)}

# Graded pharmacophore-type dissimilarity in [0, 1]. Same type = 0. Chemically
# related types (donor<->acceptor, aromatic<->hydrophobe, cation<->donor,
# anion<->acceptor) cost less than unrelated ones. This makes the "color" term
# smoother than a hard 0/1 indicator so near-miss feature swaps are penalized
# gently.
def _build_type_cost() -> np.ndarray:
    c = np.ones((len(_TYPES), len(_TYPES)), dtype=float)
    np.fill_diagonal(c, 0.0)
    related = {
        ("Donor", "Acceptor"): 0.5,
        ("Aromatic", "Hydrophobe"): 0.4,
        ("Cation", "Donor"): 0.5,
        ("Anion", "Acceptor"): 0.5,
        ("Cation", "Anion"): 1.0,
    }
    for (a, b), v in related.items():
        i, j = _TYPE_IDX[a], _TYPE_IDX[b]
        c[i, j] = c[j, i] = v
    return c


_TYPE_COST = _build_type_cost()

_FDEF = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
_FACTORY = ChemicalFeatures.BuildFeatureFactory(_FDEF)


@dataclass
class PharmacophoreCloud:
    """A molecule as a weighted cloud of pharmacophore features.

    Attributes:
        types: Integer type id per feature (index into ``_TYPES``).
        weights: Non-negative mass per feature, summing to 1.
        dist: ``(n, n)`` intra-molecular feature-feature distance matrix.
        smiles: Source SMILES, for debugging and provenance.
        families: Human-readable type name per feature.
    """

    types: np.ndarray
    weights: np.ndarray
    dist: np.ndarray
    smiles: str
    families: list[str]

    def __len__(self) -> int:
        return len(self.types)

    @staticmethod
    def from_smiles(smiles: str, **kwargs) -> "PharmacophoreCloud | None":
        """Build a cloud from a SMILES string, or None if it cannot be parsed."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return PharmacophoreCloud.from_mol(mol, smiles=smiles, **kwargs)

    @staticmethod
    def from_mol(
        mol: Chem.Mol,
        *,
        smiles: str | None = None,
        geometry: str = "bounds",
    ) -> "PharmacophoreCloud | None":
        """Build a cloud from an RDKit mol.

        Args:
            mol: Input molecule (Hs are added internally for geometry).
            smiles: Optional provenance string.
            geometry: Source of the intra-molecular distances. One of
                ``"bounds"`` (distance-geometry bounds midpoints; conformer-free,
                the default), ``"topology"`` (shortest-path bond counts), or
                ``"conformer"`` (a single ETKDG 3D conformer).

        Returns:
            A cloud, or None if no pharmacophore features are found.
        """
        if smiles is None:
            smiles = Chem.MolToSmiles(mol)
        molh = Chem.AddHs(mol)
        feats = _FACTORY.GetFeaturesForMol(molh)
        rows = []
        for f in feats:
            fam = _FAMILY_CANON.get(f.GetFamily())
            if fam is None:
                continue
            rows.append((_TYPE_IDX[fam], fam, list(f.GetAtomIds())))
        if not rows:
            return None
        types = np.array([r[0] for r in rows], dtype=int)
        families = [r[1] for r in rows]
        atom_sets = [r[2] for r in rows]
        weights = np.full(len(rows), 1.0 / len(rows))
        dist = _feature_distances(molh, atom_sets, geometry)
        return PharmacophoreCloud(types, weights, dist, smiles, families)


def _feature_distances(
    molh: Chem.Mol, atom_sets: list[list[int]], geometry: str
) -> np.ndarray:
    """Compute the pairwise distance matrix between pharmacophore features.

    Each feature is a set of atoms; the feature-feature distance is the mean of
    the underlying atom-atom distances across the two sets.
    """
    n = len(atom_sets)
    if geometry == "topology":
        dmat = Chem.GetDistanceMatrix(molh).astype(float)
    elif geometry == "bounds":
        bm = rdDistGeom.GetMoleculeBoundsMatrix(molh)
        # Symmetric midpoint of the lower/upper distance bounds.
        lower = np.tril(bm, -1)
        upper = np.triu(bm, 1)
        dmat = (upper + upper.T + lower + lower.T) / 2.0
    elif geometry == "conformer":
        from rdkit.Chem import AllChem

        conf_mol = Chem.Mol(molh)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D
        if AllChem.EmbedMolecule(conf_mol, params) != 0:
            dmat = Chem.GetDistanceMatrix(molh).astype(float)
        else:
            conf = conf_mol.GetConformer()
            pos = conf.GetPositions()
            diff = pos[:, None, :] - pos[None, :, :]
            dmat = np.sqrt((diff**2).sum(-1))
    else:
        raise ValueError(f"unknown geometry {geometry!r}")

    out = np.zeros((n, n))
    for i in range(n):
        ai = atom_sets[i]
        for j in range(i + 1, n):
            aj = atom_sets[j]
            block = dmat[np.ix_(ai, aj)]
            out[i, j] = out[j, i] = float(block.mean())
    return out


def _init_plan(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Independent-coupling initialization ``a b^T``."""
    return np.outer(a, b)


def fgw_distance(
    ca: PharmacophoreCloud,
    cb: PharmacophoreCloud,
    *,
    alpha: float = 0.5,
    epsilon: float = 0.05,
    n_iter: int = 200,
    sinkhorn_iter: int = 50,
    tol: float = 1e-9,
    return_plan: bool = False,
):
    """Entropic Fused Gromov-Wasserstein distance between two clouds.

    Minimizes ``(1 - alpha) * <M, T> + alpha * GW_quadratic(D_a, D_b, T)`` over
    couplings ``T`` with the given marginals, using entropic mirror-descent
    (Peyre et al. 2016; Vayer et al. 2019). ``M`` is the pharmacophore type-cost
    and the quadratic term is the squared distance-distortion between the two
    intra-molecular distance matrices.

    Args:
        ca: First pharmacophore cloud.
        cb: Second pharmacophore cloud.
        alpha: Weight on the structural (shape) term vs the feature (color)
            term. ``alpha=0`` is pure Wasserstein on types; ``alpha=1`` is pure
            Gromov-Wasserstein on geometry.
        epsilon: Entropic regularization strength.
        n_iter: Outer mirror-descent iterations.
        sinkhorn_iter: Inner Sinkhorn iterations per outer step.
        tol: Convergence tolerance on the transport plan.
        return_plan: If True, also return the coupling matrix.

    Returns:
        The FGW cost (float). If ``return_plan``, a ``(cost, plan)`` tuple.
    """
    a, b = ca.weights, cb.weights
    Da, Db = ca.dist, cb.dist
    # Normalize distance scales so alpha is comparable across molecule sizes.
    sa = Da[Da > 0].mean() if (Da > 0).any() else 1.0
    sb = Db[Db > 0].mean() if (Db > 0).any() else 1.0
    Da = Da / sa
    Db = Db / sb

    # Color cost matrix M from the graded type-cost lookup.
    M = _TYPE_COST[np.ix_(ca.types, cb.types)]

    # Constants for the GW quadratic-term gradient (Proposition 1, Peyre 2016)
    # with the squared-loss L(x,y) = (x - y)^2:
    #   grad = constC - 2 * Da @ T @ Db^T
    # where constC = (Da^2 a) 1^T + 1 (Db^2 b)^T.
    Da2, Db2 = Da**2, Db**2
    constC = np.outer(Da2 @ a, np.ones_like(b)) + np.outer(
        np.ones_like(a), Db2 @ b
    )

    T = _init_plan(a, b)
    log_a, log_b = np.log(a), np.log(b)
    for _ in range(n_iter):
        gw_grad = constC - 2.0 * (Da @ T @ Db.T)
        grad = (1.0 - alpha) * M + alpha * gw_grad
        # One entropic (Sinkhorn) projection step of mirror descent.
        K = -grad / epsilon
        K -= K.max()
        logK = K
        f = np.zeros_like(a)
        g = np.zeros_like(b)
        for _ in range(sinkhorn_iter):
            f = log_a - _logsumexp(logK + g[None, :], axis=1)
            g = log_b - _logsumexp(logK + f[:, None], axis=0)
        T_new = np.exp(f[:, None] + logK + g[None, :])
        if np.abs(T_new - T).sum() < tol:
            T = T_new
            break
        T = T_new

    gw_term = np.sum((constC - 2.0 * (Da @ T @ Db.T)) * T)
    color_term = np.sum(M * T)
    cost = (1.0 - alpha) * color_term + alpha * gw_term
    cost = max(cost, 0.0)
    if return_plan:
        return cost, T
    return cost


def _logsumexp(x: np.ndarray, axis: int) -> np.ndarray:
    m = x.max(axis=axis, keepdims=True)
    return (m + np.log(np.exp(x - m).sum(axis=axis, keepdims=True))).squeeze(axis)


def pots_distance(smiles_a: str, smiles_b: str, **kwargs) -> float:
    """FGW distance between two molecules given as SMILES."""
    geometry = kwargs.pop("geometry", "bounds")
    ca = PharmacophoreCloud.from_smiles(smiles_a, geometry=geometry)
    cb = PharmacophoreCloud.from_smiles(smiles_b, geometry=geometry)
    if ca is None or cb is None:
        return float("nan")
    return fgw_distance(ca, cb, **kwargs)


def pots_similarity(smiles_a: str, smiles_b: str, gamma: float = 4.0, **kwargs) -> float:
    """A bounded similarity ``exp(-gamma * FGW)`` in ``(0, 1]``."""
    d = pots_distance(smiles_a, smiles_b, **kwargs)
    if np.isnan(d):
        return float("nan")
    return float(np.exp(-gamma * d))
