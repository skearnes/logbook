"""Structural representations used as baselines and as student inputs.

Every representation here is computable from a SMILES string alone, which is the
constraint the whole experiment operates under: the student must work on
molecules that have never been assayed.

The 3D entries generate a single ETKDG conformer rather than an ensemble. That
follows the conformer ablation in the 2026-07-11 entry, where one ETKDG
conformer outperformed a distance-bounds "implicit geometry" -- ensembles cost
far more and were not what carried the signal there.
"""

from __future__ import annotations

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.rdDistGeom import EmbedMolecule, ETKDGv3

# Parse failures are expected across a million ChEMBL structures and are handled
# by returning zero vectors; the per-molecule warnings are pure noise.
RDLogger.DisableLog("rdApp.*")

N_BITS = 2048
ETKDG_SEED = 0xF00D


def _generator(kind: str):
    if kind == "ecfp4":
        return rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=N_BITS)
    if kind == "ecfp6":
        return rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=N_BITS)
    if kind == "fcfp4":
        invariants = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
        return rdFingerprintGenerator.GetMorganGenerator(
            radius=2, fpSize=N_BITS, atomInvariantsGenerator=invariants
        )
    if kind == "rdkit":
        return rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=N_BITS)
    if kind == "atompair":
        return rdFingerprintGenerator.GetAtomPairGenerator(fpSize=N_BITS)
    if kind == "topotorsion":
        return rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=N_BITS)
    raise ValueError(f"unknown fingerprint: {kind}")


def bit_fingerprints(smiles: list[str], kind: str = "ecfp4") -> np.ndarray:
    """Binary fingerprints as a packed uint8 array of shape (n, N_BITS)."""
    generator = _generator(kind)
    out = np.zeros((len(smiles), N_BITS), dtype=np.uint8)
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        DataStructs.ConvertToNumpyArray(generator.GetFingerprint(mol), out[i])
    return out


def count_fingerprints(smiles: list[str], kind: str = "ecfp4") -> np.ndarray:
    """Count fingerprints as float32, the usual input for a learned model."""
    generator = _generator(kind)
    out = np.zeros((len(smiles), N_BITS), dtype=np.float32)
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        for idx, count in generator.GetCountFingerprint(mol).GetNonzeroElements().items():
            out[i, idx] = count
    return out


# Ipc grows factorially with graph size and overflows even float64 on large
# molecules, poisoning every downstream statistic with inf/NaN.
UNSTABLE_DESCRIPTORS = {"Ipc"}

# Retained descriptors are clipped to a range that comfortably covers real
# physicochemical values, so a single pathological molecule cannot dominate
# standardisation.
DESCRIPTOR_CLIP = 1e12


def descriptors_2d(smiles: list[str]) -> np.ndarray:
    """The RDKit physicochemical descriptor block, numerically sanitised."""
    entries = [
        (name, fn) for name, fn in Descriptors.descList
        if name not in UNSTABLE_DESCRIPTORS
    ]
    # Accumulated in float64: several descriptors exceed the float32 range on
    # large molecules, and the overflow silently becomes inf.
    out = np.zeros((len(smiles), len(entries)), dtype=np.float64)
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        for j, (_, fn) in enumerate(entries):
            try:
                out[i, j] = fn(mol)
            except Exception:
                # A handful of descriptors throw on exotic valences; leaving the
                # entry at zero is preferable to dropping the whole molecule.
                pass
    out = np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)
    return np.clip(out, -DESCRIPTOR_CLIP, DESCRIPTOR_CLIP).astype(np.float32)


def _embed(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = ETKDGv3()
    params.randomSeed = ETKDG_SEED
    if EmbedMolecule(mol, params) != 0:
        return None
    return mol


def usr_descriptors(smiles: list[str], catz: bool = True) -> np.ndarray:
    """USR (12-dim) or USRCAT (60-dim) shape descriptors from one conformer.

    USRCAT extends USR with pharmacophore-typed atom subsets, so it carries
    "colour" as well as shape and is the closer analogue to a ROCS-style score.
    """
    width = 60 if catz else 12
    out = np.zeros((len(smiles), width), dtype=np.float32)
    for i, smi in enumerate(smiles):
        mol = _embed(smi)
        if mol is None:
            continue
        try:
            values = (
                rdMolDescriptors.GetUSRCAT(mol) if catz else rdMolDescriptors.GetUSR(mol)
            )
            out[i] = np.asarray(values, dtype=np.float32)
        except Exception:
            pass
    return out


def tanimoto_matrix(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Dense Tanimoto between two sets of binary fingerprints."""
    left = left.astype(np.float32)
    right = right.astype(np.float32)
    intersection = left @ right.T
    counts_left = left.sum(axis=1)[:, None]
    counts_right = right.sum(axis=1)[None, :]
    union = counts_left + counts_right - intersection
    return np.divide(
        intersection, union, out=np.zeros_like(intersection), where=union > 0
    )


def tanimoto_pairs(fps: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    """Tanimoto for specific row pairs, avoiding an all-vs-all matrix."""
    left = fps[pairs[:, 0]].astype(np.float32)
    right = fps[pairs[:, 1]].astype(np.float32)
    intersection = (left * right).sum(axis=1)
    union = left.sum(axis=1) + right.sum(axis=1) - intersection
    return np.divide(
        intersection, union, out=np.zeros_like(intersection), where=union > 0
    )


def cosine_pairs(features: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    """Cosine similarity for specific row pairs of a dense feature matrix."""
    left = features[pairs[:, 0]]
    right = features[pairs[:, 1]]
    norms = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    dots = (left * right).sum(axis=1)
    return np.divide(dots, norms, out=np.zeros_like(dots), where=norms > 0)
