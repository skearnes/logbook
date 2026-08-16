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

"""Precompute foundation-model molecular embeddings and cache them to .npy.

Runs in the torch-enabled env; the benchmark (base env) loads the cached arrays
by SMILES-content hash. Keeps the heavy dependency out of the main harness.

Models:
- ChemBERTa (``DeepChem/ChemBERTa-77M-MTR``): a RoBERTa SMILES transformer.
  We mean-pool the last hidden state over non-pad tokens.
"""

from __future__ import annotations

import hashlib
import os
import sys

import numpy as np

CACHE = os.path.join(os.path.dirname(__file__), "cache")
os.makedirs(CACHE, exist_ok=True)

MODELS = {
    # (hf_repo, trust_remote_code) -- ChemBERTa is an older RoBERTa SMILES model;
    # MoLFormer-XL is a stronger, more recent SMILES transformer.
    "chemberta": ("DeepChem/ChemBERTa-77M-MTR", False),
    "molformer": ("ibm-research/MoLFormer-XL-both-10pct", True),
}


def _key(smiles, tag):
    h = hashlib.md5(("||".join(smiles)).encode()).hexdigest()[:16]
    return os.path.join(CACHE, f"fm_{tag}_{h}.npy")


def embed(smiles: list[str], tag: str = "chemberta", batch_size: int = 64):
    """Return (and cache) mean-pooled embeddings for a list of SMILES."""
    path = _key(smiles, tag)
    if os.path.exists(path):
        return np.load(path)
    import torch
    from transformers import AutoModel, AutoTokenizer

    name, trust = MODELS[tag]
    tok = AutoTokenizer.from_pretrained(name, trust_remote_code=trust)
    model = AutoModel.from_pretrained(name, trust_remote_code=trust).eval()
    vecs = []
    with torch.no_grad():
        for i in range(0, len(smiles), batch_size):
            batch = smiles[i : i + batch_size]
            enc = tok(batch, padding=True, truncation=True, max_length=256,
                      return_tensors="pt")
            out = model(**enc).last_hidden_state  # (B, T, H)
            mask = enc["attention_mask"].unsqueeze(-1).float()
            pooled = (out * mask).sum(1) / mask.sum(1).clamp(min=1)
            vecs.append(pooled.cpu().numpy())
            print(f"  {tag}: {min(i+batch_size, len(smiles))}/{len(smiles)}",
                  flush=True)
    emb = np.concatenate(vecs).astype(np.float32)
    np.save(path, emb)
    return emb


if __name__ == "__main__":
    # Precompute for every SMILES column across all benchmark CSVs.
    import glob
    import pandas as pd

    data = os.path.join(os.path.dirname(__file__), "data")
    all_smiles = set()
    for csv in glob.glob(os.path.join(data, "*.csv")):
        try:
            df = pd.read_csv(csv)
        except Exception:
            continue
        col = "SMILES" if "SMILES" in df.columns else None
        if col:
            all_smiles.update(df[col].dropna().tolist())
    smiles = sorted(all_smiles)
    tags = sys.argv[1:] or list(MODELS)
    for tag in tags:
        print(f"Embedding {len(smiles)} unique SMILES with {tag}", flush=True)
        emb = embed(smiles, tag)
        # Save a SMILES->row lookup so the base-env harness can index by SMILES.
        np.save(os.path.join(CACHE, f"fm_{tag}_smiles.npy"),
                np.array(smiles, dtype=object), allow_pickle=True)
        np.save(os.path.join(CACHE, f"fm_{tag}_matrix.npy"), emb)
        print(f"Saved {tag} embeddings {emb.shape}", flush=True)
