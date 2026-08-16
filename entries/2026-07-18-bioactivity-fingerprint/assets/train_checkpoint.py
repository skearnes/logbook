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

"""Train the multi-task student once and save it for downstream use.

``evaluate.py`` trains its own copy in-process for the intrinsic comparison. The
transfer benchmark needs the same encoder without paying for that whole run, so
this trains and serialises it standalone.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy import sparse

import fingerprints
import student

N_WORKERS = 12
CHUNK = 10_000


def _counts_chunk(chunk: list[str]) -> np.ndarray:
    return fingerprints.count_fingerprints(chunk, "ecfp4")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eval", type=Path, required=True)
    parser.add_argument("--teacher", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--epochs", type=int, default=12)
    args = parser.parse_args()

    compounds = pd.read_parquet(args.eval / "compounds.parquet")
    profiles = sparse.load_npz(args.teacher / "profiles_raw.npz")
    molregnos = pd.read_parquet(args.teacher / "profile_molregnos.parquet")["molregno"]
    profile_row = pd.Series(np.arange(len(molregnos)), index=molregnos.to_numpy())
    profiles = profiles[compounds["molregno"].map(profile_row).to_numpy()]

    train_rows = np.flatnonzero((compounds["split"] == 0).to_numpy())
    smiles = compounds["canonical_smiles"].tolist()
    print(f"training compounds: {len(train_rows):,d}", flush=True)

    blocks = [
        [smiles[i] for i in train_rows[j:j + CHUNK]]
        for j in range(0, len(train_rows), CHUNK)
    ]
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        features = np.vstack(list(pool.map(_counts_chunk, blocks)))

    model = student.train_multitask(
        features, profiles[train_rows], n_epochs=args.epochs
    )

    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / "student_multitask.pt"
    torch.save(
        {
            "state_dict": model.state_dict(),
            "n_input": fingerprints.N_BITS,
            "n_targets": profiles.shape[1],
        },
        path,
    )
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
