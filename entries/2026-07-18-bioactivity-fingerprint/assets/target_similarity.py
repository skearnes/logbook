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

"""Build a structure-independent target-target relatedness matrix.

Most ChEMBL compounds are annotated against a single target, so a bioactivity
similarity built on raw target-set overlap is degenerate: two compounds either
share their one target or score zero. Relating the targets themselves rescues
those pairs -- two compounds hitting different but homologous kinases are more
alike than two hitting a kinase and a GPCR.

The relatedness must not be derived from ligand chemistry. SEA-style
target-target similarity is computed from the Tanimoto overlap of ligand sets,
which would smuggle chemical structure into a ground truth used to score a
structure-only model, inflating it for free. Both sources here are independent
of ligand structure:

- ``sequence`` : all-vs-all protein sequence identity from MMseqs2.
- ``family``   : depth of shared lineage in ChEMBL's curated protein
                 classification hierarchy.

The two are combined by taking the stronger signal, so homologues are related
even when the classification is coarse, and family members are related even
when sequence identity is low.
"""

from __future__ import annotations

import argparse
import shutil
import sqlite3
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

# Sequence identity below this is indistinguishable from background similarity
# between unrelated proteins, so it is treated as unrelated.
MIN_SEQUENCE_IDENTITY = 0.30

# ChEMBL protein classes run from level 1 (e.g. "Enzyme") to level 8. Shared
# lineage is scored by depth, normalised against this ceiling.
MAX_CLASS_DEPTH = 6

# Sharing only "Enzyme" or "Transporter" relates a third of all target pairs,
# which would make every compound weakly similar to every other and drown the
# signal. Three shared levels is roughly the family level (e.g.
# Enzyme -> Kinase -> Protein Kinase) and is the point at which shared lineage
# starts to imply shared ligand recognition.
MIN_CLASS_DEPTH = 3

SEQUENCE_SQL = """
SELECT tc.tid, cs.sequence
FROM target_components tc
JOIN component_sequences cs ON cs.component_id = tc.component_id
WHERE cs.sequence IS NOT NULL
"""

# protein_class_desc is the full lineage as a delimited path, so shared-prefix
# depth is read straight off the string without walking parent_id links.
CLASS_SQL = """
SELECT tc.tid, pc.protein_class_desc
FROM target_components tc
JOIN component_class cc ON cc.component_id = tc.component_id
JOIN protein_classification pc ON pc.protein_class_id = cc.protein_class_id
"""


def load_sequences(con: sqlite3.Connection, tids: set[int]) -> dict[int, str]:
    """Return the longest sequence for each requested target."""
    frame = pd.read_sql_query(SEQUENCE_SQL, con)
    frame = frame[frame["tid"].isin(tids)]
    # A handful of targets map to several components; the longest is the best
    # single proxy for the binding entity.
    frame["length"] = frame["sequence"].str.len()
    frame = frame.sort_values("length").drop_duplicates("tid", keep="last")
    return dict(zip(frame["tid"], frame["sequence"]))


def sequence_identity(sequences: dict[int, str], threads: int) -> pd.DataFrame:
    """All-vs-all sequence identity via MMseqs2, as a tid_1/tid_2/identity frame."""
    if shutil.which("mmseqs") is None:
        raise RuntimeError("mmseqs not found on PATH")

    work = Path(tempfile.mkdtemp(prefix="mmseqs_"))
    fasta = work / "targets.fasta"
    with fasta.open("w") as fh:
        for tid, seq in sequences.items():
            fh.write(f">{tid}\n{seq}\n")

    result = work / "result.tsv"
    subprocess.run(
        [
            "mmseqs", "easy-search", str(fasta), str(fasta), str(result), str(work / "tmp"),
            "--format-output", "query,target,fident",
            # Exhaustive search: the target set is small and missing a homologous
            # pair silently zeroes a similarity we specifically want to capture.
            "-s", "7.5",
            "--max-seqs", str(len(sequences)),
            "-e", "1e-3",
            "--threads", str(threads),
        ],
        check=True,
        capture_output=True,
    )

    hits = pd.read_csv(result, sep="\t", names=["tid_1", "tid_2", "identity"])
    shutil.rmtree(work, ignore_errors=True)
    hits = hits[hits["identity"] >= MIN_SEQUENCE_IDENTITY]
    return hits


def family_similarity(con: sqlite3.Connection, tids: set[int]) -> pd.DataFrame:
    """Relatedness from shared depth in the ChEMBL protein class hierarchy."""
    frame = pd.read_sql_query(CLASS_SQL, con)
    frame = frame[frame["tid"].isin(tids)]

    lineages: dict[int, list[tuple[str, ...]]] = {}
    for tid, desc in zip(frame["tid"], frame["protein_class_desc"]):
        lineages.setdefault(tid, []).append(tuple(desc.split("  ")))

    # Group targets by their level-1 class so only plausibly-related pairs are
    # compared; cross-group pairs share no lineage and score zero by definition.
    buckets: dict[str, list[int]] = {}
    for tid, paths in lineages.items():
        for path in paths:
            buckets.setdefault(path[0], []).append(tid)

    rows = []
    for members in buckets.values():
        members = sorted(set(members))
        for i, tid_1 in enumerate(members):
            for tid_2 in members[i + 1:]:
                depth = max(
                    _shared_prefix(p_1, p_2)
                    for p_1 in lineages[tid_1]
                    for p_2 in lineages[tid_2]
                )
                if depth >= MIN_CLASS_DEPTH:
                    rows.append((tid_1, tid_2, min(depth, MAX_CLASS_DEPTH) / MAX_CLASS_DEPTH))

    return pd.DataFrame(rows, columns=["tid_1", "tid_2", "identity"])


def _shared_prefix(path_1: tuple[str, ...], path_2: tuple[str, ...]) -> int:
    depth = 0
    for part_1, part_2 in zip(path_1, path_2):
        if part_1 != part_2:
            break
        depth += 1
    return depth


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--activities", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    activities = pd.read_parquet(args.activities, columns=["tid"])
    tids = set(activities["tid"].unique().tolist())
    print(f"curated targets: {len(tids):,d}")

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    sequences = load_sequences(con, tids)
    print(f"sequences: {len(sequences):,d}")

    by_sequence = sequence_identity(sequences, args.threads)
    print(f"sequence pairs >= {MIN_SEQUENCE_IDENTITY}: {len(by_sequence):,d}")

    by_family = family_similarity(con, tids)
    print(f"family pairs: {len(by_family):,d}")
    con.close()

    index = {tid: i for i, tid in enumerate(sorted(tids))}
    size = len(index)

    def to_matrix(frame: pd.DataFrame) -> sparse.csr_matrix:
        rows = frame["tid_1"].map(index).to_numpy()
        cols = frame["tid_2"].map(index).to_numpy()
        keep = ~(pd.isna(rows) | pd.isna(cols))
        return sparse.csr_matrix(
            (frame["identity"].to_numpy()[keep],
             (rows[keep].astype(int), cols[keep].astype(int))),
            shape=(size, size), dtype=np.float32,
        )

    # Stronger of the two signals wins, so neither source can veto the other:
    # homologues stay related under a coarse classification, and family members
    # stay related despite low sequence identity.
    matrix = to_matrix(by_family).maximum(to_matrix(by_sequence))
    matrix = matrix.maximum(matrix.T)
    matrix.setdiag(1.0)
    matrix.eliminate_zeros()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(args.out, matrix)
    order = pd.DataFrame({"tid": sorted(tids)})
    order.to_parquet(args.out.with_suffix(".tids.parquet"), index=False)

    density = matrix.nnz / size**2
    print(f"matrix {size}x{size}, nnz={matrix.nnz:,d} (density {density:.4f}) -> {args.out}")


if __name__ == "__main__":
    main()
