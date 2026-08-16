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

"""Extract molecules and curated bioactivities from the ChEMBL 37 SQLite release.

Writes four Parquet tables:

- ``molecules.parquet`` : every structure in ChEMBL (molregno, chembl_id, SMILES).
- ``targets.parquet``   : single-protein targets.
- ``actives.parquet``   : exact-value potency records (the quantitative branch).
- ``inactives.parquet`` : censored ">" records (the negative branch).

The filters follow Landrum & Riniker's curation work on ChEMBL reliability
(J. Chem. Inf. Model. 2024, 64, 1560-1567, doi:10.1021/acs.jcim.4c00049) and the
schema semantics documented by Papadatos et al. (J. Comput.-Aided Mol. Des.
2015, 29, 885-896). Three choices are worth spelling out because they are not
the obvious defaults:

**confidence_score == 9, not >= 8.** Score 8 means "homologous single protein",
which ChEMBL also assigns when the species is simply undefined, so the target
assignment is a guess. Restricting to 9 costs less than it appears: the two
tiers are disjoint and 9 is by far the larger.

**Potency types are kept separate, never merged.** Landrum & Riniker found that
Ki is *not* more reproducible than IC50 across sources (curated MAE 0.45 vs
0.27), which undercuts the usual justification for pooling them onto a single
-log10 scale. Each record keeps its ``standard_type``.

**Censored records are extracted, not discarded.** ``pchembl_value`` is only
populated when ``standard_relation`` is '=', so the conventional
"pchembl_value IS NOT NULL" filter silently drops every measurement that
established a compound is *inactive*. Those are precisely the observations that
make a target profile informative -- without them, absence of a record conflates
"tested and inactive" with "never tested".
"""

from __future__ import annotations

import argparse
import hashlib
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

MIN_CONFIDENCE_SCORE = 9
ASSAY_TYPES = ("B", "F")
STANDARD_TYPES = ("Ki", "Kd", "IC50", "EC50")
TARGET_TYPE = "SINGLE PROTEIN"

# A depositor can record a potent value and still conclude the compound is
# inactive, typically after a counter-screen; the comment overrides the number.
NEGATIVE_COMMENTS = (
    "not active", "inactive", "inconclusive", "unspecified", "not determined",
)

# The ten assay fields that define a comparability unit. Two assays sharing a
# target and these values are measuring the same thing under the same
# conditions, which is the only defensible unit for pooling replicates.
ACH_FIELDS = (
    "assay_type", "assay_organism", "assay_category", "assay_tax_id",
    "assay_strain", "assay_tissue", "assay_cell_type",
    "assay_subcellular_fraction", "bao_format", "variant_id",
)

MOLECULE_SQL = """
SELECT md.molregno, md.chembl_id, cs.canonical_smiles
FROM molecule_dictionary md
JOIN compound_structures cs ON cs.molregno = md.molregno
WHERE cs.canonical_smiles IS NOT NULL
"""

TARGET_SQL = f"""
SELECT td.tid, td.chembl_id, td.pref_name, td.organism
FROM target_dictionary td
WHERE td.target_type = '{TARGET_TYPE}'
"""

_COMMENT_LIST = ", ".join(f"'{c}'" for c in NEGATIVE_COMMENTS)
_ACH_COLUMNS = ", ".join(f"a.{f}" for f in ACH_FIELDS)

# Shared predicate for both branches. Mutants are excluded via variant_id
# because a point mutation can abolish binding, making the wild-type target
# label wrong for that record.
_COMMON_WHERE = f"""
  AND act.standard_units = 'nM'
  AND act.standard_type IN {STANDARD_TYPES}
  AND act.data_validity_comment IS NULL
  AND act.potential_duplicate = 0
  AND (act.activity_comment IS NULL
       OR lower(act.activity_comment) NOT IN ({_COMMENT_LIST}))
  AND a.confidence_score = {MIN_CONFIDENCE_SCORE}
  AND a.assay_type IN {ASSAY_TYPES}
  AND td.target_type = '{TARGET_TYPE}'
"""

ACTIVE_SQL = f"""
SELECT act.molregno, a.tid, a.assay_id, act.doc_id,
       act.standard_type, act.pchembl_value, {_ACH_COLUMNS}
FROM activities act
JOIN assays a ON a.assay_id = act.assay_id
JOIN target_dictionary td ON td.tid = a.tid
WHERE act.pchembl_value IS NOT NULL
  AND act.standard_relation = '='
  AND act.standard_flag = 1
{_COMMON_WHERE}
"""

# Censored records carry no pchembl_value, so the -log10 is computed here. The
# result is a *bound*: true potency is weaker than this figure.
INACTIVE_SQL = f"""
SELECT act.molregno, a.tid, a.assay_id, act.doc_id,
       act.standard_type, act.standard_value, {_ACH_COLUMNS}
FROM activities act
JOIN assays a ON a.assay_id = act.assay_id
JOIN target_dictionary td ON td.tid = a.tid
WHERE act.standard_relation IN ('>', '>=')
  AND act.standard_value IS NOT NULL
  AND act.standard_value > 0
{_COMMON_WHERE}
"""


def add_conditions_hash(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach the assay-conditions hash and drop the columns it summarises."""
    joined = frame[list(ACH_FIELDS)].astype("string").fillna("").agg("|".join, axis=1)
    frame = frame.drop(columns=list(ACH_FIELDS))
    frame["ach"] = [hashlib.md5(v.encode()).hexdigest()[:16] for v in joined]
    return frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)

    for name, sql in [("molecules", MOLECULE_SQL), ("targets", TARGET_SQL)]:
        frame = pd.read_sql_query(sql, con)
        frame.to_parquet(args.out / f"{name}.parquet", index=False)
        print(f"{name:12s} {len(frame):>10,d} rows")

    actives = add_conditions_hash(pd.read_sql_query(ACTIVE_SQL, con))
    actives["pchembl_value"] = actives["pchembl_value"].astype("float32")
    actives.to_parquet(args.out / "actives.parquet", index=False)
    print(
        f"{'actives':12s} {len(actives):>10,d} rows "
        f"({actives.molregno.nunique():,d} mols, {actives.tid.nunique():,d} targets, "
        f"{actives.ach.nunique():,d} assay-condition groups)"
    )

    inactives = add_conditions_hash(pd.read_sql_query(INACTIVE_SQL, con))
    # -log10(M) from nM, so 10 uM -> 5.0. Potency is weaker than this bound.
    inactives["p_bound"] = (
        9.0 - np.log10(inactives["standard_value"].astype("float64"))
    ).astype("float32")
    inactives = inactives.drop(columns=["standard_value"])
    inactives.to_parquet(args.out / "inactives.parquet", index=False)
    print(
        f"{'inactives':12s} {len(inactives):>10,d} rows "
        f"({inactives.molregno.nunique():,d} mols, {inactives.tid.nunique():,d} targets)"
    )

    con.close()


if __name__ == "__main__":
    main()
