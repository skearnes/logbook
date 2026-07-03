# 2026-07-02 — Reaction classification against the RXNO ontology

## Question

The ORD pipeline currently uses [Rxn-INSIGHT](https://github.com/mrodobbe/Rxn-INSIGHT)
for reaction classification and naming. I'd prefer output that maps onto the **RXNO**
(Name Reaction Ontology) — or NameRxn, but that's proprietary. What would it take to
build a classification system/model that emits RXNO categories? This entry scopes the
landscape: how Rxn-INSIGHT works today, what RXNO/NameRxn/IBM RXN actually are, whether
the newer rxnmapper helps, and the realistic options.

Follow-ups (2026-07-02) widened the scope to two questions that reshaped the answer:
(1) is `RXNMapper_v2` a genuine improvement over the original mapper, and (2) is it
acceptable to target **NameRxn** classes instead of RXNO if a good open model already
exists or enough public labels exist to train one? See the two follow-up sections and the
revised Summary.

## Summary

**RXNO is a vocabulary, not a classifier**, and NameRxn — which emits RXNO IDs directly
alongside its own `N.N.N` codes — is proprietary. My original read was that this left
only a hand-built crosswalk or a NameRxn license. The follow-up research (2026-07-02,
new sections below) found a **better open path**, and it hinges on one verified fact:

**The public Schneider-50k dataset is labeled with *real* NameRxn `N.N.N` codes.** I
fetched the actual file: `schneider50k.tsv`'s `rxn_class` column holds values like
`6.1.5`, `7.1.1`, `1.8.5`, `1.7.9` — genuine NameRxn hierarchical codes, MIT-licensed,
not anonymized integers. Those codes map straight into RXNO. So an open classifier
**already gives you real NameRxn/RXNO labels today** — the catch is a **50-leaf-class
ceiling** (9 of the 10 superclasses, 28 categories, 50 hand-picked leaves). The full
~1000-leaf space has no open dataset; that still needs proprietary NameRxn/Pistachio.

Revised recommendation, by the resolution you actually need:

1. **~50 NameRxn leaf classes is enough → train/run on Schneider-50k (recommended).**
   Pretrained **rxnfp** ships a Schneider classifier (~98%) that emits real NameRxn
   `N.N.N` codes you can run today; **DRFP** or **SynCat** (2026 SOTA) train to the same
   50 classes. All MIT. This is a clean jump from today's 10 Rxn-INSIGHT superclasses to
   50 named leaves (plus the superclass/category levels for free by truncating the code),
   and — because the labels *are* NameRxn codes — it **supersedes the original
   "crosswalk Rxn-INSIGHT names → RXNO" idea** for emitting machine-readable classes.
2. **Want fine (~1500-class), open, near-NameRxn resolution → evaluate the Schwaller
   ReactionClassifier** (2026, MIT, arXiv 2607.01061; Rxn-INSIGHT's own author is a
   coauthor). An agentic-LLM pipeline that self-expands a *verified* SMIRKS taxonomy to
   ~1,546 L3 types / 14k classes and matches NameRxn accuracy (~97.7%). It emits its
   **own** dotted codes, not NameRxn's/RXNO IDs, but at near-NameRxn resolution, fully
   open. Brand-new (1-day-old preprint); the LLM expansion layer needs a Gemini API key,
   though the deterministic classifier + released rules run offline.
3. **Need the full ~1000 NameRxn leaves with exact codes/RXNO IDs → license NameRxn.**
   Still the only turnkey path to complete leaf-level coverage.

The original RXNO-crosswalk-of-Rxn-INSIGHT-names (detailed below, now a **fallback**)
stays useful only if you specifically want to keep Rxn-INSIGHT's human-readable named
reactions; for emitting codes it's dominated by option 1.

On atom mapping: the mapper the question pointed to is **RXNMapper_v2** (a genuine
retrain, not the `rxn4chemistry/rxnmapper` 0.4.x maintenance line) — worth piloting to
improve Rxn-INSIGHT's mappings, but not the lever that unlocks classification (details
below).

## Current state in the ORD codebase

Rxn-INSIGHT is live, in `ord-schema` (not `ord-interface`):

- `ord-schema/ord_schema/orm/reaction_class.py` calls
  `Reaction(smiles, rxn_mapper=RXNMapper()).get_reaction_info()` and returns a
  `(reaction_class, reaction_name)` tuple, e.g. `("C-C Coupling", "Suzuki coupling with
  boronic acids")`.
- Results land in a dedicated `derived.reaction_classes` table (PK `reaction_id`; two
  indexed `Text` columns, `reaction_class` / `reaction_name`), populated during the
  derived post-pass via the optional `--classify_reactions` flag. Rxn-INSIGHT's
  `"OtherReaction"` sentinel is normalized to NULL.
- We already depend on **rxnmapper** transitively, via a `skearnes/Rxn-INSIGHT` fork
  pinned for modern NumPy/SciPy, with `transformers>=4.53,<5` and `setuptools<81` pins
  to keep it importable.
- `ord-interface` and `ord-app` depend on `ord-schema[orm]` but **not**
  `[reaction-class]`, so there's no downstream coupling. Classification is capped at 4
  worker shards because each worker loads a transformer model.

**Structural point that makes this cheap:** the ORD proto has no RXNO/NameRxn field.
Classification is purely derived data, so adding RXNO is a derived-table change, not a
schema change.

## What Rxn-INSIGHT actually is

Correcting a common misattribution: it's **Dobbelaere, Lengyel, Stevens & Van Geem
(Ghent University), *J. Cheminformatics* 16:37 (2024), MIT-licensed** — not a
Schwaller/EPFL project. Its only tie to Schwaller is that it *depends on* rxnmapper. It
is **rule-based, not ML**:

- **Class**: 10 superclasses + "Miscellaneous", assigned by bond-electron-matrix (BE)
  predicates evaluated in a fixed dispatch order (first match wins). The 10 come from
  **Carey et al. 2006** — the same source as NameRxn's level-1 superclasses.
- **Name**: deterministic matching against **528 curated SMIRKS** templates; the first
  template whose predicted product matches the actual product wins, else
  `"OtherReaction"`. No confidence score.
- **Output**: `get_reaction_info()` returns an 18-key dict (`CLASS`, `NAME`,
  `MAPPED_REACTION`, functional groups, rings, scaffold, tag, etc.).
- **Coverage ceiling (on USPTO)**: ~90% get a superclass, only **~51% get a specific
  name**, ~10% fall into Miscellaneous. Entirely dependent on atom-mapping quality — the
  authors flag that BE-matrix classification "becomes inaccurate when the atom mapping
  fails."

The key implication: Rxn-INSIGHT already emits **named reactions** whose vocabulary
overlaps heavily with RXNO's ~500 named reactions, and a superclass scheme that shares
RXNO/NameRxn's lineage. That's exactly what makes a crosswalk viable.

## The RXNO / NameRxn / IBM RXN landscape

- **RXNO (Name Reaction Ontology)** — RSC-maintained (Colin Batchelor), an OBO Foundry
  ontology, **CC-BY 4.0**, distributed as OWL + OBO. IDs like `RXNO:0000006`
  (Diels-Alder). ~500+ reaction classes (OLS reports ~1,019 terms including imports from
  BFO/CHEBI/MOP). It is **named-reaction-dominated**; generic mechanistic transformations
  are delegated to a companion ontology (MOP). **Ships no classifier** — it is a
  controlled vocabulary only. Current loaded release is a somewhat stale 2021-12-16.
- **NameRxn (NextMove Software)** — proprietary rule-based expert system. Given a
  reaction it emits a three-level `superclass.category.name` code **and the matching RXNO
  ID directly** (documented example: `3.1.1 Bromo Suzuki coupling → RXNO:0000140`).
  ~1000+ leaf classes; ships in the HazELNut suite. It **is** the crosswalk the RSC
  ontology relies on in practice.
- **Schneider / USPTO-50k** — the widely used reaction-class benchmarks originate from
  Schneider, Lowe, Sayle & Landrum (*JCIM* 2015) and are **labeled by NameRxn**:
  Schneider-50k carries 50 real NameRxn leaf codes, while the retrosynthesis USPTO-50k
  keeps only the 10 superclasses (table below). The follow-up section draws the (often
  muddled) distinction.
- **IBM RXN's open model = `rxnfp`** (Schwaller et al., *Nat. Mach. Intell.* 3, 2021,
  MIT) — a BERT reaction-fingerprint classifier, ~98% class accuracy, but trained on
  **NameRxn-labeled** Pistachio, so it's an ML surrogate for a proprietary scheme. Its
  one fully open label set, USPTO 1k TPL, uses **template-hash** labels, not RXNO.
  [DRFP](https://github.com/reymond-group/drfp) reaches ~99% of rxnfp's performance
  without atom mapping.

The 10 top-level classes (USPTO-50k counts):

| # | Superclass | Count |
| --- | --- | --- |
| 1 | Heteroatom alkylation and arylation | 15,204 |
| 2 | Acylation and related processes | 11,972 |
| 3 | C–C bond formation | 5,667 |
| 4 | Heterocycle formation | 909 |
| 5 | Protections | 672 |
| 6 | Deprotections | 8,405 |
| 7 | Reductions | 4,642 |
| 8 | Oxidations | 822 |
| 9 | Functional group interconversion (FGI) | 1,858 |
| 10 | Functional group addition (FGA) | 231 |

**Bottom line for RXNO:** there is no open `(SMILES → RXNO ID)` dataset. Only the
proprietary NameRxn produces RXNO IDs; open models predict into NameRxn-derived or
template-hash label spaces, neither of which is RXNO.

## Atom mapping: two different "new" rxnmappers

Two things both get called "the new rxnmapper"; they are different.

**`rxn4chemistry/rxnmapper` 0.4.3 (2026-02-13), MIT — maintenance, same model.** The
0.4.x line is the **same ALBERT model** (`albert_heads_8_uspto_all_1310k`), not a retrain:

| Version | Date | Notable content |
| --- | --- | --- |
| 0.3.1 | 2024-08-27 | Relaxed torch requirement |
| 0.4.0 | 2024-09-19 | Atom-placeholder (`*`) support; better long-reaction errors |
| 0.4.1 | 2025-03-13 | Dropped Python 3.7 from CI |
| 0.4.2 | 2025-06-01 | Removed the `torch<2.1` cap |
| 0.4.3 | 2026-02-13 | Pinned `transformers>=4.0,<5`; install-doc fixes |

Output unchanged — `{mapped_rxn, confidence}` dicts, `BatchedMapper` for bulk. Worth
upgrading only because older versions don't install cleanly on Python 3.11/3.12 (breakage
in issues #54/#67/#70 — the reason our fork carries pins); it won't improve mappings.

**`yvsgrndjn/RXNMapper_v2` — the actual retrain (what the follow-up pointed to).** A
genuinely reengineered mapper from the **Reymond group (Yves Grandjean)**, coauthors incl.
Schwaller and Genheden (AstraZeneca), MIT. Default model `alberta-uspto-2800k` (~2.8M
USPTO reactions, on HuggingFace); it also bundles the original 2021 model as a selectable
option. It **keeps the original API** (`from rxnmapper import RXNMapper`,
`get_attention_guided_atom_maps`), so it's a drop-in for Rxn-INSIGHT's `RXNMapper()` call
site. Friction:

- **Not on PyPI** — source install only. Dist name is `rxnmapper-v2` but it **imports as
  `rxnmapper`**, so it collides with the pinned `rxnmapper` package; you install it
  *instead of* it, not alongside.
- **Heavier, tighter pins:** `transformers==4.46.3` (exact), `torch>=2.4`, `numpy<1.24`,
  `rdkit>=2024.3.5` — expect resolver work against the fork's environment.
- **No published benchmarks yet.** The cited 2026 ChemRxiv preprint doesn't resolve and
  the README has no accuracy tables. Actively maintained (pushed 2026-06-29).

**Verdict:** RXNMapper_v2 targets Rxn-INSIGHT's documented weak spot (BE-matrix
classification degrades when the mapping fails), so better mappings *plausibly* lift
classification — but the gain is unquantified and the packaging/namespace friction is
real. **Pilot it in an isolated env; don't swap the production dependency until the
preprint lands with numbers (and ideally a PyPI release).** Either way, atom mapping is an
input-quality lever, not what sets the label taxonomy.

## Open NameRxn classification (follow-up 2026-07-02)

The pivotal follow-up question — *can I get NameRxn classes openly, and skip RXNO?* — has
a positive answer, with a hard granularity ceiling.

**What's public, and at what granularity:**

- **Schneider-50k** (`schneider50k.tsv`, in the rxnfp + DRFP repos, MIT) — **50 real
  NameRxn `N.N.N` leaf codes**, 1,000 reactions each (50k total, verified from the file),
  spanning 9 of the 10 superclasses and 28 categories. Human-readable names come from the
  sibling `rxnclass2name.json` (e.g. `3.1.1` = "Bromo Suzuki coupling", `10.1.5` =
  "Wohl-Ziegler bromination"); `rxnclass2id.json` maps the 50 codes to 0–49 for ML. So
  you get **code + name + RXNO-mappable label**, openly.
- **USPTO-50k** (retrosynthesis; retrosim/GLN, MIT) — only the **10 superclass integers**.
  Correction to a common claim: Schneider-50k and this USPTO-50k are **not the same
  reactions relabeled** — they're ~7%-overlapping samples from the same USPTO/Lowe corpus
  under the same NameRxn taxonomy. Use Schneider-50k when you want leaf codes.
- **USPTO-1k-TPL** — 1,000 classes, but **reaction-template hashes, not NameRxn**. Easy to
  mistake for "1000 NameRxn classes"; it isn't.
- **Pistachio / NameRxn** — the full **~967 NameRxn leaf classes** live only here,
  proprietary. No open slice at that granularity.
- **ORD itself** carries no structured NameRxn field (only a free-text `REACTION_TYPE`
  identifier), so there's nothing in our own data to bootstrap-train from.

**Open models that emit these labels:**

| Model | Emits | Classes | License | Notes |
| --- | --- | --- | --- | --- |
| **rxnfp** (Schneider classifier) | real NameRxn `N.N.N` | 50 | MIT | pretrained, ~98%, run today |
| **rxnfp** (1k TPL) | template hashes | 1000 | MIT | not NameRxn |
| **DRFP** + your MLP | whatever you train on | 50 (Schneider) | MIT | fingerprint only, 0.956 on Schneider |
| **SynCat** (2026) | NameRxn (Schneider) | 50 | code TBC; paper CC-BY-NC | GINE GNN, ~0.988 SOTA, atom-mapping-free |
| `pingzhili` ChemBERTa (HF) | 10 superclasses | 10 | MIT | 87%, reactant/product input |

**SynCat** ([phuocchung123/SynCat](https://github.com/phuocchung123/SynCat)) is worth a
closer look: it's **atom-mapping-free** (no rxnmapper dependency), and its edge over rxnfp
*grows with granularity* — on 680-class reaction-center clusters it hits 0.982 vs rxnfp's
0.802. It also ships several open finer-grained label schemes beyond the 50 NameRxn codes
— SynTemp reaction-center clusters (143 / 356 / 680 classes) and a 63-class mechanistic
set — though those are cluster/mechanism IDs, not NameRxn/RXNO codes. Flag: the paper is
**CC-BY-NC** (non-commercial); confirm the code/data license before relying on it in ORD.

**The near-NameRxn open outlier — Schwaller ReactionClassifier (2026).** arXiv 2607.01061,
MIT, `schwallergroup/ReactionClassifier`; authors incl. **Maarten Dobbelaere** (the
Rxn-INSIGHT author). An agentic-LLM (Gemini 3) pipeline that writes *verifiable*
generalized SMIRKS rules and self-expands a 68-class seed into **14,073 classes over 7
levels** (19 superclasses, 106 subclasses, ~1,546 L3 types). Deterministic verification: a
rule's template must reproduce the recorded product before a label is emitted. Reported
~97.7% on unseen reactions — matching NameRxn — and expert chemists preferred its labels
over NameRxn (82.6% vs 61.8%). Ships a ~666k-reaction labeled USPTO DB (NameRxn columns
excluded as proprietary). Caveats: its codes are its **own** taxonomy (not NameRxn/RXNO
IDs), OOD coverage on 2025 academic reactions trails NameRxn (~68% vs 89%), and it's a
1-day-old, non-peer-reviewed preprint; the LLM layer needs a Gemini key (the deterministic
classifier + released rules run offline).

**So, does rxnfp "spit out NameRxn classes"?** Its *Schneider* classifier does — real
`N.N.N` codes at 50-leaf granularity. Its 1k-TPL classifier does not. And there are
*enough public labels to train* an open NameRxn classifier only up to ~50 leaf classes;
beyond that the proprietary ceiling holds — unless you accept a self-generated taxonomy
like the ReactionClassifier's.

## Fallback: RXNO crosswalk of Rxn-INSIGHT names

Kept for the case where you retain Rxn-INSIGHT for its human-readable named reactions and
want RXNO IDs bolted on. Now **dominated by the Schneider-50k path** for emitting
machine-readable codes (that path yields real NameRxn codes → RXNO directly, no curation),
but this remains low-effort and additive if you keep Rxn-INSIGHT:

1. **Curation table** — a versioned data file mapping `rxn_insight_name → RXNO:id` (the
   528 SMIRKS names) plus `superclass → RXNO:id` (10 entries). Because both taxonomies
   trace to Carey 2006 and RXNO's ~500 named reactions overlap Rxn-INSIGHT's 528, most
   entries are direct matches. An LLM-assisted first pass against the RXNO OWL, plus
   chemist review, is roughly a day or two. This is the real work — but it's bounded and
   reviewable, not open-ended.
2. **Schema** — add a nullable `reaction_rxno_id` (and optionally a cached `rxno_label`)
   to `derived.reaction_classes`. No proto change.
3. **Population** — extend `classify_reaction_smiles` to look up the RXNO ID from the
   name it already computes. Pure dictionary lookup; no new model call.
4. **Display** — embed RXNO's CC-BY labels/definitions in the UI; link out to OLS.
5. **Honest coverage note** — named reactions get a leaf RXNO ID (~51%), classed-only
   reactions get a coarse superclass RXNO node (~39%), Miscellaneous gets nothing
   (~10%). Document this so downstream users don't overread the labels.

## Next steps

- **Prototype the Schneider-50k path**: run pretrained rxnfp (or train DRFP/SynCat) over a
  sample of our reactions, write the 50 `N.N.N` codes (+ names via `rxnclass2name.json`, +
  RXNO IDs) into a new `derived.reaction_classes` column, and compare labels against the
  current Rxn-INSIGHT output.
- **Evaluate the Schwaller ReactionClassifier** on a batch of ORD reactions — check
  coverage/quality vs Rxn-INSIGHT and whether the finer taxonomy justifies the Gemini
  dependency. If it holds up, it could replace both Rxn-INSIGHT and the crosswalk.
- **Pick the granularity target** (10 superclasses → 50 leaves → ~1,500 → full ~967
  NameRxn); that choice selects the path. Only the last needs a NameRxn license.
- Pilot **RXNMapper_v2** in an isolated env to see if better mappings lift Rxn-INSIGHT
  classification; keep it out of the production dependency until it has published numbers
  and a PyPI release.
- Fold the `rxn4chemistry/rxnmapper` 0.4.3 pin bump into the fork as routine hygiene.

## For the next session (handoff)

**Where things live.** All repos are siblings under `~/github/ord/`: `ord-schema`
(classification code + ORM), `ord-interface` (search API/UI + local-DB tooling), `ord-app`
(frontend), `ord-data` (datasets), `ord-infrastructure` (Pulumi). Python runs in the
**`ord` conda env** (`~/mambaforge/envs/ord`); its `bin/initdb` and Postgres binaries drive
the local DB, *not* the system PATH.

**How classification runs today.** Install the extra into the env:
`pip install ord-schema[reaction-class]` (pulls the `skearnes/Rxn-INSIGHT` fork pinned at
commit `eb71946`, plus rxnmapper). Then
`python -m ord_schema.orm.scripts.add_datasets --classify_reactions ...` populates
`derived.reaction_classes`. The classifier entry points are `classify_reaction_smiles` and
`update_reaction_classes` in `ord-schema/ord_schema/orm/reaction_class.py`. A local/test
Postgres is built via `setup_test_postgres`
(`ord-interface/ord_interface/client/build_database.py`).

**Open question before committing to option 1 (rxnfp path).** Verify whether the rxnfp
repo ships a *ready-to-run Schneider classifier checkpoint* or only reaction fingerprints
you must fit a head on. The ~98% figure quoted in the landscape is the *Pistachio* number;
confirm the *Schneider-50k* accuracy and the actual load/predict API
(`rxn4chemistry/rxnfp`, Zenodo weights) before scoping the work. DRFP is the fallback
(fingerprint + your own MLP — definitely train-your-own).

**RXNO-mapping caveat — don't assume it's free.** Public Schneider data has the `N.N.N`
codes and names but **not** RXNO IDs. Turning a code into `RXNO:xxxxxxx` needs a crosswalk,
most cheaply by matching the NameRxn class *name* (from `rxnclass2name.json`) against the
RXNO OWL — a small, bounded lookup step, not automatic. Skip it entirely if you only need
NameRxn codes and don't actually need RXNO IDs.

**Decision still open (owner: Steven).** The granularity target — 10 superclasses / 50
NameRxn leaves / ~1,500 (ReactionClassifier) / full ~967 (license NameRxn) — is unmade,
and that choice selects the path. Nothing below the full-967 tier needs a paid license.

## References

- Rxn-INSIGHT: paper <https://pmc.ncbi.nlm.nih.gov/articles/PMC10980627/>, repo
  <https://github.com/mrodobbe/Rxn-INSIGHT>
- RXNO: OBO Foundry <https://obofoundry.org/ontology/rxno.html>, repo
  <https://github.com/rsc-ontologies/rxno>, browser <https://www.ebi.ac.uk/ols4/ontologies/rxno>
- NameRxn (NextMove): <https://www.nextmovesoftware.com/namerxn.html>
- Schneider et al., *JCIM* 2015 (50k scheme, NameRxn-labeled):
  <https://nextmovesoftware.com/blog/2015/02/05/paper-on-reaction-fingerprints-now-out/>
- rxnfp: repo <https://github.com/rxn4chemistry/rxnfp>, paper
  <https://doi.org/10.1038/s42256-020-00284-w>; 50k results
  <https://rxn4chemistry.github.io/rxnfp/results_classification_50k/>
- Schneider-50k data (real NameRxn codes): `schneider50k.tsv` +
  `rxnclass2name.json` <https://github.com/rxn4chemistry/rxnfp/tree/master/data>;
  paper <https://doi.org/10.1021/ci5006614>
- USPTO-50k (10 superclasses, retrosynthesis): retrosim
  <https://github.com/connorcoley/retrosim>, GLN <https://github.com/Hanjun-Dai/GLN>
- DRFP: repo <https://github.com/reymond-group/drfp>, paper
  <https://doi.org/10.1039/D1DD00006C>
- SynCat (2026): <https://doi.org/10.1039/D5DD00367A>
- HuggingFace 10-class model:
  <https://huggingface.co/pingzhili/chemberta-v2-finetuned-uspto-50k-classification>
- Schwaller ReactionClassifier (2026): paper <https://arxiv.org/abs/2607.01061>, repo
  <https://github.com/schwallergroup/ReactionClassifier>
- NameRxn / Pistachio (proprietary, ~967 classes):
  <https://www.nextmovesoftware.com/namerxn.html>,
  <https://www.nextmovesoftware.com/pistachio.html>
- rxnmapper: repo <https://github.com/rxn4chemistry/rxnmapper>, PyPI
  <https://pypi.org/project/rxnmapper/>, paper <https://doi.org/10.1126/sciadv.abe4166>
- RXNMapper_v2: repo <https://github.com/yvsgrndjn/RXNMapper_v2>, model
  <https://huggingface.co/yvsgrndjn/alberta-uspto-2800k>
- Current ORD usage: `ord-schema/ord_schema/orm/reaction_class.py`,
  `ord-schema/ord_schema/orm/derived_mappers.py` (`ReactionClasses`)
