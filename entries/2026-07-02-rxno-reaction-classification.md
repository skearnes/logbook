# 2026-07-02 — Reaction classification against the RXNO ontology

## Question

The ORD pipeline currently uses [Rxn-INSIGHT](https://github.com/mrodobbe/Rxn-INSIGHT)
for reaction classification and naming. I'd prefer output that maps onto the **RXNO**
(Name Reaction Ontology) — or NameRxn, but that's proprietary. What would it take to
build a classification system/model that emits RXNO categories? This entry scopes the
landscape: how Rxn-INSIGHT works today, what RXNO/NameRxn/IBM RXN actually are, whether
the newer rxnmapper helps, and the realistic options.

## Summary

**RXNO is a vocabulary, not a classifier.** No open tool reads a reaction and emits an
RXNO ID; the only turnkey engine that does is NextMove's **proprietary NameRxn** (it
outputs an RXNO ID alongside its own `N.N.N` code). Every large RXNO/NameRxn-labeled
corpus in existence (Pistachio, the Schneider-50k set) was labeled by NameRxn. There is
**no open `(reaction SMILES → RXNO ID)` dataset** to train on.

So "output RXNO categories" has three honest shapes, and the recommended one is the
middle:

1. **Crosswalk Rxn-INSIGHT names → RXNO IDs (recommended, open).** We already run
   Rxn-INSIGHT and store `reaction_name`; RXNO is CC-BY and freely embeddable; and — the
   load-bearing fact — Rxn-INSIGHT's taxonomy and RXNO/NameRxn's superclasses **both
   descend from Carey et al. 2006**, so the alignment is real, not forced. This is an
   additive derived-data change, no new model or infra.
2. **License NameRxn (turnkey, proprietary).** Highest accuracy and coverage, hands us
   RXNO IDs directly. Costs money and raises a label-redistribution question for an open
   database.
3. **Train an open ML classifier to emit RXNO (not worth it).** Blocked by the missing
   open RXNO-labeled dataset; you'd need NameRxn to generate labels first (circular).

**The trigger that flips me from option 1 to option 2:** if we need leaf-level RXNO
granularity at high coverage. NameRxn names essentially everything into 1000+ classes;
Rxn-INSIGHT assigns a *specific* name to only ~51% of reactions. If coarse-but-open is
acceptable, the crosswalk wins.

The rxnmapper upgrade is worth taking for dependency hygiene but is **not** the lever
that unlocks RXNO (details below).

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
- **Schneider / USPTO-50k** — the widely used 10-superclass scheme originates from
  Schneider, Lowe, Sayle & Landrum, *JCIM* 2015, **labeled by NameRxn**. Its 10
  superclasses *are* NameRxn's level-1 classes (see table).
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

## rxnmapper — the new version

Latest is **0.4.3 (2026-02-13)**, MIT. Reality check: the 0.4.x line is the **same
ALBERT model** (`albert_heads_8_uspto_all_1310k`) — compatibility/maintenance work, not a
retrain.

| Version | Date | Notable content |
| --- | --- | --- |
| 0.3.1 | 2024-08-27 | Relaxed torch requirement |
| 0.4.0 | 2024-09-19 | Atom-placeholder (`*`) support; better long-reaction errors |
| 0.4.1 | 2025-03-13 | Dropped Python 3.7 from CI |
| 0.4.2 | 2025-06-01 | Removed the `torch<2.1` cap |
| 0.4.3 | 2026-02-13 | Pinned `transformers>=4.0,<5`; install-doc fixes |

Output is unchanged — a list of `{mapped_rxn, confidence}` dicts, with `BatchedMapper`
for bulk use. Atom mapping quality matters for classification because a correct mapping
defines the **reaction center**, which drives template/SMARTS extraction (as in RDChiral)
— but note that rxnfp and DRFP classify **without** atom mapping at all.

**Verdict:** worth upgrading because older versions don't install cleanly on Python
3.11/3.12 (documented breakage in issues #54/#67/#70, which is exactly why our fork
carries those pins), but it won't materially improve mappings. Treat it as routine
dependency hygiene, decoupled from the RXNO decision.

## What the recommended path (option 1) involves

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

- Decide option 1 vs. option 2 based on whether leaf-level coverage is a hard
  requirement.
- If option 1: draft the RXNO crosswalk starter table by matching Rxn-INSIGHT's 528
  SMIRKS names against the RXNO OWL, then get it chemist-reviewed.
- Independently, take the rxnmapper 0.4.3 bump into the fork's pins as dependency
  hygiene.

## References

- Rxn-INSIGHT: paper <https://pmc.ncbi.nlm.nih.gov/articles/PMC10980627/>, repo
  <https://github.com/mrodobbe/Rxn-INSIGHT>
- RXNO: OBO Foundry <https://obofoundry.org/ontology/rxno.html>, repo
  <https://github.com/rsc-ontologies/rxno>, browser <https://www.ebi.ac.uk/ols4/ontologies/rxno>
- NameRxn (NextMove): <https://www.nextmovesoftware.com/namerxn.html>
- Schneider et al., *JCIM* 2015 (50k scheme, NameRxn-labeled):
  <https://nextmovesoftware.com/blog/2015/02/05/paper-on-reaction-fingerprints-now-out/>
- rxnfp: repo <https://github.com/rxn4chemistry/rxnfp>, paper
  <https://doi.org/10.1038/s42256-020-00284-w>; DRFP <https://github.com/reymond-group/drfp>
- rxnmapper: repo <https://github.com/rxn4chemistry/rxnmapper>, PyPI
  <https://pypi.org/project/rxnmapper/>, paper <https://doi.org/10.1126/sciadv.abe4166>
- Current ORD usage: `ord-schema/ord_schema/orm/reaction_class.py`,
  `ord-schema/ord_schema/orm/derived_mappers.py` (`ReactionClasses`)
