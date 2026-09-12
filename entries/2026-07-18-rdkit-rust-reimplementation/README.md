# 2026-07-18 — Reimplementing RDKit in Rust: scoping and a specs-first pivot

- **Date:** 2026-07-18
- **Author:** Steven Kearnes
- **Acknowledgments:** Prepared with [Claude Code](https://claude.com/claude-code) (Claude Opus 4.8)
- **Status:** final (scoping; hydrogen spec in progress)
- **Tags:** rdkit, rust, reimplementation, scoping, hydrogen-model,
  specification, sanitization
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

## Summary

Scoping conversation on what it would take to reimplement RDKit in Rust. The
question moved through three framings — full parity, a useful subset, and
specs-first clean-room design — and only the third is worth pursuing. Includes a
worked hydrogen-handling spec
([HYDROGEN_MODEL_SPEC.md](assets/HYDROGEN_MODEL_SPEC.md)) as a test of whether
the method produces anything real.

Follows on from [2026-07-15 — RDKit maintenance planning
docs](../2026-07-15-rdkit-maintenance-plans/README.md), particularly the API
ergonomics findings; this is the "what if we redesigned instead of patching"
branch of that work.

## Framing 1: full parity — no

Measured from the checkout:

| Scope | Lines |
| --- | --- |
| C++ in `Code/` | ~547k |
| — tests | ~180k |
| — Python/Java wrappers | ~119k |
| — actual library logic | ~250k |
| Python in `rdkit/` | ~118k |

Largest modules: `MolDraw2D` 52k, `FileParsers` 51k, `ChemReactions` 23k,
`MolStandardize` 21k, `SmilesParse` 21k.

Roughly 50–100 person-years, against a target that moves three releases a year.
The blocker is not volume, it is that the semantics are undocumented:
sanitization ordering, valence rules, four aromaticity models, and stereo
perception exist as accumulated bug fixes, not as a spec. Parity means
bug-for-bug compatibility, so the real work is reverse-engineering rather than
writing. Verdict: use the existing `rdkit-sys` FFI bindings instead.

## Framing 2: the useful subset — better, but mis-decomposed

Proposed scope was molecule I/O, fingerprints, and tautomers/protonation, on the
grounds that most of RDKit is niche. The tail does drop cleanly (`MolDraw2D`,
`ChemReactions`, `RGroupDecomposition`, `FMCS`, `DistGeom`). But the three
targets differ enormously in difficulty, roughly inverted from their size:

- **Fingerprints (8k lines) — genuinely easy.** Published algorithms, and
  correctness is decidable: differential-test bit-for-bit against RDKit over
  ChEMBL.
- **Molecule I/O (`SmilesParse` 8k) — a trojan horse for the graph core.** A
  SMILES string is not a molecule until sanitization runs, which pulls in the 43k
  lines at GraphMol's root: `Chirality.cpp` alone is 4k, plus ~10k lines of
  chirality *tests*, a direct readout of the edge-case count. Canonical SMILES is
  worse — "correct" means byte-identical to RDKit, because those strings are
  database keys, and a ranking difference silently breaks joins rather than
  erroring.
- **Tautomers/protonation (`MolStandardize` 6.3k + ~1k lines of rule catalogs) —
  smallest and hardest.** The rules are SMIRKS, so porting them requires a SMARTS
  matcher with RDKit-identical aromaticity-aware matching, which lands back in
  sanitization. And the output is heuristic by construction: tautomer
  canonicalization is enumerate-then-score with a hand-tuned function, and
  `acid_base_pairs.in` is ordered such that the order *is* the behavior.

These targets do not decompose along the API surface. They share one substrate —
sanitized-graph semantics plus SMARTS matching — and that substrate holds nearly
all of the unspecified behavior. The niche code being cut is the part that was
cleanly separable; the core is not. Realistic scope ~55–70k lines, 3–5
person-years, mostly differential-testing archaeology.

## Framing 3: specs-first — the one worth doing

Dropping the parity constraint changes the economics, because RDKit's release
cadence stops mattering. Two things make this better than a port:

- **The spec has value independent of any implementation.** A rigorous account of
  the hydrogen model, sanitization ordering, and one aromaticity model is a
  contribution whether or not Rust ships — a far better risk profile than a
  rewrite that is worthless at 70%.
- **It fixes the epistemics of testing.** Differential-testing against RDKit
  alone makes RDKit definitionally correct, so "our bug" and "their bug" are
  indistinguishable. With a spec, a disagreement is classifiable: spec says X,
  RDKit does Y, we do Z. Disagreements become information.

### Triage: which areas admit specs

- **Tier 1 — internal design warts. Spec pays off fully, observable behavior can
  stay identical.** Sanitization as an explicit ordered pipeline, ring-info
  caching and invalidation, valence/charge bookkeeping, and the untyped property
  bag on atoms (where much of the cross-module coupling hides).
- **Tier 2 — genuine choices. Spec-able, but the spec is a decision that costs
  compatibility.** Aromaticity is the type case: RDKit has four models because
  there is no right answer. Pick one, define it rigorously, make alternatives
  explicit plugins. Same for stereo perception and canonical ranking. The cost is
  that canonical output diverges, breaking anything using it as a key.
- **Tier 3 — resists specification.** Tautomer scoring, protonation rules,
  normalization transforms. These encode empirical judgment; no spec discipline
  produces a right answer. Best move is to make the heuristic *data* — versioned,
  declarative rule sets with explicit evaluation and scoring semantics. RDKit
  half-does this with the `.in` catalogs; the improvement is making ordering and
  scoring first-class rather than build-time data files with undocumented order
  dependence.

This inverts the Framing 2 advice: tautomers stay least tractable, but because
there is nothing to discover, not because of archaeology.

The hydrogen model started in Tier 1 and moved partly into Tier 2 once it became
clear the change is observable through SMARTS (below).

### The Rust argument, restated

The strongest case for Rust is not memory safety. It is that a large fraction of
RDKit footguns have the shape "did you run the right preprocessing pass?" —
`getNumImplicitHs()` carries a literal runtime assertion for "you forgot to
sanitize" (`Atom.cpp:303`). Typestate makes that entire bug class unrepresentable
at zero runtime cost. This only becomes available once different abstractions are
permitted.

## The worked example: hydrogens

Full draft in [HYDROGEN_MODEL_SPEC.md](assets/HYDROGEN_MODEL_SPEC.md). The design
that survived:

**Hydrogens are always graph nodes.** The implicit/explicit distinction is a
property of *iteration* — `mol.atoms()` versus `mol.heavy_atoms()`, chosen per
call site — not a state the molecule is in. This settles where properties live
with one answer: on the hydrogen node, always, because it always exists.

**Queries and transforms are still authored implicit-H.** Patterns stay written
as `[CH3]` and match against the target's heavy view, so MolStandardize's 145
tautomer transforms port unchanged. Queries need a separate type, though, since a
SMARTS `[CH3]` is a *constraint* and cannot be materialized as three nodes.

**Suppression is derived at write time.** A hydrogen is written as a node iff it
carries information a count cannot hold — isotope, charge, atom map, degree ≠ 1,
stereo/wedge role, SGroup membership. That predicate is a *floor*, not a rule:
`Minimal`, `All`, `Polar`, and `AsLoaded` become one writer parameter, and
because the floor is always enforced, policy varies style but cannot vary
correctness.

**Coordinates are optional per atom.** `Option<Point3>`, where absence states
that the position is unknown. Nothing is fabricated at load, and the type forces
every reader to handle the missing case.

Findings worth keeping regardless of whether any Rust is written:

- RDKit spreads hydrogen state across four independent pieces, and
  `getTotalNumHs()` sums only two of them with `includeNeighbors` defaulting to
  `false`. Since `AddHs` zeroes the explicit count (`AddHs.cpp:608`), a methyl
  carbon reports **zero hydrogens** after `AddHs`.
- The fifteen knobs in `RemoveHsParameters` are not fifteen policies. Each marks
  a case where the count representation lost information and a flag was bolted on
  — an ad-hoc enumeration of one predicate.
- **The implicit/explicit duality is in the Daylight SMARTS standard**, not just
  in RDKit. `h` is defined as "implicit-H-count" and `D` as "*n* explicit
  connections", so RDKit implements both faithfully and `[CD1]` losing a methyl
  carbon after `AddHs` is spec-conformant, not a defect. Always-explicit
  therefore requires deliberate divergence on two primitives rather than a
  cleanup — which is what moved this from Tier 1 to partly Tier 2.
- Measured across 1717 pattern lines in RDKit's shipped catalogs: `D` appears
  170 times, `h` **zero** times. So deleting `h` is nearly free, and since those
  catalogs are all authored against implicit-H molecules, rebinding `D` to heavy
  degree *preserves all 170* rather than threatening them.
- RDKit cannot express a missing coordinate. `Conformer` offers only
  `GetAtomPosition`, so an unresolved position must be encoded as a missing
  *atom*: a PDB lysine with an unresolved side chain loads as **four atoms** —
  chemically a lysine, structurally a fragment. Chemical identity and
  observational completeness are conflated.
- RDKit declines to decide whether D/T are hydrogens, keeping
  `queryAtomNonHydrogenDegree` and `queryAtomHeavyAtomDegree` side by side with
  the disagreement recorded in comments (`QueryOps.h:89`, `:102`). Matching is
  isotope-blind while Morgan fingerprints are isotope-aware — and both defaults
  are defensible, which is why a global rule is the wrong shape.

## What the process suggests about the method

Nine rounds of review on the *small* problem, chosen because it was tractable
enough to finish. Three positions were retracted (materialization as molecule
state, conformers as heavy-atom prefix arrays, place-and-flag hydrogen
coordinates), each recorded in §8 of the spec.

Two observations that bear on the cost estimate:

- **The errors were design-discipline failures, not chemistry.** Over-applying a
  principle across a domain boundary — banning partial *coordinates* because the
  spec bans partial *representation* — or inventing a mechanism for a property
  the model already expressed. None of them needed RDKit archaeology to catch,
  which was the cost driver the estimate was built around.
- **Retractions were the productive part.** Revision 1's collapsibility predicate
  was deleted as unnecessary, then returned two rounds later as the load-bearing
  primitive — relocated from the data model to the serializer, where it now does
  triple duty (serialization floor, transform reconciliation, isotope handling).
  The predicate was right; its layer was wrong.

A meaningful fraction of what looks like cruft turns out to be load-bearing, and
the reverse also holds. Finding out which is which *is* the work, and it is not
parallelizable across people — though a second reviewer catches self-consistency
errors far faster than the author does.

## Cost estimate

Writing real specs for hydrogens, sanitization ordering, aromaticity, and stereo
is itself a research project — roughly 6–12 months of one experienced person,
most of it spent reading RDKit to learn *why* each wart exists. On the evidence
above that is optimistic, but for a different reason than expected: the
bottleneck looks like self-consistency of the design rather than archaeology.

## Next steps

- Validate the hydrogen spec in Python over ChEMBL before writing any Rust (§7 of
  the spec). Days of work, no Rust required, and a high disagreement rate is the
  signal to stop early.
- If it converges, write the valence-model spec next — the hydrogen spec assumes
  one exists, and it is the immediate blocker.
- Prototype the `Mol`/`Query` split to find out how much it duplicates. If
  matching, traversal, and serialization each need two implementations, the fork
  is expensive and the always-explicit design gets materially less attractive.
- Survey real-world pattern corpora outside RDKit for `D` and `h` usage — the
  population the `D` rebinding actually risks. (The RDKit-internal half is done.)
- State a hydrogen reconciliation policy, replay the tautomer catalog under it,
  and diff against RDKit.
- Benchmark always-explicit with heavy-first partitioning against RDKit on
  fingerprinting, to test the claim that traversal cost is recoverable and only
  memory is paid.
- Design a migration story for opaque atom indices before adopting heavy-first
  storage; the breakage is silent, which is the worst kind.
- Survey prior art properly, especially Richard Apodaca's `chemcore`/`purr` work,
  specifically for *where it stopped* — that is an empirical cost estimate.
- Independently: the `getTotalNumHs()` default-argument footgun is worth
  reporting upstream regardless of whether any of this proceeds. It fits the
  additive/documentary fixes recommended in the
  [2026-07-15 API ergonomics plan](../2026-07-15-rdkit-maintenance-plans/assets/API_ERGONOMICS_PLAN.md).
