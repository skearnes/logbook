# 2026-07-18 — Reimplementing RDKit in Rust: scoping and a specs-first pivot

## Summary

Scoping conversation on what it would take to reimplement RDKit in Rust. The
question moved through three framings — full parity, a useful subset, and
specs-first clean-room design — and only the third is worth pursuing. Includes a
worked example spec for hydrogen handling
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
  SMILES string is not a molecule until sanitization runs, which pulls in the
  43k lines at GraphMol's root: `Chirality.cpp` alone is 4k, plus ~10k lines of
  chirality *tests*, which is a direct readout of the edge-case count. Canonical
  SMILES is worse: "correct" means byte-identical to RDKit, because those
  strings are database keys, and a ranking difference silently breaks joins
  rather than erroring.
- **Tautomers/protonation (`MolStandardize` 6.3k + ~1k lines of rule catalogs) —
  smallest and hardest.** The rules are SMIRKS, so porting them requires a
  SMARTS matcher with RDKit-identical aromaticity-aware matching, which lands
  back in sanitization. And the output is heuristic by construction: tautomer
  canonicalization is enumerate-then-score with a hand-tuned function, and
  `acid_base_pairs.in` is ordered such that the order *is* the behavior.

The lesson: these targets do not decompose along the API surface. They share one
substrate — sanitized-graph semantics plus SMARTS matching — and that substrate
holds nearly all of the unspecified behavior. The niche code being cut is the
part that was cleanly separable; the core is not. Realistic scope ~55–70k lines,
3–5 person-years, mostly differential-testing archaeology.

## Framing 3: specs-first — the one worth doing

Dropping the parity constraint changes the economics, because RDKit's release
cadence stops mattering. Two things make this better than a port:

- **The spec has value independent of any implementation.** A rigorous account
  of the hydrogen model, sanitization ordering, and one aromaticity model is a
  contribution whether or not Rust ships — a far better risk profile than a
  rewrite that is worthless at 70%.
- **It fixes the epistemics of testing.** Differential-testing against RDKit
  alone makes RDKit definitionally correct, so "our bug" and "their bug" are
  indistinguishable. With a spec, a disagreement is classifiable: spec says X,
  RDKit does Y, we do Z. Disagreements become information.

### Triage: which areas admit specs

- **Tier 1 — internal design warts. Spec pays off fully, observable behavior can
  stay identical.** Hydrogen model, sanitization as an explicit ordered
  pipeline, ring-info caching and invalidation, valence/charge bookkeeping, and
  the untyped property bag on atoms (where much of the cross-module coupling
  hides).
- **Tier 2 — genuine choices. Spec-able, but the spec is a decision that costs
  compatibility.** Aromaticity is the type case: RDKit has four models because
  there is no right answer. Pick one, define it rigorously, make alternatives
  explicit plugins. Same for stereo perception and canonical ranking. The cost
  is that canonical output diverges, breaking anything using it as a key.
- **Tier 3 — resists specification.** Tautomer scoring, protonation rules,
  normalization transforms. These encode empirical judgment; no spec discipline
  produces a right answer. Best move is to make the heuristic *data* —
  versioned, declarative rule sets with explicit evaluation and scoring
  semantics. RDKit half-does this with the `.in` catalogs; the improvement is
  making ordering and scoring first-class rather than build-time data files with
  undocumented order dependence.

Note this inverts the Framing 2 advice: tautomers stay least tractable, but
because there is nothing to discover, not because of archaeology.

### The Rust argument, restated

The strongest case for Rust is not memory safety. It is that a large fraction of
RDKit footguns have the shape "did you run the right preprocessing pass?" —
`getNumImplicitHs()` carries a literal runtime assertion for "you forgot to
sanitize" (`Atom.cpp:303`). Typestate makes that entire bug class
unrepresentable at zero runtime cost. This only becomes available once different
abstractions are permitted.

## The worked example: hydrogens

Full draft in [HYDROGEN_MODEL_SPEC.md](assets/HYDROGEN_MODEL_SPEC.md). Headline
results:

- RDKit spreads hydrogen state across four independent pieces
  (`d_numExplicitHs`, `d_implicitValence`, `df_noImplicit`, and H atoms as graph
  nodes).
- `getTotalNumHs()` sums only the first two and defaults `includeNeighbors` to
  `false`, while `AddHs` zeroes the explicit count (`AddHs.cpp:608`). So after
  `AddHs`, a methyl carbon reports zero hydrogens by default.
- The original hunch was "three models should be one." That is nearly right, but
  `numExplicitHs` is load-bearing: it encodes **provenance** (a count *asserted*
  by the input versus one *inferred* from valence rules), which matters for
  lossless round-tripping and for the query case where "unspecified" differs
  from zero. Confirmed in the parsers — `smiles.yy:379`, `MolFileParser.cpp:824`
  and `:3090`, `smarts.yy:50-56` all set the count and the flag together.
- So the clean model is one invariant (total H count) plus two orthogonal bits:
  provenance (asserted/inferred/unspecified) and materialization (node vs
  count).
- Best find: the fifteen boolean knobs in `RemoveHsParameters`
  (`MolOps.h:302-331`) are not fifteen policies. They are an ad-hoc enumeration
  of one predicate — *does this hydrogen carry information beyond its own
  existence?* Stating that predicate once derives all fifteen, and makes the
  `_isotopicHs` side-channel property unnecessary.

The general lesson, which is the main reason to keep going: a meaningful
fraction of what looks like cruft turns out to be load-bearing for a case that
would otherwise break. Finding out which is which *is* the work, and it is not
parallelizable across people.

## Cost estimate

Writing real specs for hydrogens, sanitization ordering, aromaticity, and stereo
is itself a research project — roughly 6–12 months of one experienced person,
most of it spent reading RDKit to learn *why* each wart exists.

## Next steps

- Validate the hydrogen spec in Python over ChEMBL before writing any Rust
  (§7 of the spec). Days of work, no Rust required, and a high disagreement rate
  is the signal to stop early.
- If it converges, write the valence-model spec next — the hydrogen spec assumes
  one exists, and it is the immediate blocker.
- Resolve the typestate encoding question (§6) before it metastasizes across
  every signature; two type parameters is already awkward.
- Survey prior art properly, especially Richard Apodaca's `chemcore`/`purr`
  work, specifically for *where it stopped* — that is an empirical cost
  estimate.
- Independently: the `getTotalNumHs()` default-argument footgun is worth
  reporting upstream regardless of whether any of this proceeds. It fits the
  additive/documentary fixes recommended in the
  [2026-07-15 API ergonomics plan](../2026-07-15-rdkit-maintenance-plans/assets/API_ERGONOMICS_PLAN.md).
