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
- The fifteen boolean knobs in `RemoveHsParameters` (`MolOps.h:302-331`) are not
  fifteen policies. They are an ad-hoc enumeration of one predicate — *does this
  hydrogen carry information beyond its own existence?* Each flag marks a case
  where the count representation lost information and a flag was bolted on to
  prevent it, which is the strongest available evidence that counts are the
  wrong default.

### Revision 2: always-explicit hydrogens

The spec went through two revisions in one sitting, and the second is the
interesting one.

Revision 1 modeled materialization as a *state of the molecule* (materialized /
dematerialized / mixed) with a partial `dematerialize` operation, and derived
the fifteen flags from a single "collapsibility predicate." Better than
enumerating them — but it asserted that materialization is a view rather than a
property of the molecule, then immediately contradicted itself by making it a
molecule state with a lossy transform.

Revision 2 takes the claim seriously: **hydrogens are always graph nodes, and
the implicit/explicit distinction is a property of *iteration*** —
`mol.atoms()` versus `mol.heavy_atoms()`, chosen per call site. Consequences:

- "Where do properties live" gets one answer: on the hydrogen node, always.
  Isotope, atom map, charge, wedge role, SGroup membership all have exactly one
  home, because the node always exists.
- The collapsibility predicate becomes *unnecessary* rather than
  well-factored. The fifteen flags do not need deriving; there is no lossy
  transform for them to guard. `_isotopicHs` goes away for the same reason.
- Provenance survives, but splits: on `Mol` it is a serialization annotation
  (`[CH3]` and `C` give identical graphs), while on `Query` the
  specified/unspecified distinction is core matching semantics. RDKit's
  `numExplicitHs`/`noImplicit` pair was carrying both at once, which is why it
  reads as incoherent.
- The typestate problem partly answers itself: no `HExplicit`/`HImplicit`
  parameter is needed, just `Mol<Sanitized>`.
- Stereo should simplify — every stereocenter has four real neighbors, every
  wedge bond a real endpoint, so the phantom-neighbor special-casing in
  `Chirality.cpp` should shrink. Unmeasured; recorded as a hypothesis.

Costs, in order of severity: **queries cannot use this representation** (a
SMARTS `[CH3]` is a constraint, not an assertion; `[C]` has no materialized
form), forcing a `Mol`/`Query` type split — arguably correct anyway, since
RDKit stores query atoms inside `RWMol` and that ambiguity leaks everywhere.
Then memory, roughly 2x atoms and more than 2x bonds — traversal cost is
recoverable by storing atoms heavy-first so `heavy_atoms()` is a contiguous
slice, but that breaks atom-index stability, which downstream code depends on
and which fails silently.

### The view abstraction extends to authoring

The initial read was that the `Mol`/`Query` split would be expensive because the
SMIRKS catalogs assume implicit-H semantics and would need rewriting. That was
wrong. Views apply to *authoring* as well as iteration: patterns stay written as
`[CH3]`, matching against the target's heavy view, so the 145 tautomer
transforms port unchanged. Matches return heavy-atom handles, with hydrogens
reachable on demand.

This requires the heavy-atom property set to be complete — `total_h`,
`heavy_degree`, `total_degree`, valences, ring properties — which is what makes
the representation invisible to rule authors.

**The dual representation is in the SMARTS standard, not just in RDKit.**
SMARTS already carries both views as distinct primitives (`QueryOps.h:83-120`) —

| SMARTS | Daylight definition | RDKit implementation | View |
| --- | --- | --- | --- |
| `D` | "*n* explicit connections" | `getDegree()` | graph neighbors |
| `X` | "*n* total connections" | `getTotalDegree()` | graph + counts |
| `h` | "*n* implicit hydrogens" | `getTotalNumHs(false)` | counts only |
| `H` | "*n* attached hydrogens" | `getTotalNumHs(true)` | graph + counts |

`getDegree()` counts real graph neighbors, so `[CD1]` matches a methyl carbon
while hydrogens are implicit and silently stops matching after `AddHs`. `H` and
`X` are representation-independent; `D` and `h` are not.

The first read was that this is an RDKit artifact — the query language wanting a
model the data model failed to supply. **Checking the Daylight spec corrected
that.** `h` ("implicit-H-count") and `D` ("explicit connections") are Daylight
primitives defined in exactly those terms, and RDKit implements them faithfully.
The `[CD1]` instability is spec-conformant, not a defect.

So always-explicit is not a cleanup here. It requires deliberately diverging
from a published standard on two primitives: `D` must be rebound to heavy degree
(read literally it would become equivalent to `X`, breaking every `D` pattern
ever written), and `h` must be deleted or aliased to `H` for lack of a referent.
Rebinding `D` preserves what authors *mean* by `[CD1]` and makes it stable, so
the trade looks right — but it is a semantic change to a standardized language
and has to be argued, with a migration note, rather than sold as a free win.

**This reclassifies the hydrogen model from Tier 1 to partly Tier 2** in the
triage above: not an invisible internal change, but an observable semantic one
with a real compatibility cost.

A smaller find in the same vein: in Daylight's vocabulary a bracket hydrogen
like `[CH3]` is *implicit* — not a graph node. RDKit stores that count in a
field called `numExplicitHs`, which means the opposite of the standard's word
for it. Probably a real source of confusion, and independent support for reading
that field as tracking *assertion* rather than materialization.

**Measured** on RDKit 2026.03.4 against `CC(=O)O`, before and after `AddHs`:
`[CD1]` and `[Ch3]` both change their answer; `[CX4]`, `[CH3]`, `[CH0]`, and
`[Cv4]` do not. Across 1717 non-comment pattern lines in RDKit's shipped
catalogs, `D` appears 170 times (37 of 38 lines in `FunctionalGroups.txt`, 84 of
116 in `patty_rules.txt`) and `h` appears zero times.

So the two decisions are sized very differently. Deleting `h` is close to free.
`D` is pervasive — but those catalogs are all authored against implicit-H
molecules, so rebinding `D` to heavy degree *preserves all 170*. The rebinding
is what keeps them working, not a threat to them. Only patterns authored against
explicit-H molecules invert, and RDKit's corpus contains none. Caveat: one
curated corpus is not a representative sample of real-world SMARTS.

### Rejected: a per-node is_implicit flag

The natural middle path — keep every H in the graph so properties have a home,
but tag each node implicit or explicit — restores referents for Daylight `h` and
`D` and would remove the divergence entirely. It should still be rejected. The
governing test:

> A per-node annotation is acceptable if and only if nothing but the serializer
> reads it. Once matching reads it, the dual representation is back.

The flag fails by construction: for `h` and `D` to work it must be live during
matching, so `[CD1]` again answers differently depending on how flags happen to
be set. **It does not fix the `D` instability, it re-implements it** — buying
spec compliance by preserving the defect. It also revives the §2.2 problem, since
`RemoveHs` becomes "flip explicit to implicit" and an isotope-labeled H marked
implicit cannot be written as `[CH3]`, forcing the writer to override the flag.
And algorithms face three views instead of two.

What the flag is *right* about: something per-hydrogen genuinely is needed for
round-tripping, and per-heavy-atom provenance can't supply it — writing `[2H]C`
means folding three hydrogens and materializing one. But that is **derived at
write time, not stored**: materialize an H iff it carries information count
notation can't hold (isotope, charge, atom map, degree ≠ 1, stereo/wedge role,
SGroup membership). That is revision 1's collapsibility predicate, relocated
from the data model to the serializer — which is what makes it correct. As a
graph operation it was lossy and induced molecule states; as a writer policy it
destroys nothing and tracks no partial state.

Rather than diverging from Daylight unconditionally, the divergence can also be
made a **matching dialect**: keep the model pure, and let a Daylight-dialect
matcher derive an implicit/explicit assignment at match time via the same
predicate. Spec-exact semantics on request, no flag in the data model, and
legacy pattern files get legacy semantics.

### Where leaning in actually strains: coordinates

Geometry, not chemistry, is the real cost — though the first analysis of it was
wrong and got corrected.

The initial claim was that conformers should become prefix arrays over heavy
atoms with an optional hydrogen extension, so a hydrogen-less 3D file needs no
fabricated geometry. That invents partial state to solve a mostly nonexistent
problem. **A conformer covers every atom, hydrogens included** — which is what
RDKit already does (`Conformer` is sized to `numAtoms`; verified on 2026.03.4,
where nine-atom explicit-H ethanol gives nine positions and `Compute2DCoords`
assigns real coordinates to the hydrogens), and what the physics forces, since
3D embedding requires hydrogens.

So the cost does not fall on embedding. It falls on **loading**. PDB files,
crystal structures, and vendor SDFs routinely carry 3D heavy-atom coordinates
and no hydrogens — X-ray below ~1.2 Å does not resolve them.

The next draft said the loader should therefore place idealized hydrogens and
record "these were idealized" as provenance. **That was also wrong**, and the
correction is the better idea: coordinates are simply optional per atom.

### Coordinates are optional per atom

```rust
conformer.position(atom) -> Option<Point3>
```

Absence *is* the provenance. Nothing is fabricated at load, and the type forces
every reader to handle the missing case. The place-and-flag design fabricates
geometry and then depends on consumers checking a flag to learn it isn't real —
the same failure shape as `getTotalNumHs`'s default argument, where the obvious
reading silently returns something wrong.

The stated objection to optional positions had been that they make conformers
"partial," reintroducing the partial state the design removed. That was a
category error: the ban on partial state is about partial *representation*
(a half-materialized molecule is ambiguous about what it means), whereas a
missing coordinate is **missing data** — an ordinary thing to state and a
dishonest thing to paper over.

**The decisive argument is that this isn't a hydrogen question.** Experimental
structures lack positions for heavy atoms too — disordered regions, unresolved
side chains. RDKit cannot express that: `Conformer` offers only
`GetAtomPosition` with no notion of absence, so a missing position must be
encoded as a missing *atom*. Verified on 2026.03.4, a PDB lysine with an
unresolved side chain loads as **four atoms** — chemically a lysine, structurally
a fragment. Chemical identity and observational completeness are conflated, and
the molecule is silently wrong rather than explicitly incomplete. `Option<Point3>`
handles hydrogens and heavy atoms uniformly and fixes it.

Consequences: hydrogen placement becomes an explicit opt-in operation rather than
a hidden step in parsing; provenance demotes to annotating deliberate placement;
and completeness gives the typestate machinery a real second use, since force
fields, RMSD, and shape comparison need full geometry
(`conformer.complete() -> Option<CompleteConformer>`).

It also **retracts the strip-after-embed breakage**. The workflow becomes "drop
hydrogen positions, keep hydrogen nodes," which is expressible and recovers most
of the win, since a multi-conformer library stores the graph once and coordinates
N times. The residual cost is one molecule's worth of hydrogen nodes per library
— far smaller than claimed.

One other follow-on stands: **depictions are not conformers.** RDKit stores 2D
layouts as `Conformer` with `is3D=False`, but a depiction is a rendering
artifact, and conflating them is why "must a conformer have hydrogen
coordinates" looked ambiguous — the answer differs for the two things sharing a
type.

This also retracts the earlier "partitioning solves two problems" claim.
Heavy-first storage still buys contiguous heavy iteration; it is not needed for
conformers.

What survives of the objection is narrower and sits in transforms rather than
queries: transforms *write*, so when a rule changes a heavy atom's hydrogen
count, some H node must actually move, and the engine needs a stated policy for
**which one**. Irrelevant for equivalent hydrogens; decisive when one is
deuterium, which an implicit-H-authored rule cannot express. RDKit has the same
ambiguity and hides it in the `_isotopicHs` side channel; this design forces it
into the open as an explicit policy. Real work, but far smaller than rewriting
the catalogs.

The general lesson, and the main reason to keep going: a meaningful fraction of
what looks like cruft turns out to be load-bearing for a case that would
otherwise break — and the reverse also holds, since revision 1's carefully
derived predicate turned out to be solving a problem the design did not need to
have. Finding out which is which *is* the work, and it is not parallelizable
across people.

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
- Prototype the `Mol`/`Query` split to find out how much it duplicates. If
  matching, traversal, and serialization each need two implementations, the fork
  is expensive and the always-explicit design gets materially less attractive.
- ~~Confirm the `[CD1]` instability empirically and count `D`/`h` usage in
  shipped catalogs.~~ Done — see above. Remaining: survey real-world pattern
  corpora outside RDKit, which are the population the `D` rebinding risks.
- State a hydrogen reconciliation policy, replay the tautomer catalog under it,
  and diff against RDKit to see how often the isotope ambiguity is reached.
- Benchmark always-explicit with heavy-first partitioning against RDKit on
  fingerprinting, to test the claim that traversal cost is recoverable and only
  memory is paid.
- Design a migration story for opaque atom indices before adopting heavy-first
  storage; the breakage is silent, which is the worst kind.
- Survey prior art properly, especially Richard Apodaca's `chemcore`/`purr`
  work, specifically for *where it stopped* — that is an empirical cost
  estimate.
- Independently: the `getTotalNumHs()` default-argument footgun is worth
  reporting upstream regardless of whether any of this proceeds. It fits the
  additive/documentary fixes recommended in the
  [2026-07-15 API ergonomics plan](../2026-07-15-rdkit-maintenance-plans/assets/API_ERGONOMICS_PLAN.md).
