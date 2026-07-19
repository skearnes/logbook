# Hydrogen model specification (draft)

A worked example of the "specify, then implement" approach to a clean-room
cheminformatics core. The goal is not to describe what RDKit does, but to state
what a hydrogen model *should* be, with RDKit's current behavior as evidence for
which distinctions are real and which are incidental.

Status: draft for discussion, revision 2 (see §8). Line references are to the
RDKit checkout as of 2026-07-18.

## 1. Why hydrogens first

Hydrogen handling is the smallest problem in the core that still exhibits every
pathology worth fixing:

- State that is only valid after an unenforced preprocessing step.
- One concept spread across multiple fields.
- A transform (`AddHs`/`RemoveHs`) that is documented as an inverse pair but is
  not actually invertible.
- A default argument that silently returns a wrong answer.

If the method works here it should generalize to sanitization ordering, ring
info, and stereo perception. If it does not, we learn that cheaply.

## 2. What RDKit does today

An atom's hydrogen state is spread across four independent pieces:

| State | Where | Meaning |
| --- | --- | --- |
| `d_numExplicitHs` | `Atom.h:408` | H count asserted on the atom, not in the graph |
| `d_implicitValence` | `Atom.cpp:298` | Cached inferred count; `-1` until valence is calculated |
| `df_noImplicit` | `Atom.h:407` | Suppresses inference entirely |
| H atoms as graph nodes | the molecular graph | Hydrogens materialized as real atoms |

`getTotalNumHs()` sums only the first two:

```cpp
// Atom.cpp:287
unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
  int res = getNumExplicitHs() + getNumImplicitHs();
  if (includeNeighbors && dp_mol) { /* count H neighbors */ }
  return res;
}
```

`includeNeighbors` defaults to `false`. `AddHs` zeroes the explicit count when
it materializes hydrogens (`AddHs.cpp:608`). Therefore, after `AddHs`, a methyl
carbon reports `getTotalNumHs() == 0` — the default reading of "how many
hydrogens does this atom have" is wrong precisely when hydrogens are most
present.

`getNumImplicitHs()` additionally carries a runtime precondition:

```cpp
// Atom.cpp:303
PRECONDITION(d_implicitValence > -1,
             "getNumImplicitHs() called without preceding call to "
             "calcImplicitValence()");
```

This is a runtime assertion for "you forgot to sanitize."

### 2.1 What the parsers reveal

The parsers show that `numExplicitHs` and `noImplicit` are always set together,
and always mean the same thing — *this count came from the input; do not infer
it*:

- SMILES bracket atoms: `smiles.yy:379`, `smiles.yy:385`, and the `H_TOKEN`
  productions at `smiles.yy:402-408`.
- Molfile `HCOUNT` fields: `MolFileParser.cpp:824`, `:3090`, `:3201`, each
  paired with a `setNoImplicit(true)`.
- SMARTS: `smarts.yy:50-56`, `:684-695` — same pairing, but here the count is a
  *query*, and "unspecified" is a meaningful third state distinct from zero.

So the two fields encode one concept with three values, not two concepts.

Note the SMARTS row: it is the same syntax carrying a fundamentally different
meaning. That divergence drives the central design decision in §3.

### 2.2 What `RemoveHsParameters` reveals

`MolOps.h:302-331` defines fifteen boolean knobs governing which hydrogens
`RemoveHs` will not remove: `removeDegreeZero`, `removeHigherDegrees`,
`removeOnlyHNeighbors`, `removeIsotopes`, `removeAndTrackIsotopes`,
`removeDummyNeighbors`, `removeDefiningBondStereo`, `removeWithWedgedBond`,
`removeWithQuery`, `removeMapped`, `removeInSGroups`, `removeHydrides`,
`removeNontetrahedralNeighbors`, and two deprecated aliases.

Read as a whole, these are not fifteen independent policies. They are an ad-hoc
enumeration of a single predicate: *does this hydrogen carry information beyond
its own existence?* An isotope, an atom map, a wedge bond, a query, an SGroup
role, a second bond — each is a piece of data that a bare integer count cannot
hold.

The presence of this list is the strongest available evidence that
count-representation is the wrong default: every entry is a case where the count
representation lost information and a flag was added to prevent it.

## 3. The model

### 3.1 Hydrogens are always graph nodes

Every hydrogen is a node in the molecular graph. There is no count
representation, no materialize/dematerialize transform, and no state a molecule
can be in with respect to hydrogens.

What varies is not the molecule but the **iteration view**, chosen at each call
site:

```rust
mol.atoms()        // all atoms, hydrogens included
mol.heavy_atoms()  // hydrogens skipped
mol.bonds()        // all bonds
mol.heavy_bonds()  // bonds between heavy atoms
```

This is the load-bearing decision in the spec, and it settles the "where do
properties live" question with a single answer: **on the hydrogen node, always**.
Isotope, atom map, charge, wedge role, SGroup membership, and query attachment
have exactly one home, because the node they belong to always exists.

Everything in §2.2 dissolves. There is no lossy transform, so there is nothing
for the fifteen flags to protect against, and the `_isotopicHs` side-channel
property that RDKit uses to smuggle isotope labels across `RemoveHs`/`AddHs` has
no reason to exist.

The alternative — treating materialization as a molecule state with a partial
dematerialize operation — was revision 1 of this spec and is strictly worse; see
§8.

### 3.2 Total hydrogen count

Every atom has a total hydrogen count. It is derived, not stored:

```text
total_h(atom) = |{n in neighbors(atom) : n.atomic_num == 1}|
```

Always defined, never cached, never dependent on a preceding call. This single
derivation replaces `getNumExplicitHs`, `getNumImplicitHs`, and both readings of
`getTotalNumHs`, and it eliminates the `Atom.cpp:303` precondition by
construction rather than by discipline.

### 3.3 Provenance splits in two

Provenance — did this count come from the input, or was it inferred? — does not
disappear, but it stops being one concept. It splits along the `Mol`/`Query`
boundary of §3.4:

- **On `Mol`, provenance is a serialization concern.** `[CH3]` and `C` produce
  identical graphs; the difference is only that the first must be written back
  as `[CH3]`. This is an annotation consulted by writers, not part of the
  chemical model. That is a demotion, and a welcome one.
- **On `Query`, the distinction is core.** "Unspecified H count" versus
  "constrained to 3" is a matching-semantics difference, not a formatting one.

RDKit's `numExplicitHs`/`noImplicit` pair was carrying both of these at once,
which is why it reads as incoherent.

There is one residual case: an *asserted* count on `Mol` also means "do not
recompute this if the atom is edited." Under always-explicit that concern
largely evaporates, because H counts never change implicitly at all (§5).

### 3.4 Queries are a separate type

A SMARTS pattern cannot use this representation. `[CH3]` in SMARTS is a
constraint over hydrogen count, not an assertion of it, and a constraint cannot
be materialized as three nodes. `[C]` — no hydrogen constraint at all — has no
materialized form whatsoever.

So the type hierarchy forks:

- `Mol` — concrete, always-explicit hydrogens, no query state.
- `Query` — hydrogen counts as constraints; may also contain explicit `[H]`
  nodes where the pattern genuinely calls for one.

Substructure matching is then `match(query: &Query, target: &Mol)`, operating
across two representations by design rather than by accident.

Critically, **queries and transforms are still authored without explicit
hydrogens**. `[CH3]` stays `[CH3]`; nobody writes out three `[H]` nodes. The
pattern matches over the target's heavy view, evaluating `H` against the derived
`total_h` of §3.2. Matches therefore return heavy-atom handles only, with the
attached hydrogens reachable via `atom.hydrogens()` when a caller wants them.

The always-explicit representation is thus invisible to rule authors. This is
what makes the fork affordable: MolStandardize's 145 tautomer transforms and 68
normalizations are authored against implicit-H semantics and stay that way (but
see §3.7 for the transform-side caveat).

The fork is still a genuine improvement on its own merits: RDKit stores query
atoms inside `RWMol`, so every consumer of a molecule must consider whether it
might be holding a query, and that ambiguity leaks throughout the codebase.

### 3.4.1 The heavy-atom property set

For queries to be authored implicit-H, every property they can test must be
available on a heavy atom without the caller seeing hydrogen nodes. RDKit
already defines this set — and it already contains both views:

| SMARTS | RDKit implementation | View |
| --- | --- | --- |
| `D` | `getDegree()` (`QueryOps.h:83`) | graph neighbors |
| `X` | `getTotalDegree()` (`QueryOps.h:86`) | graph neighbors + counts |
| `h` | `getTotalNumHs(false)` (`QueryOps.h:118`) | counts only |
| `H` | `getTotalNumHs(true)` (`QueryOps.h:115`) | graph neighbors + counts |

These primitives are defined against *storage state* rather than chemistry.
Measured on RDKit 2026.03.4 against `CC(=O)O`, before and after `AddHs`:

| Pattern | Implicit H | After `AddHs` | |
| --- | --- | --- | --- |
| `[CD1]` | `((0,),)` | `()` | changed |
| `[Ch3]` | `((0,),)` | `()` | changed |
| `[CX4]` | `((0,),)` | `((0,),)` | stable |
| `[CH3]` | `((0,),)` | `((0,),)` | stable |
| `[CH0]` | `((1,),)` | `((1,),)` | stable |
| `[Cv4]` | `((0,), (1,))` | `((0,), (1,))` | stable |

The same pattern gives different answers for the same molecule depending on an
unrelated preprocessing call. `H`, `X`, and `v` are representation-independent;
`D` and `h` are not.

**This is not an RDKit choice.** The Daylight SMARTS specification defines these
primitives in exactly these terms:

- `H<n>` — "total-H-count", "*n* attached hydrogens"
- `h<n>` — "implicit-H-count", "*n* implicit hydrogens"
- `D<n>` — "degree", "*n* explicit connections"
- `X<n>` — "connectivity", "*n* total connections"

RDKit implements the spec faithfully. The dual representation is baked into the
*query language standard*: `h` and `D` only have meaning in a model where a
hydrogen may be either a graph node or a count. The `[CD1]` instability above is
spec-conformant behavior, not a defect in RDKit.

Note also the terminology inversion: in Daylight's vocabulary a bracket hydrogen
such as `[CH3]` is *implicit* — it is not a graph node. RDKit stores that count
in a field named `numExplicitHs`, which means the opposite of the specification's
word for it. This is an independent argument for §3.3's reading that the field
is really tracking *assertion*, not materialization.

### 3.4.2 Divergence from the Daylight specification

Because the dual representation is in the standard, always-explicit cannot adopt
Daylight SMARTS unchanged. Two primitives must be redefined, and this is a
deliberate incompatibility rather than a cleanup:

- **`D` is rebound to heavy degree.** Read literally, "explicit connections"
  under always-explicit would count hydrogen nodes, making `D` equivalent to `X`
  and breaking essentially every `D` pattern ever written. Binding it to
  `heavy_degree` instead preserves the *intent* of existing patterns — an author
  writing `[CD1]` means methyl-like — and makes the primitive stable across
  preprocessing. It is nonetheless a semantic change to a standardized language.
- **`h` is deleted**, or aliased to `H`. Under one representation there are no
  implicit hydrogens for it to count.

Usage across RDKit's own shipped catalogs (1717 non-comment pattern lines in
`Data/` and the MolStandardize catalogs) sizes the two decisions very
differently:

| Primitive | Occurrences | Notes |
| --- | --- | --- |
| `D` | 170 | 37 of 38 lines in `FunctionalGroups.txt`; 84 of 116 in `patty_rules.txt` |
| `h` | 0 | unused throughout |

Deleting `h` is therefore close to free, at least by this sample. `D` is heavily
used — but every one of these catalogs is authored against implicit-hydrogen
molecules, which is how RDKit is used in practice, so rebinding `D` to heavy
degree *preserves all 170*. The rebinding is not a risk to existing patterns; it
is what keeps them working. Only patterns authored against explicit-hydrogen
molecules would invert, and there appear to be none in RDKit's own corpus.

The trade is defensible: the spec-conformant behavior is itself the footgun, and
the redefinition makes patterns mean what their authors meant. But it must be
documented as intentional divergence from a published standard, with a migration
note, rather than presented as a free improvement. Note also that RDKit's
shipped catalogs are one curated corpus, not a representative sample of
real-world SMARTS; the survey in §7 still stands.

This also reclassifies the hydrogen model in the project's own triage. It is not
purely a Tier 1 internal wart with identical observable behavior; it is partly a
Tier 2 decision, because SMARTS matching semantics change in a way callers can
observe.

The remaining properties are stable under always-explicit and need no
redefinition:

- `total_h(atom)` — count of hydrogen neighbors. Backs SMARTS `H`.
- `heavy_degree(atom)` — count of non-hydrogen neighbors. Backs SMARTS `D`,
  which now means the same thing permanently.
- `total_degree(atom)` — all neighbors; equals `heavy_degree + total_h`. Backs
  SMARTS `X`.
- `total_valence(atom)` — summed bond orders over all neighbors. Backs `v`.
- `heavy_valence(atom)` — summed bond orders to heavy neighbors.
- Ring membership and ring-connection counts, evaluated over the heavy view.

SMARTS `h` — "implicit hydrogen count" — has no referent in this model and
should be deleted. It is the one primitive that exists purely as an artifact of
the dual representation.

### 3.4.3 Transforms need a hydrogen reconciliation policy

Queries only read, so §3.4 fully covers them. Transforms write, and that is
where always-explicit imposes a real obligation.

A rule authored implicit-H — a 1,3 keto-enol shift, say — changes heavy-atom
bond orders and expects hydrogen counts to follow from valence. Against an
always-explicit graph, some hydrogen node must physically move from carbon to
oxygen. The engine therefore has to:

1. Match over the heavy view.
2. Apply the authored connectivity and bond-order changes to heavy atoms.
3. Recompute the hydrogen count each affected atom should now carry.
4. Reconcile: add, remove, or move actual hydrogen nodes to match.

Steps 3–4 reintroduce valence inference — but *scoped to the transform engine*
rather than pervading the data model, which is the right place for it.

Step 4 has a genuine ambiguity: **which hydrogen moves?** When the candidates are
chemically equivalent it does not matter. When one is deuterium it matters a
great deal, and a rule authored in implicit-H terms has no way to express the
choice. This model does not create that problem — RDKit has it too, and handles
it with the `_isotopicHs` side channel — but it does force it into the open,
where it must be settled as a stated policy (for example: prefer moving a
default-isotope hydrogen; move an isotope-labeled one only when it is the sole
candidate). An explicit policy is an improvement over emergent behavior, but it
is unavoidable work rather than a free win.

### 3.4.4 Rejected: a per-node `is_implicit` flag

An obvious middle path is to keep every hydrogen in the graph — so properties
have a home — but tag each hydrogen node `is_implicit` or `is_explicit`. Its
attraction is that it restores referents for Daylight `h` and `D` and removes
the divergence of §3.4.2.

It should be rejected. The governing test:

> A per-node annotation is acceptable if and only if nothing but the serializer
> reads it. Once matching reads it, the dual representation is back.

The flag fails that test by construction. For `h` and `D` to work, the flag must
be live during matching, so every molecule carries an implicit/explicit
assignment that changes match results — and `[CD1]` again answers differently
depending on how the flags happen to be set. The flag does not fix the `D`
instability; it re-implements it, buying specification compliance by preserving
the defect.

The rest of §2.2 returns with it. `RemoveHs` becomes "flip explicit to
implicit," which revives the question *may this hydrogen be flipped?* — an
isotope-labeled hydrogen marked implicit cannot be written as `[CH3]`, so the
writer must override the flag. Advisory state silently overridden by a
downstream layer is the same footgun family as `getTotalNumHs`'s default
argument. And algorithms face three views (all hydrogens, explicit-flagged,
implicit-flagged) rather than two.

### 3.4.5 What the flag was right about: write-time derivation

The underlying intuition is nonetheless correct: something per-hydrogen *is*
needed for faithful round-tripping, and the per-heavy-atom provenance of §3.3
does not supply it. Writing `[2H]C` requires folding three hydrogens into a
count while materializing the fourth — a per-node decision that provenance
cannot express.

That decision is **derived at write time, not stored**. A hydrogen is written as
a node if and only if it carries information that count notation cannot hold:

- a non-default isotope;
- a non-zero formal charge;
- an atom map number;
- degree other than one;
- a stereo or wedge role on its bond;
- membership in a substance group.

This is revision 1's collapsibility predicate, relocated from the data model to
the serializer. The relocation is what makes it correct. As a graph operation it
was lossy and induced molecule states; as a writer policy it destroys nothing,
tracks no partial state, and expresses precisely what belongs in a serializer —
what the target format can represent.

One case resists derivation: "this molfile drew its hydrogens explicitly and
byte-fidelity is wanted." That is a genuine stored annotation, but it is
per-file in practice and, critically, is still read only by the writer. It
passes the test.

### 3.4.6 Daylight compatibility as a matching dialect

Rather than diverging from the specification unconditionally (§3.4.2), make it
selectable. The in-memory model stays pure; a matcher running in Daylight
dialect derives an implicit/explicit assignment at match time using the §3.4.5
predicate and evaluates `h` and `D` against it.

```rust
match(query, target, SmartsDialect::Daylight)  // h, D per the specification
match(query, target, SmartsDialect::Modern)    // h absent, D is heavy degree
```

Specification-exact semantics are then available on request, the data model
stays uncontaminated, and legacy pattern files — the case that actually wants
legacy semantics — get them. The divergence becomes a default rather than an
absolute.

### 3.5 Representation

The obvious objection to always-explicit is cost. The NCI 5K sample has a median
of 15 heavy atoms and a mean of 16.4; drug-like hydrogen-to-heavy ratios run
roughly 1:1 to 1.3:1, so this roughly doubles atom count and more than doubles
bond count. At screening-library scale that is not free.

Traversal cost, however, is recoverable by **partitioning storage**:

- Atoms are stored heavy-first: heavy atoms occupy `[0, n_heavy)`, hydrogens
  occupy `[n_heavy, n_atoms)`.
- Neighbor lists are likewise ordered heavy-first.

`heavy_atoms()` is then a contiguous slice scan with no per-atom branch and no
cache pollution from hydrogen records, and `heavy_neighbors()` is a prefix slice
of the neighbor list. The common case — heavy-only traversal, which is most
algorithms — costs what the count representation costs. The remaining penalty is
memory, not time.

The price is that atom indices no longer correspond to input file order and are
not stable across edits. RDKit consumers depend on both. The intended resolution
is to make indices opaque handles and expose input order as an explicit
property, which is better design but a real migration hazard: it breaks
downstream code silently rather than loudly. See §6.

### 3.5.1 Coordinates are where the model actually strains

Chemistry is not where always-explicit costs the most; geometry is. If every
hydrogen is a node, every conformer needs a coordinate for every hydrogen. A
molfile carrying 3D heavy-atom coordinates and no hydrogens then forces an
unpleasant choice: fabricate hydrogen positions on load, or admit a
"no coordinate" sentinel — and a sentinel makes conformers partial,
reintroducing exactly the partial state §3.1 removed. Two-dimensional depiction
has the same problem from the other direction, since it deliberately does not
place most hydrogens.

Heavy-first partitioning resolves this. With heavy atoms occupying
`[0, n_heavy)`, a conformer is a **prefix array** over the heavy atoms, with an
optional hydrogen extension covering `[n_heavy, n_atoms)`:

```text
conformer.heavy   -> coordinates for [0, n_heavy)        always present
conformer.all     -> coordinates for [0, n_atoms)        present iff extended
```

No sentinel, no fabricated geometry, and no partial state — a conformer either
has the hydrogen extension or it does not, which is a property of the conformer
rather than of individual atoms.

This is the second thing partitioning buys, after contiguous heavy iteration
(§3.5). Two independent problems resolved by one storage decision is reason to
treat heavy-first ordering as load-bearing in the design rather than as a
performance optimization to be applied later.

### 3.6 What the type system still carries

Because hydrogens are no longer a molecule state, the typestate encoding gets
simpler than revision 1 assumed. There is no `HExplicit`/`HImplicit` parameter:

```rust
fn embed_conformer(mol: &Mol<Sanitized>) -> Conformer;
fn morgan_fingerprint(mol: &Mol<Sanitized>) -> Fingerprint;
```

Sanitization state remains worth encoding, for the reason `Atom.cpp:303` exists:
"did you run the required preprocessing pass" is a static property currently
checked at runtime, and typestate makes that class of bug a compile error at no
runtime cost. This remains the strongest argument for a clean-room
implementation in an expressive type system — stronger than memory safety.

## 4. Worked examples

| Input | Graph | Provenance | Notes |
| --- | --- | --- | --- |
| `C` (SMILES) | C + 4 H nodes | inferred | writes back as `C` |
| `[CH3-]` | C(-1) + 3 H nodes | asserted | writes back bracketed |
| Molfile C, no `HCOUNT` | C + 4 H nodes | inferred | |
| Molfile C, `HCOUNT=3` | C + 3 H nodes | asserted | |
| `[2H]C` | C + 3 H + 1 D node | inferred | isotope lives on the node |
| `[CH3]` in SMARTS | `Query`, H-count constraint = 3 | n/a | not materialized |
| `[C]` in SMARTS | `Query`, no H constraint | n/a | no materialized form |

The `[2H]C` row is the case RDKit handles with `removeIsotopes`,
`removeAndTrackIsotopes`, and the `_isotopicHs` side-channel property. Here it
needs no special handling at all: the deuterium is a node, and its isotope label
lives on it permanently.

## 5. What this buys

- One derived accessor instead of three, with no correctness-relevant default
  argument and no cache.
- `noImplicit` disappears; provenance replaces it, split cleanly across
  `Mol` (serialization) and `Query` (matching semantics).
- The fifteen removal flags disappear rather than being derived.
- The `_isotopicHs` side channel disappears.
- Every property has exactly one home, because the node always exists.
- Queries and transforms are still authored implicit-H, so the SMIRKS catalogs
  port unchanged (modulo the reconciliation policy of §3.4.3).
- SMARTS `D` becomes representation-independent: `[CD1]` matches a methyl carbon
  unconditionally, rather than silently ceasing to match after `AddHs`. Note
  this is a deliberate divergence from the Daylight specification, not a bug
  fix — see §3.4.2.
- Stereo perception simplifies: a tetrahedral center always has four real
  neighbors and a wedge bond always has a real endpoint, so the phantom-neighbor
  special-casing threaded through `Chirality.cpp` (4k lines) should shrink.
  Unmeasured — treat as a hypothesis to check, not a claim.
- Hydrogen counts never change implicitly. In RDKit, editing a molecule silently
  recomputes them on the next sanitize; here, adding or removing hydrogen is a
  deliberate graph edit.
- "Forgot to sanitize" becomes a compile error.

## 6. Open questions

- **Atom index stability.** Heavy-first partitioning breaks input-order
  correspondence and edit stability. Opaque handles plus an explicit input-order
  property is the intended answer, but the failure mode for ported code is
  silent. Needs a migration story before adoption.
- **How much does `Query` duplicate?** If `Query` and `Mol` share little, the
  fork is cheap; if matching, traversal, and serialization all need two
  implementations, it is expensive. Worth prototyping before committing.
- **How much existing SMARTS does the `D` redefinition break?** The rebinding in
  §3.4.2 preserves intent for patterns meaning "methyl-like", but any pattern
  authored against an explicit-H molecule means the opposite. Needs a survey of
  real-world pattern corpora, not just RDKit's shipped catalogs.
- **The hydrogen reconciliation policy (§3.4.3).** Which hydrogen node moves
  when an implicit-H-authored transform changes a heavy atom's hydrogen count?
  Irrelevant for equivalent hydrogens, decisive when one is isotope-labeled.
  Needs to be stated, and the stated policy needs checking against RDKit's
  emergent behavior on the tautomer catalog to see how often they diverge.
- **Memory at scale.** Doubling atom count matters for billion-compound
  libraries. Does a compressed on-disk form that inflates on load suffice, or
  does the in-memory representation itself need a count-based variant for
  screening workloads? Note that reintroducing one would bring back most of §2.2.
- **Conformer hydrogen extensions (§3.5.1).** Is "has hydrogen coordinates" a
  property of the conformer, or does a molecule with several conformers need
  them to agree? Disagreement is representable and probably should be rejected,
  but that is an invariant to state rather than discover.
- **Query hydrogen constraints.** SMARTS allows ranges and negation, so the
  constraint type is richer than an optional integer.
- **Valence model coupling.** Inferred counts depend on the valence model, which
  is itself unspecified in RDKit. This spec assumes one exists; that is the next
  document, and it is the immediate blocker.

## 7. Validation plan

The spec is only worth anything if it can be checked against reality:

1. Implement `total_h` as a *Python* shim over RDKit first, before writing any
   Rust: add hydrogens, count H neighbors.
2. Run it over ChEMBL. For every atom, compare to
   `getTotalNumHs(includeNeighbors=True)` on the original molecule.
3. Every disagreement is a finding. Classify each as (a) spec is wrong, (b)
   RDKit is wrong, or (c) legitimate design divergence.
4. Round-trip every molfile in the corpus and confirm asserted counts survive
   `AddHs` and writing.
5. Enumerate the cases where RDKit's default `RemoveHs` declines to remove a
   hydrogen, and confirm each corresponds to information that the always-explicit
   model simply retains.
6. Measure the memory and traversal cost of always-explicit with heavy-first
   partitioning against RDKit on a fingerprinting benchmark, to test the §3.5
   claim that traversal cost is recoverable.
7. **Done** (§3.4.1, §3.4.2). `D` and `h` confirmed unstable across `AddHs` on
   RDKit 2026.03.4; `H`, `X`, and `v` stable. Shipped catalogs use `D` 170 times
   and `h` zero times. Remaining: survey real-world pattern corpora outside
   RDKit, which are the population the `D` rebinding actually risks.
8. Replay the tautomer catalog under a stated reconciliation policy (§3.4.2) and
   diff against RDKit's output, to find how often the isotope ambiguity is
   actually reached.

Steps 1–5 are deliberate: they cost days in Python, require no Rust, and a high
disagreement rate is the signal to stop before committing to an implementation.

## 8. Design history

**Revision 1** modeled materialization as a *state of the molecule* — a molecule
could be materialized, dematerialized, or mixed — with `dematerialize` as a
partial operation. That required a "collapsibility predicate" stating when a
hydrogen may be reduced to a count (no isotope, no atom map, no query, degree
one, no stereo or wedge role, no SGroup membership), from which RDKit's fifteen
`RemoveHsParameters` flags could be derived.

Deriving fifteen flags from one rule was an improvement over enumerating them,
but revision 1 asserted that "materialization is a view, not a property of the
molecule" and then immediately contradicted itself by making it a molecule state
with a lossy transform between values.

**Revision 2** takes the original claim seriously: hydrogens are always nodes,
and the view is a property of *iteration* rather than of the molecule. The
collapsibility predicate is then unnecessary rather than merely well-factored,
and the fifteen flags do not need deriving because there is no lossy transform
for them to guard. The cost is that queries no longer fit the representation,
forcing the `Mol`/`Query` split of §3.4 — which is arguably correct on its own
merits, but is now a load-bearing part of the design rather than an incidental
cleanup.

Revision 2 initially treated that split as expensive, on the assumption that
SMIRKS catalogs would need rewriting against explicit-H targets. That was wrong:
the view abstraction extends to *authoring*, not just iteration, so patterns stay
implicit-H and match against the target's heavy view (§3.4). What survives of
the objection is narrower and lives in transforms rather than queries — the
reconciliation policy of §3.4.2.

Investigating that led to the §3.4.1 finding: SMARTS *already* carries both
views as separate primitives (`D`/`X`, `h`/`H`), but binds two of them to
storage state rather than to chemistry.

That was initially recorded as evidence the design is right — the language
wanting a model the data model failed to supply. Checking the Daylight
specification corrected it. `h` ("implicit-H-count") and `D` ("*n* explicit
connections") are Daylight primitives, defined in those terms by the standard;
RDKit implements them faithfully. The dual representation is therefore in the
*query language standard*, not merely in RDKit's data model.

This makes the finding sharper but less comfortable. Always-explicit is not
cleaning up an implementation artifact — it requires deliberately diverging from
a published specification on two primitives (§3.4.2), which moves the hydrogen
model out of Tier 1 (invisible internal change) and partly into Tier 2
(observable semantic change with a compatibility cost). The trade still looks
right, because the spec-conformant behavior is the footgun, but it must be
argued rather than assumed.
