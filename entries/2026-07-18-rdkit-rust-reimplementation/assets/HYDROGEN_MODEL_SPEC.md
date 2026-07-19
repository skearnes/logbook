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

This fork is a real cost — every SMARTS-driven subsystem sits on the seam,
including MolStandardize's 145 tautomer transforms and 68 normalizations. But it
is also a genuine improvement: RDKit stores query atoms inside `RWMol`, so every
consumer of a molecule must consider whether it might be holding a query, and
that ambiguity leaks throughout the codebase.

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
- **Porting SMIRKS transforms.** The tautomer and normalization catalogs are
  authored against implicit-H semantics. Rewriting them for always-explicit
  targets is real work with behavioral risk, in exactly the Tier 3 area where no
  spec exists to check against.
- **Memory at scale.** Doubling atom count matters for billion-compound
  libraries. Does a compressed on-disk form that inflates on load suffice, or
  does the in-memory representation itself need a count-based variant for
  screening workloads? Note that reintroducing one would bring back most of §2.2.
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
