# Hydrogen model specification (draft)

A worked example of the "specify, then implement" approach to a clean-room
cheminformatics core. The goal is not to describe what RDKit does, but to state
what a hydrogen model *should* be, with RDKit's current behavior as evidence for
which distinctions are real and which are incidental.

Status: draft for discussion. Line references are to the RDKit checkout as of
2026-07-18.

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

So the two fields encode one concept with three values, not two concepts. That
is the first simplification.

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
hold. This is the second simplification, and it is the more valuable one.

## 3. The model

### 3.1 Invariant

Every atom has a **total hydrogen count**. It is a property of the chemical
graph, always defined, never dependent on a flag or on a preceding
preprocessing call.

```text
total_h(atom) = count of hydrogens bonded to atom,
                whether or not they are materialized as nodes
```

This single accessor replaces `getNumExplicitHs`, `getNumImplicitHs`, and both
readings of `getTotalNumHs`.

### 3.2 Bit one: provenance

Orthogonal to the count is where it came from:

- **Asserted** — the input stated it (`[CH3]`, a molfile `HCOUNT` field). It is
  data, and it must survive round-tripping unchanged.
- **Inferred** — derived from the valence model. It is a computed view and may
  legitimately change if the valence model changes.
- **Unspecified** — valid only in query contexts, where "no constraint on
  hydrogen count" differs from "constrained to zero."

Provenance is what `noImplicit` is really tracking. Naming it directly means the
"three hydrogen models" collapse to one count plus one enum, and it also means a
round-tripped molfile can be checked for fidelity rather than hoped at.

### 3.3 Bit two: materialization

Also orthogonal: is a given hydrogen a node in the graph, or a number on its
neighbor? This is a property of the *representation*, not of the molecule.
Materialization is a view, and the two views describe the same chemistry.

The key point is that this view is **not total in both directions**.
Materialization always succeeds. Dematerialization does not.

### 3.4 The collapsibility predicate

Rather than fifteen flags, state the condition once. A materialized hydrogen is
**collapsible** — representable as an increment to its neighbor's count —
if and only if it carries no information beyond its existence:

- atomic number 1, and no isotope label;
- formal charge zero;
- no atom map number;
- no query attached;
- degree exactly one;
- its single bond is a plain single bond carrying no stereo or wedge role;
- it is not a member of any substance group;
- its neighbor is not a stereocenter whose configuration the hydrogen helps
  define (including non-tetrahedral cases).

Every one of RDKit's fifteen flags is a special case of this predicate.
Deriving them from one stated rule, rather than accumulating them as they were
discovered, is the actual deliverable of this spec.

### 3.5 Laws

```text
materialize(materialize(m))     == materialize(m)          -- idempotent
dematerialize(materialize(m))   == m,  if m is fully dematerialized
materialize(dematerialize(m))   ≅ m,   iff every H in m is collapsible
```

The third law is conditional, and that condition is the honest content of the
model. `dematerialize` is a **partial** operation: it collapses what it can and
leaves the rest as nodes. A molecule may therefore be in a mixed state, and the
type system should say so rather than pretending otherwise.

This is a real departure from RDKit, where `RemoveHs` silently leaves hydrogens
behind and the caller has no type-level indication that it happened.

### 3.6 Encoding the states

The precondition at `Atom.cpp:303` is a runtime check for a static property.
Preprocessing requirements belong in the type:

```rust
// Sketch; the exact encoding is an open question (see §6).
fn embed_conformer(mol: &Mol<Sanitized, HExplicit>) -> Conformer;
fn morgan_fingerprint(mol: &Mol<Sanitized, AnyH>) -> Fingerprint;
```

A caller cannot pass an unsanitized molecule, or one without materialized
hydrogens, to a function that requires them. Satisfying the bound is written
`mol.with_explicit_h()`, which is a no-op when already satisfied. There is no
runtime cost; the entire "did you remember to call X first" bug class becomes a
compile error.

This generalizes well beyond hydrogens, and is a stronger argument for a
clean-room implementation in a language with an expressive type system than
memory safety is.

## 4. Worked examples

| Input | Total H | Provenance | Materialized |
| --- | --- | --- | --- |
| `C` (SMILES) | 4 | inferred | no |
| `[CH3-]` | 3 | asserted | no |
| `[CH4]` | 4 | asserted | no |
| Molfile C, no `HCOUNT` | 4 | inferred | no |
| Molfile C, `HCOUNT=3` | 3 | asserted | no |
| `[CH3]` in SMARTS | 3 | asserted (query) | no |
| `[C]` in SMARTS | — | unspecified | no |
| `C` after `AddHs` | 4 | asserted | yes |
| `[2H]C` after `RemoveHs` | 3 inferred + 1 node | mixed | partial |

The last row is the case RDKit handles with `removeIsotopes` and
`removeAndTrackIsotopes` (plus an `_isotopicHs` side-channel property to restore
the isotope on re-materialization). Under the collapsibility predicate it needs
no special casing: the deuterium is not collapsible because it carries an
isotope label, so it stays a node. No side-channel property, no flag.

## 5. What this buys

- One accessor instead of three, with no correctness-relevant default argument.
- `noImplicit` disappears as a concept; provenance replaces it and is more
  informative.
- Fifteen removal flags become one predicate, with the flags recoverable as a
  policy override for callers who genuinely need non-default behavior.
- Round-trip fidelity becomes checkable: asserted counts must survive
  unchanged, inferred counts may not.
- The `_isotopicHs` side channel goes away.
- "Forgot to sanitize" and "forgot to add hydrogens" become compile errors.

## 6. Open questions

- **Typestate shape.** Two type parameters (`Mol<S, H>`) is already awkward and
  will get worse as more preprocessing states appear (aromaticity perception,
  stereo assignment, ring info). Alternatives: a single state parameter with
  sealed marker traits, or a capability-set encoding. This needs a real design
  pass before it is adopted; getting it wrong makes every signature noisy.
- **Mixed materialization in the type system.** Does `HPartial` need to be a
  distinct state, or is it enough for `AnyH` to be the permissive supertype?
- **Query hydrogen counts.** SMARTS allows constraints richer than a single
  integer (ranges, negation). "Unspecified" may need to generalize to a
  constraint type rather than being an enum variant.
- **Does provenance need to be per-atom?** Probably yes for correct molfile
  round-tripping, but it costs a byte per atom; worth measuring against a
  molecule-level flag.
- **Valence model coupling.** Inferred counts depend on the valence model, which
  is itself unspecified in RDKit. This spec assumes a valence model exists; that
  is the next document.

## 7. Validation plan

The spec is only worth anything if it can be checked against reality:

1. Implement `total_h` and the collapsibility predicate as a *Python* shim over
   RDKit first, before writing any Rust.
2. Run it over ChEMBL. For every atom, compare the spec's answer to
   `getTotalNumHs(includeNeighbors=True)`.
3. Every disagreement is a finding. Classify each as (a) spec is wrong, (b)
   RDKit is wrong, or (c) legitimate design divergence.
4. Round-trip every molfile in the corpus and confirm asserted counts survive.
5. Compare the collapsibility predicate against `RemoveHs` with default
   parameters; each divergence should map to a named flag or be a genuine bug.

Step 1 is deliberate. Validating the spec against RDKit costs days in Python and
requires no Rust at all — and if the disagreement rate is high, that is the
signal to stop before committing to an implementation.
