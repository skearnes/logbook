# Hydrogen model specification (draft)

- **Date:** 2026-07-18
- **Author:** Steven Kearnes
- **Acknowledgments:** Prepared with [Claude Code](https://claude.com/claude-code) (Claude Opus 4.8)
- **Status:** draft for discussion
- **Tags:** rdkit, rust, hydrogen-model, specification, sanitization,
  clean-room
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

A worked example of the "specify, then implement" approach to a clean-room
cheminformatics core. The goal is not to describe what RDKit does, but to state
what a hydrogen model *should* be, using RDKit's behavior as evidence for which
distinctions are real and which are incidental.

Status: draft for discussion. Design history and rejected alternatives are in §8
and §9. Line references are to the RDKit checkout as of 2026-07-18; measurements
are against RDKit 2026.03.4.

## 1. Why hydrogens first

Hydrogen handling is the smallest problem in the core that still exhibits every
pathology worth fixing:

- State that is only valid after an unenforced preprocessing step.
- One concept spread across multiple fields.
- A transform (`AddHs`/`RemoveHs`) documented as an inverse pair but not actually
  invertible.
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

`includeNeighbors` defaults to `false`, and `AddHs` zeroes the explicit count
when it materializes hydrogens (`AddHs.cpp:608`). After `AddHs`, a methyl carbon
therefore reports `getTotalNumHs() == 0`: the default reading of "how many
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

`numExplicitHs` and `noImplicit` are always set together, and always mean the
same thing — *this count came from the input; do not infer it*:

- SMILES bracket atoms: `smiles.yy:379`, `:385`, and the `H_TOKEN` productions at
  `:402-408`.
- Molfile `HCOUNT` fields: `MolFileParser.cpp:824`, `:3090`, `:3201`, each paired
  with a `setNoImplicit(true)`.
- SMARTS: `smarts.yy:50-56`, `:684-695` — same pairing, but here the count is a
  *query*, and "unspecified" is a meaningful third state distinct from zero.

So the two fields encode one concept with three values, not two concepts. The
SMARTS row carries a different meaning under the same syntax, and that divergence
drives the central decision in §3.4.

### 2.2 What `RemoveHsParameters` reveals

`MolOps.h:302-331` defines fifteen boolean knobs governing which hydrogens
`RemoveHs` will not remove: `removeDegreeZero`, `removeHigherDegrees`,
`removeOnlyHNeighbors`, `removeIsotopes`, `removeAndTrackIsotopes`,
`removeDummyNeighbors`, `removeDefiningBondStereo`, `removeWithWedgedBond`,
`removeWithQuery`, `removeMapped`, `removeInSGroups`, `removeHydrides`,
`removeNontetrahedralNeighbors`, and two deprecated aliases.

These are not fifteen independent policies. They are an ad-hoc enumeration of a
single predicate: *does this hydrogen carry information beyond its own
existence?* An isotope, an atom map, a wedge bond, a query, an SGroup role, a
second bond — each is data that a bare integer count cannot hold.

Every entry is a case where the count representation lost information and a flag
was added to prevent it, which is the strongest available evidence that
count-representation is the wrong default.

## 3. The model

### 3.1 Hydrogens are always graph nodes

Every hydrogen is a node in the molecular graph. There is no count
representation, no materialize/dematerialize transform, and no state a molecule
can be in with respect to hydrogens.

What varies is the **iteration view**, chosen at each call site:

```rust
mol.atoms()        // all atoms, hydrogens included
mol.heavy_atoms()  // hydrogens skipped
mol.bonds()        // all bonds
mol.heavy_bonds()  // bonds between heavy atoms
```

This settles the question of where properties live with a single answer: **on the
hydrogen node, always**. Isotope, atom map, charge, wedge role, SGroup
membership, and query attachment have exactly one home, because the node they
belong to always exists.

Everything in §2.2 dissolves. There is no lossy transform, so nothing for the
fifteen flags to protect against, and the `_isotopicHs` side-channel property
RDKit uses to smuggle isotope labels across `RemoveHs`/`AddHs` has no reason to
exist.

### 3.2 Total hydrogen count

Every atom has a total hydrogen count, derived rather than stored:

```text
total_h(atom) = |{n in neighbors(atom) : n.atomic_num == 1}|
```

Always defined, never cached, never dependent on a preceding call. This single
derivation replaces `getNumExplicitHs`, `getNumImplicitHs`, and both readings of
`getTotalNumHs`, and eliminates the `Atom.cpp:303` precondition by construction
rather than by discipline.

### 3.3 Provenance splits in two

Provenance — did this count come from the input, or was it inferred? — survives,
but stops being one concept. It splits along the `Mol`/`Query` boundary:

- **On `Mol`, provenance is a serialization concern.** `[CH3]` and `C` produce
  identical graphs; the only difference is that the first must be written back
  bracketed. This is an annotation consulted by writers, not part of the chemical
  model.
- **On `Query`, the distinction is core.** "Unspecified H count" versus
  "constrained to 3" is a matching-semantics difference, not a formatting one.

RDKit's `numExplicitHs`/`noImplicit` pair carries both at once, which is why it
reads as incoherent. The naming compounds it: in Daylight's vocabulary a bracket
hydrogen like `[CH3]` is *implicit* — not a graph node — so RDKit's field name
means the opposite of the standard's word for the same thing.

### 3.4 Queries are a separate type

A SMARTS pattern cannot use this representation. `[CH3]` in SMARTS is a
constraint over hydrogen count, not an assertion of it, and a constraint cannot
be materialized as three nodes. `[C]` — no hydrogen constraint at all — has no
materialized form whatsoever.

So the type hierarchy forks:

- `Mol` — concrete, always-explicit hydrogens, no query state.
- `Query` — hydrogen counts as constraints; may also contain explicit `[H]` nodes
  where the pattern genuinely calls for one.

Substructure matching is then `match(query: &Query, target: &Mol)`, operating
across two representations by design rather than by accident.

Queries and transforms are still **authored without explicit hydrogens**. `[CH3]`
stays `[CH3]`; nobody writes out three `[H]` nodes. The pattern matches over the
target's heavy view, evaluating `H` against the derived `total_h` of §3.2, and
matches return heavy-atom handles with attached hydrogens reachable via
`atom.hydrogens()`. The representation is therefore invisible to rule authors,
which is what makes the fork affordable: MolStandardize's 145 tautomer transforms
and 68 normalizations stay as written (subject to §3.4.4).

The fork is also an improvement on its own merits. RDKit stores query atoms
inside `RWMol`, so every consumer of a molecule must consider whether it might be
holding a query, and that ambiguity leaks throughout the codebase.

#### 3.4.1 The heavy-atom property set

For queries to be authored implicit-H, every property they can test must be
available on a heavy atom without the caller seeing hydrogen nodes. RDKit already
defines this set, and it already contains both views:

| SMARTS | Daylight definition | RDKit implementation |
| --- | --- | --- |
| `D` | "*n* explicit connections" | `getDegree()` (`QueryOps.h:83`) |
| `X` | "*n* total connections" | `getTotalDegree()` (`QueryOps.h:86`) |
| `h` | "*n* implicit hydrogens" | `getTotalNumHs(false)` (`QueryOps.h:118`) |
| `H` | "*n* attached hydrogens" | `getTotalNumHs(true)` (`QueryOps.h:115`) |

`D` and `h` are defined against storage state rather than chemistry. Measured
against `CC(=O)O`, before and after `AddHs`:

| Pattern | Implicit H | After `AddHs` | |
| --- | --- | --- | --- |
| `[CD1]` | `((0,),)` | `()` | changed |
| `[Ch3]` | `((0,),)` | `()` | changed |
| `[CX4]` | `((0,),)` | `((0,),)` | stable |
| `[CH3]` | `((0,),)` | `((0,),)` | stable |
| `[CH0]` | `((1,),)` | `((1,),)` | stable |
| `[Cv4]` | `((0,), (1,))` | `((0,), (1,))` | stable |

The same pattern gives different answers for the same molecule depending on an
unrelated preprocessing call.

**This is not an RDKit choice.** The Daylight SMARTS specification defines `h` as
"implicit-H-count" and `D` as "*n* explicit connections"; RDKit implements both
faithfully. The dual representation is baked into the *query language standard* —
`h` and `D` only have meaning in a model where a hydrogen may be either a graph
node or a count. The `[CD1]` instability is spec-conformant behavior, not a
defect.

#### 3.4.2 Divergence from the Daylight specification

Because the dual representation is in the standard, always-explicit cannot adopt
Daylight SMARTS unchanged. Two primitives must be redefined, and this is a
deliberate incompatibility rather than a cleanup:

- **`D` is rebound to heavy degree.** Read literally, "explicit connections"
  under always-explicit would count hydrogen nodes, making `D` equivalent to `X`
  and breaking essentially every `D` pattern ever written. Binding it to
  `heavy_degree` preserves the *intent* of existing patterns — an author writing
  `[CD1]` means methyl-like — and makes the primitive stable across
  preprocessing.
- **`h` is deleted**, or aliased to `H`. Under one representation there are no
  implicit hydrogens for it to count.

Usage across RDKit's shipped catalogs (1717 non-comment pattern lines in `Data/`
and the MolStandardize catalogs) sizes the two decisions very differently:

| Primitive | Occurrences | Notes |
| --- | --- | --- |
| `D` | 170 | 37 of 38 lines in `FunctionalGroups.txt`; 84 of 116 in `patty_rules.txt` |
| `h` | 0 | unused throughout |

Deleting `h` is close to free. `D` is heavily used, but every one of these
catalogs is authored against implicit-hydrogen molecules, so rebinding `D` to
heavy degree *preserves all 170* — the rebinding is what keeps them working
rather than a risk to them. Only patterns authored against explicit-hydrogen
molecules invert, and RDKit's corpus contains none. One curated corpus is not a
representative sample, so the survey in §7 still stands.

The trade is defensible: the spec-conformant behavior is itself the footgun, and
the redefinition makes patterns mean what their authors meant. But it must be
documented as intentional divergence from a published standard, with a migration
note, rather than presented as a free improvement. It also reclassifies the
hydrogen model in the project's triage — not a Tier 1 internal wart with
identical observable behavior, but partly a Tier 2 decision, since SMARTS
matching semantics change observably.

The remaining properties need no redefinition:

- `total_h(atom)` — count of hydrogen neighbors, isotope-blind (§3.6). Backs `H`.
- `heavy_degree(atom)` — count of neighbors with atomic number greater than one.
  Backs `D`.
- `total_degree(atom)` — all neighbors; equals `heavy_degree + total_h`. Backs
  `X`.
- `total_valence(atom)` — summed bond orders over all neighbors. Backs `v`.
- `heavy_valence(atom)` — summed bond orders to heavy neighbors.
- Ring membership and ring-connection counts, evaluated over the heavy view.

#### 3.4.3 Daylight compatibility as a matching dialect

Rather than diverging unconditionally, make it selectable. The in-memory model
stays pure; a matcher running in Daylight dialect derives an implicit/explicit
assignment at match time using the §3.5.1 predicate and evaluates `h` and `D`
against it.

```rust
match(query, target, SmartsDialect::Daylight)  // h, D per the specification
match(query, target, SmartsDialect::Modern)    // h absent, D is heavy degree
```

Specification-exact semantics are then available on request, the data model stays
uncontaminated, and legacy pattern files — the case that actually wants legacy
semantics — get them.

#### 3.4.4 Transforms need a hydrogen reconciliation policy

Queries only read. Transforms write, and that is where always-explicit imposes a
real obligation.

A rule authored implicit-H — a 1,3 keto-enol shift, say — changes heavy-atom bond
orders and expects hydrogen counts to follow from valence. Against an
always-explicit graph, some hydrogen node must physically move from carbon to
oxygen. The engine has to:

1. Match over the heavy view.
2. Apply the authored connectivity and bond-order changes to heavy atoms.
3. Recompute the hydrogen count each affected atom should now carry.
4. Reconcile: add, remove, or move actual hydrogen nodes to match.

Steps 3–4 reintroduce valence inference, but scoped to the transform engine
rather than pervading the data model, which is the right place for it.

Step 4 carries a genuine ambiguity: **which hydrogen moves?** Irrelevant when
candidates are chemically equivalent, decisive when one is deuterium, and a rule
authored in implicit-H terms cannot express the choice. RDKit has the same
ambiguity and handles it with the `_isotopicHs` side channel. Here it must be
settled as a stated policy, which §3.6 largely supplies: transforms relocate only
collapsible hydrogens, so an isotope-labeled hydrogen is never silently moved.

### 3.5 Serialization

#### 3.5.1 The collapsibility predicate

Round-tripping needs a per-hydrogen decision that per-atom provenance cannot
express: writing `[2H]C` means folding three hydrogens into a count while
materializing the fourth.

That decision is **derived at write time, not stored**. A hydrogen is written as
a node if and only if it carries information that count notation cannot hold:

- a non-default isotope;
- a non-zero formal charge;
- an atom map number;
- degree other than one;
- a stereo or wedge role on its bond;
- membership in a substance group.

As a writer policy this destroys nothing, tracks no partial state, and expresses
exactly what belongs in a serializer — what the target format can represent.

#### 3.5.2 The predicate is a floor, not a rule

Because a writer needs this logic regardless, it should expose the choice rather
than hard-code one answer. The predicate determines the *minimum* set of
hydrogens that must appear as nodes; anything above that is style:

```text
required(H)  ⊆  written(H)  ⊆  all(H)
```

A writer policy selects within that range:

```rust
enum HydrogenOutput {
    Minimal,   // required only
    All,       // every hydrogen node
    Polar,     // required, plus hydrogens on N, O, S
    AsLoaded,  // required, plus whatever the source materialized
}
```

**Policy cannot produce lossy output.** The floor is enforced whatever the caller
asks for, so a request for `Minimal` still writes an isotope-labeled hydrogen as
a node. Policy varies style; it cannot vary correctness.

The policy space is richer than minimal-versus-all, which is itself the argument
for making it explicit. `Polar` is widespread — PDB convention, docking
preparation, force-field setup — and depiction often wants only the hydrogens
that define stereo. Neither is derivable from the floor; both are chemistry-driven
choices belonging to the caller. `AsLoaded` covers byte-fidelity round-tripping,
taking a stored per-file annotation as its input.

#### 3.5.3 The floor is format-dependent

`required` is parameterized by target format, because formats differ in what
their count mechanism can carry. A molfile `HCOUNT` field cannot hold an isotope;
SMILES cannot atom-map an individual hydrogen, since `[CH3:1]` maps the carbon.
The shape of the predicate is constant, its result is not.

Formats also fall into two classes that differ in kind:

- **With a count mechanism** (SMILES, molfile): suppression folds hydrogens into
  a count and loses nothing.
- **Without one** (PDB, XYZ): there is no way to state that a carbon bears three
  hydrogens, so suppression *discards* them. The floor is effectively
  all-or-nothing, and minimal output is lossy by construction.

The second class should be flagged as lossy at the API level rather than silently
accepted, since "write a PDB" quietly meaning "discard all hydrogen information"
is the kind of default this spec exists to eliminate.

Writing hydrogens as nodes also requires them to have coordinates, so
`HydrogenOutput::All` against a molfile whose conformer lacks hydrogen positions
(§3.7.2) is unsatisfiable — the caller must place them first. The type system can
carry that constraint rather than surfacing it at runtime.

### 3.6 Hydrogen isotopes

Deuterium is hydrogen. It is also non-collapsible. These are two independent
facts, the model already tracks both, and they do not interact:

| Fact | Mechanism | Consequence |
| --- | --- | --- |
| D is hydrogen | atomic number 1 | counts in `total_h`, not in `heavy_degree`; the degree identity holds; isotope-blind search keeps its recall |
| D is non-collapsible | carries a non-default isotope | always written as a node; never relocated by a transform; visible to any isotope-aware algorithm |

No new rule is needed for either — "carries a non-default isotope" is already the
first entry in the §3.5.1 predicate.

Treating D and T as heavy atoms instead would be a second mechanism for a
property the model already expresses, and an actively harmful one. It breaks the
degree identity `total_degree = heavy_degree + total_h`: for the methyl carbon of
deuterated acetic acid, D-heavy-and-still-counted-as-H gives `4 + 3 = 7` against
a `total_degree` of 4. So "D is a heavy atom" and "`[CH3]` matches CD3" cannot
both hold, and RDKit takes the second branch — `[CH3]` matches the CD3 carbon,
because `getTotalNumHs` counts neighbors by atomic number without regard to
isotope.

That default is load-bearing. Substructure search is isotope-blind today, so an
acetic acid query matches the deuterated compound. Under a global D-is-heavy rule
it would not, and a database search would silently lose deuterated analogues —
false negatives with no error raised, on exactly the compounds a deuterated-drug
program cares about.

RDKit itself declines to pick, maintaining two notions of not-hydrogen side by
side with the disagreement recorded in comments (`QueryOps.h:89`, `:102`):

```cpp
//! D and T are treated as "non-hydrogen" here
queryAtomNonHydrogenDegree  // nbr->getAtomicNum() != 1 || nbr->getIsotope() > 1

//! D and T are not treated as heavy atoms here
queryAtomHeavyAtomDegree    // nbr->getAtomicNum() > 1
```

The inconsistency runs deeper: substructure matching is isotope-*blind* while
Morgan fingerprints are isotope-*aware*. Both defaults are defensible — search
wants recall, fingerprints want discrimination, since a deuterated drug really is
a different molecule — which is evidence that the correct answer is
subsystem-dependent and that any global rule re-creates the problem for whichever
callers lose.

One presentational option remains open: exposing the collapsible partition as an
iteration view, say `skeleton_atoms()` for "heavy atoms plus isotope-labeled
hydrogens," so fingerprints and descriptors can be isotope-aware without each
re-deriving the predicate. That is a convenience API over an existing rule, and
can be added or dropped without affecting anything else.

### 3.7 Representation

The obvious objection to always-explicit is cost. The NCI 5K sample has a median
of 15 heavy atoms and a mean of 16.4; drug-like hydrogen-to-heavy ratios run
roughly 1:1 to 1.3:1, so this roughly doubles atom count and more than doubles
bond count.

#### 3.7.1 Heavy-first storage

Traversal cost is recoverable by partitioning storage:

- Atoms are stored heavy-first: heavy atoms occupy `[0, n_heavy)`, hydrogens
  occupy `[n_heavy, n_atoms)`.
- Neighbor lists are likewise ordered heavy-first.

`heavy_atoms()` is then a contiguous slice scan with no per-atom branch and no
cache pollution from hydrogen records, and `heavy_neighbors()` is a prefix slice
of the neighbor list. The common case — heavy-only traversal, which is most
algorithms — costs what the count representation costs. The remaining penalty is
memory, not time.

The price is that atom indices no longer correspond to input file order and are
not stable across edits, and RDKit consumers depend on both. The intended
resolution is to make indices opaque handles and expose input order as an
explicit property: better design, but a real migration hazard, since it breaks
downstream code silently rather than loudly.

#### 3.7.2 Coordinates are optional per atom

A conformer maps atoms to *optional* positions:

```rust
conformer.position(atom) -> Option<Point3>
```

Absence states that the position is unknown. Nothing is fabricated at load time,
and every reader is forced by the type to handle the missing case.

This is not a hydrogen question. Experimental structures routinely lack positions
for heavy atoms too — disordered regions, unresolved side chains, low-occupancy
alternates — and RDKit has no way to express it. `Conformer` exposes only
`GetAtomPosition`, with no notion of absence, so a missing position must be
encoded as a missing *atom*: a PDB lysine whose side chain is unresolved loads as
four atoms, a residue that is chemically a lysine but structurally a fragment.
Chemical identity and observational completeness are conflated, and the molecule
is silently wrong rather than explicitly incomplete.

`Option<Point3>` handles hydrogens and heavy atoms uniformly. The hydrogen case
is then the most common instance of a general problem rather than a special rule.

Consequences:

- **Hydrogen placement is an explicit operation**, not something loading does
  silently: `mol.place_hydrogens(HPlacement::Idealized)`.
- **Placement quality remains a real problem.** Idealized bond lengths and angles
  determine most hydrogen positions well, but rotatable hydrogens (hydroxyl,
  thiol, amine) depend on the hydrogen-bonding network, which is why tools such
  as `reduce` exist. It is now an opt-in routine whose quality the caller can
  reason about, rather than a hidden step in file parsing.
- **Provenance survives as secondary annotation.** Once positions can be
  deliberately placed, a consumer may want to know whether present coordinates
  were measured or generated. Absence carries the primary signal.
- **Completeness becomes a typestate.** Force fields, RMSD, and shape comparison
  require full geometry:

  ```rust
  conformer.complete() -> Option<CompleteConformer>
  fn mmff_optimize(c: &CompleteConformer) -> Energy;
  ```

Storage note: `Option<Point3>` over `f64` has no niche and would cost eight bytes
of padding per atom. A dense coordinate array plus a presence bitset is the
obvious representation; the API stays `Option`.

#### 3.7.3 Depictions are not conformers

RDKit stores 2D layouts as a `Conformer` with `is3D` false. A depiction is a
rendering artifact, not a conformation, and conflating them is why "must a
conformer have hydrogen coordinates?" looks ambiguous — the answer differs for
the two things sharing one type.

Separating them resolves it: `Conformer` is 3D molecular geometry; `Depiction` is
a 2D layout for rendering, heavy-atom by construction since hydrogens are mostly
not drawn. Each gets a clean invariant instead of a shared weak one.

#### 3.7.4 Conformer libraries

A common workflow embeds with hydrogens, then removes them to shrink a conformer
library. Always-explicit cannot drop the hydrogen nodes, but optional coordinates
make the workflow expressible as "drop the hydrogen *positions*, keep the
hydrogen nodes."

That recovers most of the benefit: a multi-conformer library stores the graph
once and the coordinates N times, so coordinates dominate, and discarding roughly
half of them is close to what stripping the atoms achieved. What remains is one
molecule's worth of hydrogen nodes per library rather than per conformer. Worth
measuring rather than asserting.

### 3.8 What the type system carries

Hydrogens are not a molecule state, so no `HExplicit`/`HImplicit` parameter is
needed:

```rust
fn embed_conformer(mol: &Mol<Sanitized>) -> Conformer;
fn morgan_fingerprint(mol: &Mol<Sanitized>) -> Fingerprint;
```

Sanitization state remains worth encoding, for the reason `Atom.cpp:303` exists:
"did you run the required preprocessing pass" is a static property currently
checked at runtime, and typestate makes that class of bug a compile error at no
runtime cost. This is the strongest argument for a clean-room implementation in
an expressive type system — stronger than memory safety.

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
`removeAndTrackIsotopes`, and the `_isotopicHs` side channel. Here it needs no
special handling: the deuterium is a node, and its isotope label lives on it
permanently.

## 5. What this buys

- One derived accessor instead of three, with no correctness-relevant default
  argument and no cache.
- `noImplicit` disappears; provenance replaces it, split across `Mol`
  (serialization) and `Query` (matching semantics).
- The fifteen removal flags and the `_isotopicHs` side channel disappear.
- Every property has exactly one home, because the node always exists.
- Hydrogen output becomes a writer policy over a correctness floor, so `Minimal`,
  `All`, and `Polar` are one parameter rather than three code paths, and none can
  emit lossy output by accident.
- Queries and transforms are still authored implicit-H, so the SMIRKS catalogs
  port unchanged (subject to §3.4.4).
- `[CD1]` matches a methyl carbon unconditionally rather than silently ceasing to
  match after `AddHs` — a deliberate divergence from Daylight (§3.4.2), not a bug
  fix.
- Missing coordinates are expressible, so an unresolved side chain no longer
  silently becomes a different molecule.
- Hydrogen counts never change implicitly. In RDKit, editing a molecule
  recomputes them on the next sanitize; here, adding or removing hydrogen is a
  deliberate graph edit.
- "Forgot to sanitize" becomes a compile error.
- Stereo perception should simplify: a tetrahedral center always has four real
  neighbors and a wedge bond always has a real endpoint, so the phantom-neighbor
  special-casing threaded through `Chirality.cpp` (4k lines) should shrink.
  Unmeasured — a hypothesis to check, not a claim.

## 6. Open questions

- **Atom index stability.** Heavy-first partitioning breaks input-order
  correspondence and edit stability. Opaque handles plus an explicit input-order
  property is the intended answer, but the failure mode for ported code is
  silent. Needs a migration story before adoption.
- **How much does `Query` duplicate?** If `Query` and `Mol` share little the fork
  is cheap; if matching, traversal, and serialization each need two
  implementations it is expensive. Worth prototyping before committing.
- **How much existing SMARTS does the `D` rebinding break?** It preserves intent
  for patterns meaning "methyl-like", but any pattern authored against an
  explicit-H molecule means the opposite. Needs a survey of real-world corpora,
  not just RDKit's catalogs.
- **The reconciliation policy (§3.4.4).** §3.6 settles the isotope case; the
  remaining choice among equivalent hydrogens still needs stating, and the stated
  policy needs checking against RDKit's emergent behavior on the tautomer
  catalog.
- **Memory at scale.** Doubling atom count matters for billion-compound
  libraries. Does a compressed on-disk form that inflates on load suffice, or
  does the in-memory representation need a count-based variant for screening?
  Reintroducing one would bring back most of §2.2.
- **Hydrogen placement quality (§3.7.2).** A fixed idealized-geometry rule is
  cheap and wrong for rotatable hydrogens; a network-aware pass in the spirit of
  `reduce` is better and much more expensive. Possibly both, selected by the
  caller.
- **How far does optionality propagate?** If positions are optional, so are
  derived geometric quantities — bond lengths, angles, torsions, centroids. Does
  every geometry accessor return `Option`, or does `CompleteConformer` become the
  ordinary currency with partial conformers confined to loading and
  serialization? The second is likely right, but needs deciding once rather than
  case by case.
- **Conformer library size (§3.7.4).** Does dropping hydrogen positions recover
  what stripping hydrogen atoms used to save? Should be close, since coordinates
  dominate, but measure it.
- **Query hydrogen constraints.** SMARTS allows ranges and negation, so the
  constraint type is richer than an optional integer.
- **Valence model coupling.** Inferred counts depend on the valence model, which
  is itself unspecified in RDKit. This spec assumes one exists; that is the next
  document, and the immediate blocker.

## 7. Validation plan

1. Implement `total_h` as a *Python* shim over RDKit before writing any Rust: add
   hydrogens, count H neighbors.
2. Run it over ChEMBL, comparing every atom to
   `getTotalNumHs(includeNeighbors=True)` on the original molecule.
3. Classify every disagreement as (a) spec is wrong, (b) RDKit is wrong, or (c)
   legitimate design divergence.
4. Round-trip every molfile in the corpus and confirm asserted counts survive
   `AddHs` and writing.
5. Enumerate the cases where RDKit's default `RemoveHs` declines to remove a
   hydrogen, and confirm each corresponds to information the always-explicit
   model simply retains.
6. Measure memory and traversal cost of always-explicit with heavy-first
   partitioning against RDKit on a fingerprinting benchmark, to test the §3.7.1
   claim that traversal cost is recoverable.
7. Survey real-world pattern corpora outside RDKit for `D` and `h` usage — the
   population the rebinding actually risks. (The RDKit-internal half of this is
   done; see §3.4.2.)
8. Replay the tautomer catalog under a stated reconciliation policy and diff
   against RDKit's output.

Steps 1–5 are deliberate: they cost days in Python, require no Rust, and a high
disagreement rate is the signal to stop before committing to an implementation.

## 8. Design history

The model reached its current form through three retracted positions. Each is
recorded because the retractions are the useful part — they show which
distinctions were real.

**Materialization as a molecule state.** The first draft modeled a molecule as
materialized, dematerialized, or mixed, with `dematerialize` as a partial
operation, and derived the fifteen `RemoveHsParameters` flags from a
collapsibility predicate. Deriving fifteen flags from one rule beat enumerating
them, but the draft claimed materialization was a view rather than a property of
the molecule and then contradicted itself by making it a molecule state with a
lossy transform. Taking the original claim seriously gives §3.1: hydrogens are
always nodes, the view belongs to iteration, and the predicate is unnecessary in
the data model. It reappears in §3.5.1 as a *writer* policy, which is where it
was always correct.

**The `Mol`/`Query` split as expensive.** This assumed SMIRKS catalogs would need
rewriting against explicit-H targets. Wrong: the view abstraction extends to
authoring, so patterns stay implicit-H (§3.4). What survives is narrower and
lives in transforms rather than queries.

**`D` and `h` as RDKit artifacts.** The §3.4.1 finding — that SMARTS already
carries both views as separate primitives — was first recorded as evidence the
design was right, the language wanting a model the data model failed to supply.
Checking the Daylight specification corrected it: both are standard primitives,
faithfully implemented. The dual representation is in the query language
standard, so always-explicit requires deliberate divergence (§3.4.2) rather than
cleanup, which moves the hydrogen model partly out of Tier 1 and into Tier 2.

**Conformers as heavy-atom prefix arrays.** Proposed so that a hydrogen-less 3D
file would need no fabricated geometry, with an optional hydrogen extension. This
invented partial state to solve a problem that mostly does not exist, since
RDKit's `Conformer` already covers every atom and embedding requires hydrogens
regardless. The replacement — the loader places idealized hydrogens and records
the fabrication as provenance — was also wrong, for a sharper reason: it
fabricates geometry and then relies on consumers checking a flag to discover it
is not real, the same failure shape as `getTotalNumHs`'s default argument.
Optional coordinates (§3.7.2) invert it. The earlier rejection of optional
positions rested on a category error: this spec bans partial *representation*,
where a half-materialized molecule is ambiguous about what it denotes, whereas a
missing coordinate is missing data.

**D and T as heavy atoms.** Considered and rejected in §3.6, which also explains
why the question dissolves rather than needing an answer.

## 9. Rejected alternatives

### A per-node `is_implicit` flag

Keep every hydrogen in the graph, so properties have a home, but tag each node
implicit or explicit. The attraction is that it restores referents for Daylight
`h` and `D` and removes the divergence of §3.4.2.

The governing test:

> A per-node annotation is acceptable if and only if nothing but the serializer
> reads it. Once matching reads it, the dual representation is back.

The flag fails by construction. For `h` and `D` to work it must be live during
matching, so every molecule carries an implicit/explicit assignment that changes
match results, and `[CD1]` again answers differently depending on how the flags
happen to be set. The flag does not fix the `D` instability; it re-implements it,
buying specification compliance by preserving the defect.

The rest of §2.2 returns with it. `RemoveHs` becomes "flip explicit to implicit,"
reviving the question *may this hydrogen be flipped?* — an isotope-labeled
hydrogen marked implicit cannot be written as `[CH3]`, so the writer must
override the flag. Advisory state silently overridden by a downstream layer is
the same footgun family as `getTotalNumHs`'s default argument. And algorithms
face three views rather than two.

What the flag is right about — that round-tripping needs a per-hydrogen decision
provenance cannot express — is handled by §3.5.1, deriving it at write time
instead of storing it.
