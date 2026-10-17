# 2026-10-17 — Next-generation RDKit: the molecule API

- **Date:** 2026-10-17
- **Author:** Steven Kearnes
- **Acknowledgments:** Prepared with [Claude Code](https://claude.com/claude-code) (Claude Opus 5.5)
- **Status:** draft (working)
- **Tags:** rdkit, next-generation, api-design, cpp20, nanobind, proposal
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

## Summary

Maps RDKit's molecule API, C++ and Python, and proposes a clean-sheet
replacement for `rdkit.v3` beside it. The scope is the molecule and what sits
next to it: atoms, bonds, conformers, stereo and substance groups, ring info,
properties, and the `MolOps` functions called on a molecule. Queries, parsers,
and writers are left for later passes. The proposal:

- A mutable, unchecked `MolBuilder` becomes an immutable, checked `Mol` in one
  `build()` that reports every problem it finds. Checks can be relaxed
  explicitly, but every value on a `Mol` stays defined.
- Atoms and bonds are views into the `Mol`; per-atom data and coordinates are
  contiguous, so Python gets NumPy views. Hydrogens are explicit atoms or an
  implicit count, nothing else.
- Algorithms are free functions that return new values, configured by option
  structs, with no process-wide switches.
- Python is a mechanical binding of the C++ API, which is shaped for nanobind.

The design assumes `ROMol` will be replaced. Wrapping it is the migration path,
and conformer storage has to break from it at the start.

Follows the
[2026-09-12 next-generation entry](../2026-09-12-rdkit-next-generation/README.md),
which settled shape and mechanics and left the API itself for later. Line
references are to upstream RDKit `master` at `7f3d37d` (2026-10-16); bare file
names resolve as described in [CURRENT_API.md](assets/CURRENT_API.md).

## Principles

1. **Thin Python.** Everything is written in C++, and Python is a
   near-mechanical binding of it: one Python name per C++ name, both snake_case,
   with no Python-only behavior beyond notebook display. The C++ API is shaped
   so that nanobind binds it with almost no glue (see [Python](#python)).
2. **Clean sheet, clean break.** The design assumes `ROMol` will be replaced
   and doesn't bend to what wrapping it can deliver. Wrapping existing code
   stays the migration path (Decision 2 of the
   [2026-09-12 entry](../2026-09-12-rdkit-next-generation/README.md)), and
   wherever the design and wrapping conflict, the design wins and the
   implementation behind that piece is replaced
   ([where wrapping falls short](#where-wrapping-romol-falls-short)).
3. **Building is separate from using.** A mutable, unchecked builder becomes a
   checked, immutable molecule in one step that can fail. Unsanitized molecules
   stay reachable through an explicit escape hatch.
4. **`const` means unchanged.** A finished molecule has no state that changes
   after it is built. Operations such as adding hydrogens return a different
   molecule.
5. **Two kinds of hydrogen.** A hydrogen is either an atom in the graph
   (explicit) or part of a count on its heavy atom (implicit). This takes from
   the July
   [hydrogen model spec](../2026-07-18-rdkit-rust-reimplementation/assets/HYDROGEN_MODEL_SPEC.md),
   which is otherwise deferred, its finding that today's four pieces of
   hydrogen state encode fewer concepts than they appear to.
6. **A small type, free functions, no switches.** The molecule holds data;
   algorithms are free functions that take it and return values. Options are
   named structs with written-out defaults, and nothing reads a process-wide
   switch.
7. **Safe to share.** A finished molecule can be read from many threads at
   once, so Python can release the GIL around any call that only reads it.
8. **Arrays where arrays fit.** Coordinates and per-atom columns are
   contiguous, so Python sees NumPy views of them without copying.
9. **Typed errors.** Failures are exceptions from a small hierarchy that maps
   one-to-one onto Python exceptions, and a failed check reports every problem
   with its atoms, not only the first.

The Boost Graph Library may be used inside the implementation; its types stay
out of the public API because nanobind can't bind them without glue.

## Today's API

The full inventory, by concern and with line references, is
[CURRENT_API.md](assets/CURRENT_API.md). This section summarizes it.

### Size

| Class | C++ public members | Python methods (nanobind) |
| --- | --- | --- |
| `ROMol` / `Mol` | 108 (68 names) | 64 (48 names) |
| `RWMol` | +22 | +18 |
| `Atom` | 59 | 73 |
| `Bond` | 45 | 49 |
| `Conformer` | 17 | 32 |
| `RingInfo` | 39 | 26 |
| `MolOps.h` / `rdmolops` | 83 functions (62 names) | 109 (100 names) |

Overloads count separately. `ROMol`, `Atom`, `Bond`, `Conformer`, and
`SubstanceGroup` also inherit `RDProps`' 12 property members.

### Shape

- A molecule is a Boost Graph Library `adjacency_list` of heap-allocated
  `Atom *` and `Bond *` (`ROMol.h:57-58`). Atoms, bonds, conformers, and
  substance groups each point back to their molecule.
- `RWMol` adds editing to `ROMol` and declares no data members of its own.
- Sanitization is a separate 14-step function (`MolOps.cpp:609-719`) that
  callers may skip. Valences, ring info, aromaticity, conjugation,
  hybridization, CIP labels, and matrices live on the molecule as caches,
  fields, or properties that `MolOps` functions fill.
- Python mirrors C++ in CamelCase with about a dozen renames (`GetSSSR` for
  `findSSSR`, `Cleanup` for `cleanUp`), and its module functions mutate a `Mol`
  in place.

### What the proposal responds to

Lines marked *observed* were also run, in a nanobind build of the author's
working branch; the rest were read in source.

| Principle or choice | Today |
| --- | --- |
| 3. Building is separate from using | Valence-derived reads throw until `updatePropertyCache` runs (`Atom.cpp:300-335`), and `needsUpdatePropertyCache` asks only whether the cache was ever filled (`Atom.cpp:744-747`). Ring info isn't reset by `addAtom`, `addBond`, `insertMol`, `replaceAtom`, or `replaceBond`; *observed*: closing hexane into a ring with `AddBond` leaves 0 rings and `IsInRing()` false. `sanitizeMol`'s header lists 11 of its 14 steps (`MolOps.h:541-552`). |
| 4. `const` means unchanged | `setName`, `setProp`, `clearProp`, and `clearComputedProps` are `const`; `getRingInfo() const` and every `getOwningMol() const` hand out mutable access. `findSSSR`, `fastFindRings`, `findRingFamilies`, the matrix functions, `assignAtomCIPRanks`, and `MolToSmiles` write to a `const ROMol &`. Python's `Atom.IsInRing()` runs `findSSSR` (`nb/Atom.cpp:78-90`). |
| 1. Thin Python | Thirteen Python functions, `SanitizeMol` and `Kekulize` among them, mutate a `Mol` through `static_cast<RWMol &>` (`nb/MolOps.cpp:417-498,1097-1105`), while `AddHs` and `RemoveHs` return new molecules. Defaults drift: `Kekulize(clearAromaticFlags=False)` inverts C++ `markAtomsBonds=true`, and `SetAtomMapNum` defaults `strict` the other way. Glue has its own bugs; *observed*: `RWMol(pickle, AllProps)` drops every property, because the flags reach `quickCopy` through two implicit conversions (`nb/Mol.cpp:1067-1081`, `ROMol.h:538`). |
| Views that can't dangle | `reference_internal` keeps the molecule alive, not the element, so atom, bond, conformer, and substance-group wrappers outlive `RemoveAtom`, `RemoveConformer`, or a vector reallocation. Sequence holders fix their range when created (`nb/seqholders.hpp:31-32`), and `Mol.__setstate__` destroys and rebuilds the live object (`Wrap_nb.h:170-171`). |
| Plain values, one reference model | Stereo groups and bookmarks hold pointers; substance groups, stereo atoms, and ring info hold indices; haptic bonds keep endpoints in a string property. `replaceBond` repoints none of the stereo groups that hold the bond it deletes (`RWMol.cpp:277-325`). Copies treat the owner pointer three ways (`Atom.cpp:201`, `Bond.cpp:27,52`, `Conformer.h:58-76`), and `QueryAtom`'s move assignment calls itself (`QueryAtom.h:79-86`). |
| 5. Two kinds of hydrogen | Four pieces of state; `getTotalNumHs(includeNeighbors = false)`; 16 fields in `RemoveHsParameters` (`MolOps.h:302-331`); `removeHs` runs a full sanitization by default. |
| 6. No switches, options as structs | Legacy stereo, non-tetrahedral stereo, and legacy ring finding are environment variables, and any value but `"0"` reads as true (`utils.h:96-106`). Pickling drops properties unless a global is changed (`MolPickler.cpp:78`). Boolean runs such as `assignStereochemistry(mol, bool, bool, bool)`; out-parameters on `sanitizeMol` and `getMolFrags`; and seven functions that mutate an `RWMol &` but copy a `const ROMol &` under one name. |
| 7. Safe to share | No core molecule file takes a lock, and the `const` writes above make concurrent reads of one molecule a race. |
| 8. Arrays where arrays fit | Conformer positions are `Point3D`s with vtable pointers (`Conformer.h:157`), so `GetPositions()` copies; `getBondWithIdx` walks the edge list (`ROMol.cpp:321-334`). |
| 9. Typed errors | Preconditions throw `Invar::Invariant`, which Python sees as `RuntimeError`, after logging (`Invariant.h:108-114`); *observed* for an out-of-range `GetAtomWithIdx`. Sanitization raises on the first problem. *Observed*: `AtomValenceException.GetType()` returns the message. |

## Proposal

Declarations are in [mol.h](assets/mol.h), with calls against them in
[example.cpp](assets/example.cpp); [check.sh](assets/check.sh) compiles both as
C++20 with `-Wall -Wextra -Wpedantic -Werror` and passes with Apple Clang 21 and
GCC 16.2. Names live in `rdkit::v3` and `rdkit.v3` (Decision 6 of the
[2026-09-12 entry](../2026-09-12-rdkit-next-generation/README.md#6-rdkitv3)).

### At a glance

| Today | Proposed | Change |
| --- | --- | --- |
| `RWMol` | `MolBuilder` | Unchecked; holds atom and bond specs, not atoms |
| `ROMol`, sanitized or not | `Mol` | Always checked; immutable; copies share storage |
| `Atom`, `Bond` objects owned by the molecule | `Atom`, `Bond` views | An index and the molecule; no state of their own |
| `atoms()` ranges, Python `_AtomSeqHolder1` | `AtomRange`, `BondRange` | Random-access; Python sequences that keep the `Mol` alive |
| `RingInfo` filled on demand | `RingInfo` built with the `Mol` | Read-only |
| Chiral tags, bond stereo, bond directions, `StereoGroup` | `StereoElement`, `StereoGroup` | Elements name their reference atoms |
| `Conformer` of `Point3D` | `Conformer` of `double` | Contiguous; no owning-molecule pointer |
| `SubstanceGroup` | `SubstanceGroup` | Plain value of index lists |
| `AtomMonomerInfo` subclasses | `ResidueInfo` | Plain optional value on the atom |
| `Dict` of `RDValue` | `Properties` | Closed set of value types; read-only on a `Mol` |
| `MolOps` in-place functions | Free functions returning new values | See [Operations](#operations) |
| `MolSanitizeException` family, `Invar::Invariant` | `Error`, `BuildError`, standard exceptions | See [Errors](#errors) |

### Building, and the escape hatch

- `MolBuilder` is mutable and checks nothing. It holds an `AtomSpec` per atom
  and a `BondSpec` per bond, can start from a finished `Mol`, and can be read
  at any point. Indices stay stable while it is edited: a removed atom leaves a
  gap, and `build()` keeps surviving atoms in order and closes the gaps
  (`built_index()` reports where an atom landed).
- `build(BuildOptions)` covers what today's sanitization does: cleanup
  normalizations, valence, implicit hydrogens and radicals, Kekulé bond orders,
  aromaticity under a chosen model, conjugation, hybridization, rings, and
  stereo validation. It returns a `Mol` or throws a `BuildError` listing every
  problem it found.
- Parsers return a `Mol` and take `BuildOptions`. For input that won't build,
  each parser has a form that returns the `MolBuilder` instead, so the input
  can be read and repaired before building. Parser names are out of scope here.
- **The escape hatch** is `BuildOptions::relax`, a set of checks to skip:
  `VALENCE`, `KEKULIZATION`, `AROMATICITY` (keep the given flags), `CLEANUP`,
  and `STEREO` (keep stereo on atoms that aren't stereocenters). A relaxed
  check never leaves a value undefined: an atom above its allowed valence gets
  zero implicit hydrogens, and a bond that can't be kekulized keeps
  `BondOrder::AROMATIC`. The `Mol` records its relaxations, and an operation
  that can't accept one (anything needing Kekulé orders, say) throws `Error`
  naming it. That is the only runtime check of a molecule's state, and only an
  explicit opt-in reaches it.
- Rings, conjugation, and hybridization can't be relaxed; they're defined for
  any graph.

> [!NOTE]
> **Alternative: a type per state.** `Mol` for checked molecules and a second
> frozen type for unchecked ones, or `Mol<Checked>` as a template.
>
> - **Enables.** "Needs Kekulé orders" checked at compile time.
> - **Not now.** Every operation would need an overload or a template per
>   state, and Python a class per state, against principle 1. Relaxed
>   molecules are rare and opt-in.
> - **Adopt when.** Relaxed molecules become common enough that runtime errors
>   on them are a real cost.

### Atoms, bonds, and columns

- `Atom` and `Bond` are views: a pointer to the `Mol` and an index. Atoms
  exist only inside a molecule, so there is no detached atom and no
  `getOwningMol()`. Because a `Mol` never changes, a view can't go stale.
- In Python every binding that returns a view or range keeps its `Mol` alive
  (nanobind's `keep_alive`, applied by one rule). Today the same annotation
  keeps the molecule alive while an edit deletes the atom under the wrapper;
  here the `Mol` can't change, so keeping it alive is enough.
- Indices are `AtomIndex` and `BondIndex`, distinct C++ types that can't be
  mixed up; Python sees plain `int`.
- Per-atom and per-bond columns (`atomic_numbers()`, `formal_charges()`,
  `implicit_h_counts()`, `bond_atoms()`, `bond_orders()`) return `std::span`,
  which Python sees as read-only NumPy views. `bond_atoms()` is the
  `(num_bonds, 2)` edge list graph models take as input.
- `Atom::neighbor_indices()` is a span into compressed adjacency storage, so
  walking the graph allocates nothing.

### Hydrogens

| Today | Proposed |
| --- | --- |
| Hydrogen atoms in the graph | Explicit hydrogens: `explicit_h_count()` counts them |
| `numExplicitHs` + `noImplicit` (bracket and molfile counts) | `AtomSpec::implicit_h_count` set: the count is given |
| Implicit valence cache, filled by sanitization | `AtomSpec::implicit_h_count` empty: inferred at build |
| `getTotalNumHs(includeNeighbors = false)` | `total_h_count()`: implicit plus explicit, no argument |

- On a `Mol`, `implicit_h_count()` is always defined. Whether it was given or
  inferred survives only as `implicit_h_count_was_given()`, which writers read
  to decide on brackets.
- `add_hs` and `remove_hs` return new molecules. By default `remove_hs` keeps
  any hydrogen that carries something a count can't: an isotope, a map number,
  properties, stereo, a substance group, a charge, or a bond to anything other
  than one heavy atom. An isotopic hydrogen therefore stays explicit, and no
  side channel carries its isotope.
- The word *explicit* changes meaning. Today `GetNumExplicitHs()` returns the
  count stored on the atom, which this proposal calls a given implicit count;
  here *explicit* means a graph atom, as in Daylight's usage. Migration docs
  need to say so.

### Bond orders and aromaticity

- A built bond has a Kekulé `order()` and a separate `is_aromatic()` flag, both
  always set. Kekulization happens once, in `build()`, and picks one valid
  structure deterministically; writers choose aromatic or Kekulé output. There's
  no in-place `Kekulize` or `SetAromaticity`.
- `with_aromaticity(mol, model)` returns a molecule perceived under another
  model.
- `BondOrder` has seven values: `ZERO`, `SINGLE`, `DOUBLE`, `TRIPLE`,
  `QUADRUPLE`, `DATIVE`, and `AROMATIC`, which appears only under
  `Relaxation::KEKULIZATION`. Which of today's other bond types survive is
  open.

### Rings and other derived data

- `build()` computes everything derived that the `Mol` exposes: implicit
  hydrogens, Kekulé orders, aromaticity, conjugation, hybridization, ring
  membership, and stereo validation. A `Mol` holds no lazy cache, so it needs
  no lock to be shared.
- `rings()` returns the symmetrized SSSR. Other ring sets, and anything costly
  that only some callers need (CIP labels, ring families), are free functions
  returning values.

### Stereo

- `stereo()` lists `StereoElement`s, each naming its own reference atoms:
  `TetrahedralStereo` gives its center, four neighbor slots (an empty slot is
  an implicit hydrogen or lone pair), and `CLOCKWISE` or `ANTICLOCKWISE`;
  `DoubleBondStereo` gives its bond, one reference neighbor at each end, and
  `CIS` or `TRANS`. The order atoms or bonds are stored in never changes what
  an element means. Atropisomers and non-tetrahedral centers join the variant
  as new types.
- `StereoGroup`s (`ABSOLUTE`, `AND`, `OR`) refer to elements by index.
- Bond directions (`/` and `\`) are SMILES syntax: parsers turn them into
  elements and writers regenerate them. Wedges are a depiction concern (see
  [Open questions](#open-questions)).
- `build()` runs one perception algorithm, the newer of today's two, with no
  switch to select the legacy one. `possible_stereo()` lists unspecified
  candidates, `cip_labels()` returns a label per element, and
  `stereo_from_conformer()` assigns stereo from coordinates.

### Properties

- `Properties` maps names to a closed set of value types (`bool`, `int64`,
  `double`, string, and vectors of each), each with one Python type, so Python
  reads them as a `dict` with no typed getters.
- Molecules, atoms, bonds, and substance groups have properties. They are set
  on the builder; `Mol::with_properties()` returns a molecule with new
  molecule-level properties and shares everything else.
- Properties hold only what callers put there. Derived data has typed homes,
  so there are no computed properties and no private `_`-prefixed properties.
  Chemistry stored in properties today gets typed fields too: atom map numbers
  (`molAtomMapNumber`), the non-tetrahedral permutation (a stereo element), and
  CIP labels (`cip_labels()`).
- `to_bytes()` writes every property; there's no global setting for which
  ones pickle.

### Conformers

- A `Conformer` is an id, a dimension (2 or 3), and `num_atoms × 3` doubles in
  one block; a 2D conformer has zero `z`. It has no pointer back to a
  molecule.
- Python's `positions` is a read-only `(num_atoms, 3)` NumPy view that keeps
  the `Mol` alive. Changing coordinates means building new conformers
  and calling `Mol::with_conformers()`.

### Substance groups and residues

- `SubstanceGroup` is a plain value: a type, atom and bond index lists, parent
  indices, and properties.
- `ResidueInfo` is an optional plain value on each atom, in place of a
  polymorphic object the atom owns.
- Bookmarks, which reaction code uses to tag atoms, leave the public API.

### Operations

Each operation takes a `const Mol &` and returns a value. The last row is the
builder step that takes over today's in-place perception functions.

| Proposed | Replaces |
| --- | --- |
| `add_hs(mol, AddHsOptions)` | `addHs` |
| `remove_hs(mol, RemoveHsOptions)` | `removeHs` and its parameter struct |
| `with_aromaticity(mol, model)` | `setAromaticity`, `Kekulize` |
| `renumber_atoms(mol, order)` | `renumberAtoms` |
| `combine(a, b)` | `combineMols` |
| `fragments(mol)`, `fragment_atoms(mol)` | `getMolFrags` overloads |
| `stereo_from_conformer(mol, id)` | `assignStereochemistryFrom3D` and relatives |
| `without_stereo(mol)` | `removeStereochemistry` |
| `possible_stereo(mol)` | `findPotentialStereo` |
| `cip_labels(mol)` | `CIPLabeler::assignCIPLabels`, `_CIPCode` properties |
| `total_formal_charge(mol)` | `getFormalCharge` |
| `adjacency_matrix(mol)`, `distance_matrix(mol)` | `getAdjacencyMatrix`, `getDistanceMat` |
| `MolBuilder::build(options)` | `sanitizeMol`, `cleanUp`, `assignRadicals`, `setConjugation`, `setHybridization`, `findSSSR`, `assignStereochemistry`, `detectChemistryProblems` |

### Errors

| C++ | Python | When |
| --- | --- | --- |
| `std::out_of_range` | `IndexError` | An index past the end |
| `std::invalid_argument` | `ValueError` | An argument that's never valid |
| `rdkit::v3::Error` | `rdkit.v3.Error` | Anything else, including an operation refusing a relaxation |
| `rdkit::v3::BuildError` | `rdkit.v3.BuildError` | `build()` failed; `problems` lists each with its atoms and bonds |

nanobind translates the two standard exceptions itself, and the two RDKit ones
need one registration each. Parsers raise rather than return `None`.

### Python

The binding is a set of mechanical rules plus a few generic casters:

| C++ | Python |
| --- | --- |
| `snake_case` function or method | Same name |
| Accessor with no arguments | Read-only property |
| `std::optional<T>` | `T` or `None` |
| `std::span` of numbers, or of `std::array`s of numbers | Read-only NumPy view, 1D or 2D, that keeps the `Mol` alive |
| `std::span` of objects, `AtomRange`, `BondRange` | Sequence that keeps the `Mol` alive |
| `std::variant` | The held alternative |
| `AtomIndex`, `BondIndex` | `int` |
| Spec and options structs | Classes with keyword constructors |
| `enum class` | nanobind enum, same member names; `Relaxation` is a flag |
| `Properties` | `dict` |
| `to_bytes()`, `from_bytes()` | Pickling |

- nanobind 2.15 has no `std::span` caster, so the binding layer provides
  generic casters for spans, the two ranges, and the index types. Everything
  else is one `def` per function.
- Every call that only reads a `Mol` releases the GIL; small accessors may skip
  it if measurement shows the release costs more than the call. Today four
  nanobind wrapper files release it (distance geometry, fingerprint
  generators, force fields, alignment) and the molecule core none.
- Free-threaded Python also needs `FREE_THREADED` passed to
  `nanobind_add_module`, which RDKit's modules don't pass
  (`RDKitUtils.cmake:230`).
- Enum members are UPPER_CASE in both languages, as today's Python enums are
  (`BondType.SINGLE`). CIP labels keep their case (`R`, `r`).
- Spec and options structs give every member a default initializer, so
  designated initializers that skip members compile under GCC's `-Wextra`.

The same task, today and proposed:

```python
# Today
from rdkit import Chem

mol = Chem.MolFromSmiles("OC(=O)c1ccccc1O")
mol_hs = Chem.AddHs(mol)
for atom in mol_hs.GetAtoms():
    print(atom.GetSymbol(), atom.GetTotalNumHs(includeNeighbors=True))
num_rings = mol.GetRingInfo().NumRings()
edges = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]

pentavalent = Chem.MolFromSmiles("C(C)(C)(C)(C)C", sanitize=False)
pentavalent.UpdatePropertyCache(strict=False)
Chem.SanitizeMol(
    pentavalent, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES
)
```

```python
# Proposed (parser name is a placeholder)
from rdkit import v3

mol = v3.parse_smiles("OC(=O)c1ccccc1O")
mol_hs = v3.add_hs(mol)
for atom in mol_hs.atoms:
    print(atom.symbol, atom.total_h_count)
num_rings = mol.rings.num_rings
edges = mol.bond_atoms  # (num_bonds, 2) NumPy view

pentavalent = v3.parse_smiles(
    "C(C)(C)(C)(C)C",
    v3.BuildOptions(relax=v3.Relaxation.VALENCE),
)
```

### Converting to and from existing molecules

- `rdkit/v3/convert.h` is the one header that names existing types:
  `from_rdkit(const RDKit::ROMol &, BuildOptions)` copies into a builder and
  builds, and `to_rdkit(const Mol &)` returns a new `RWMol`. Python has
  `rdkit.v3.from_rdkit()` and `rdkit.v3.to_rdkit()`.
- Converting in maps chiral tags and bond stereo onto elements,
  `numExplicitHs` with `noImplicit` onto given implicit counts, and aromatic
  bond types onto Kekulé orders with flags; it drops computed and private
  properties. Converting out reverses each mapping.

### Where wrapping ROMol falls short

A `Mol` can start as a wrapper holding a private, built `ROMol`. Where the
design asks for something `ROMol` can't give, the wrapper either keeps extra
state beside it or that piece gets its own implementation:

| Proposed | `ROMol` today | Over a wrapped `ROMol` |
| --- | --- | --- |
| Many threads read one `Mol` | `const` routines write: `getRingInfo() const` returns a mutable `RingInfo *` (`ROMol.h:773`), properties are a `mutable Dict` (`RDProps.h:19`), `MolToSmiles` stores output order on its `const` input (`SmilesWrite.cpp:494,747`), and no core file takes a lock | Every wrapped routine that writes through `const` runs on a private copy, or its output becomes a return value; break when copies cost too much |
| Column spans | Atoms and bonds are heap objects behind `Atom *` and `Bond *` in a BGL `adjacency_list` (`ROMol.h:57-58`) | Columns kept beside the `ROMol`, filled at build; safe because neither changes |
| Contiguous conformers | `std::vector<Point3D>` (`Conformer.h:157`), each point carrying a vtable pointer (`point.h:28-53`) | Can't be a view; conformers get their own storage from the start, copied at the boundary for wrapped routines such as embedding |
| Kekulé order and aromatic flag | One `BondType`, where `AROMATIC` replaces the order | Kekulé orders kept in a side column |
| Stereo elements | Chiral tags relative to bond storage order, bond stereo, bond directions | Derived at build and kept beside the tags that wrapped routines read |
| Two kinds of hydrogen | `numExplicitHs`, `noImplicit`, the implicit valence cache, and graph hydrogens | Normalized at build; the original fields stay for wrapped routines |
| No process-wide switches | Stereo perception reads an environment variable that `setUseLegacyStereoPerception` writes with `setenv` (`Chirality.cpp:848-858,2656,2917`) | Needs the switch turned into a parameter (the 09-12 entry's strict alternative to Decision 2) before stereo perception can be wrapped |

Conformer storage breaks first, at the start. Thread safety is the next
pressure: it holds only while no wrapped routine writes to the shared
`ROMol`.

## Open questions

- **Which of today's bond types survive.** `ONEANDAHALF` and its relatives,
  `IONIC`, `HYDROGEN`, `THREECENTER`, the directional dative types, `OTHER`,
  and `UNSPECIFIED` each have some input format or query that uses them.
- **Where wedges and other depiction data live:** on a 2D conformer, on a
  separate depiction value, or as bond properties.
- **Queries.** SMARTS molecules and query atoms and bonds are a separate type
  in this design, left for their own pass; substructure functions would take a
  `Mol` and a query.
- **Which read-only operations accept a `MolBuilder`.** Writers that echo
  input as given are the likely ones.
- **Property values.** Whether the closed set covers what callers store today,
  and whether per-atom properties need columnar storage to stay cheap.
- **Zero-argument accessors as Python properties.** It reads naturally but
  makes the two languages' spellings differ by `()`.
- **Ring sets.** Symmetrized SSSR as the one built-in set, with ring families
  and others as free functions.
- **Sanitization steps beyond the core.** How organometallic cleanup and
  atropisomer cleanup fit `BuildOptions`.
- **Typed homes for the rest of today's property-borne chemistry:** dummy-atom
  labels (`dummyLabel`), R-group labels and aliases, and haptic bond endpoints
  (`_MolFileBondEndPts`).

## Next steps

- Prototype `Mol`, `MolBuilder`, and their nanobind bindings over a wrapped
  `ROMol` for a small subset, and count the binding code per bound function to
  test principle 1.
- Map the parsers and writers the same way, since they produce and consume
  `Mol` and `MolBuilder`.
- Settle the open questions on bond types and depiction data, which the
  parsers depend on.
- Pick the next family to map: queries and substructure search, or
  conformers and embedding.
