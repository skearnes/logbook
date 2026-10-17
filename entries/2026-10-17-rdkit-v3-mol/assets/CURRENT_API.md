# RDKit's molecule API today

- **Date:** 2026-10-17
- **Author:** Steven Kearnes
- **Acknowledgments:** Prepared with [Claude Code](https://claude.com/claude-code) (Claude Opus 5.5)
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

An inventory of RDKit's molecule API, C++ and Python, by concern, followed by
the hazards a clean-sheet design would address. It supports the
[2026-10-17 entry](../README.md).

Conventions:

- Source is upstream RDKit `master` at `7f3d37d16` (2026-10-16). Bare file
  names resolve within `Code/`: `utils.h`, `types.h`, `Dict.h`, `RDValue.h`,
  `RDProps.h`, `Exceptions.h`, and `Invariant.h` are in `Code/RDGeneral/`,
  `Wrap_nb.h` is in `Code/RDBoost/`, `details.h` is `Code/GraphMol/details.h`,
  and other bare names are in `Code/GraphMol/` or its subdirectories.
- `nb/X.cpp` is `Code/GraphMol/nbWrap/X.cpp` (nanobind); `bp/X.cpp` is
  `Code/GraphMol/Wrap/X.cpp` (Boost.Python). "Python" means the nanobind
  bindings.
- Everything was read in the source. **(observed)** marks a behavior also
  run in a nanobind build of the author's working branch (15 commits past
  `23378a7`), **(inferred)** a runtime consequence that follows from the cited
  code but wasn't run, and **(unverified)** what couldn't be confirmed.
- Method counts are one per declaration, so overloads count separately.
  Constructors, destructors, and assignment operators are excluded unless a
  row says otherwise.

---

## 1. Construction and copying

### C++

| Member | Location | Notes |
|---|---|---|
| `ROMol()` | ROMol.h:518 | `initMol()` allocates an empty `RingInfo` (ROMol.cpp:168-171) |
| `ROMol(const ROMol&, bool quickCopy=false, int confId=-1)` | ROMol.h:531-536 | also the copy constructor; body ROMol.cpp:75-166 |
| `ROMol(const std::string &pickle)` | ROMol.h:538; ROMol.cpp:60-65 | **not `explicit`**; unpickles with AllProps (MolPickler.h:200-203) |
| `ROMol(const std::string&, unsigned propertyFlags)` | ROMol.h:540; ROMol.cpp:67-73 | not `explicit` |
| `ROMol(ROMol&&) noexcept` / `operator=(ROMol&&)` | ROMol.h:542-605 | re-points atom/bond/conformer/sgroup owners; moved-from `dp_ringInfo` set to `nullptr` (565, 576) |
| `operator=(const ROMol&) = delete` | ROMol.h:607-608 | copy-construct allowed, copy-assign not |
| `RWMol()`, `RWMol(const ROMol&, bool, int)`, `RWMol(const RWMol&)`, `operator=(const RWMol&)`, move ops | RWMol.h:35-54; RWMol.cpp:135-142 | `RWMol(const ROMol&, …)` is a non-explicit converting constructor |
| `RWMol::clear()` | RWMol.h:204-209 | `destroy()` + `initMol()`; does not reset batch-edit bitsets |
| `Atom()`, `Atom(unsigned)`, `Atom(const std::string&)` | Atom.h:124-129 | |
| `Atom(const Atom&)`, `operator=` | Atom.h:130-131; Atom.cpp:198-233 | `initFromOther`: drops owner/index, copies props, monomer info, **cached valences**, flags |
| `Atom(Atom&&) = default`, `operator=(Atom&&) = default` | Atom.h:134-135 | see H-C4 |
| `virtual Atom *Atom::copy() const` | Atom.h:143; Atom.cpp:254-257 | caller owns |
| `Bond()`, `Bond(BondType)`, copy, move | Bond.h:109-153; Bond.cpp:20-71 | copy ctor drops owner (Bond.cpp:27); copy-assign keeps owner (52) |
| `virtual Bond *Bond::copy() const` | Bond.h:160 | caller owns |
| `Conformer()`, `Conformer(unsigned)`, copy `= default`, move | Conformer.h:51-76 | copy and move both keep `dp_mol` |
| `QueryAtom(const Atom&)`, `QueryBond(const Bond&)` | QueryAtom.h:34-54; QueryBond.h:36-37 | explicit; build a query from the plain atom/bond |

What the copy constructor copies (ROMol.cpp:75-166):

| Item | full copy | `quickCopy=true` | Line |
|---|---|---|---|
| atoms, bonds (each with its own props) | yes | yes | 85-94 |
| RingInfo | copied | copied | 97-102 |
| StereoGroups | yes | yes | 105-119 |
| pending batch-edit bitsets `dp_delAtoms/dp_delBonds` | yes | yes | 121-130 |
| conformers | yes (only `confId` if `>= 0`) | no | 134-138 |
| SubstanceGroups | yes | no | 141-143 |
| mol props | yes | reset to empty `__computedProps` | 145 / 159-161 |
| bookmarks | yes | no | 148-157 |

Pickling: `MolPickler::pickleMol` / `molFromPickle` (MolPickler.h:168-216). The default pickle property set is a process-wide static, `defaultProperties = PicklerOps::NoProps` (MolPickler.cpp:78), with get/set at 270-283. Its mutex exists only under `RDK_BUILD_THREADSAFE_SSS`. Optional boost::serialization `save/load` (ROMol.h:1005-1014; ROMol.cpp:696-720).

Parser entry points (names only):

- `RDKit::v2::SmilesParse::MolFromSmiles` / `MolFromSmarts` return `std::unique_ptr<RWMol>` (SmilesParse.h:69, 72).
- The v1 inline `SmilesToMol` / `SmartsToMol` return a raw `RWMol*` (SmilesParse.h:118, 171, 189, 217).
- Python: `MolFromSmiles` (nb/rdmolfiles.cpp:1618, 1651), `MolFromSmarts` (1683, 1708).

### Python

| Name | Binding | Semantics |
|---|---|---|
| `Mol()`, `Mol(bytes)`, `Mol(bytes, propertyFlags)`, `Mol(str)`, `Mol(str, propertyFlags)`, `Mol(mol, quickCopy=False, confId=-1)` | nb/Mol.cpp:430-458 | 6 ctors |
| `Mol.__copy__`, `__deepcopy__` | nb/Mol.cpp:459-460; Code/RDBoost/Wrap_nb.h:304-326 | full C++ copy, plus the Python `__dict__` (class has `nb::dynamic_attr`, nb/Mol.cpp:426-428) |
| `Mol.__getstate__` / `__setstate__` | nb/Mol.cpp:1027-1028; Wrap_nb.h:149-180 | pickles with the global default props plus `__dict__`. `__setstate__` runs `self->~T(); new (self) T(pkl)` on the live object (Wrap_nb.h:170-171) |
| `Mol.ToBinary()`, `ToBinary(propertyFlags)` | nb/Mol.cpp:1014-1019 | |
| `PropertyPickleOptions`, `GetDefaultPickleProperties`, `SetDefaultPickleProperties` | nb/Mol.cpp:223-242 | process-wide |
| `RWMol(...)`, 6 ctors | nb/Mol.cpp:1059-1086 | `RWMol(pkl, propertyFlags)` quirk, see H-B1 |
| `RWMol.__copy__/__deepcopy__/__getstate__/__setstate__` | nb/Mol.cpp:1088-1089, 1130-1132 | |
| `RWMol.GetMol()` | nb/Mol.cpp:1116-1117, 167-170 | new `Mol` copy, `take_ownership` |
| `EditableMol(m)`; `GetMol()` | nb/EditableMol.cpp:104, 134; class 23-78 | holds a private `RWMol*` copy; `GetMol` copies again |
| `Atom(what: str)`, `Atom(other)`, `Atom(num: int)` | nb/Atom.cpp:176-179 | no zero-arg ctor |
| `Atom.__copy__` | nb/Atom.cpp:184-185 | `Atom::copy`, `take_ownership`; no `__deepcopy__` |
| `Bond` | nb/Bond.cpp:104-107 | not constructible from Python; no copy |
| `Conformer()`, `Conformer(numAtoms)`, `Conformer(other)` | nb/Conformer.cpp:107-110 | the copy keeps the owner pointer (H-C7) |

Boost difference: `bp/Mol.cpp` pickles through `def_pickle(mol_pickle_suite)` with `getinitargs` (bp/Mol.cpp:68-72, 932, 1042), not `__getstate__`, using the same default-props semantics.

---

## 2. Graph access and iteration

### C++

- **Storage**
  - `MolGraph = boost::adjacency_list<vecS, vecS, undirectedS, Atom*, Bond*>` (ROMol.h:57-58), held in private `d_graph` (1017).
  - `numBonds` is a protected counter (1034).
- **Index access**
  - `getAtomWithIdx` ×4 (ROMol.h:624-636): O(1), with `URANGE_CHECK` (ROMol.cpp:201-218).
  - `getBondWithIdx` ×4 (647-659): walks the edge list, O(idx) (ROMol.cpp:321-334).
  - `getBondBetweenAtoms` ×4 (661-675): `nullptr` if absent (ROMol.cpp:342-354).
  - The template overloads go through `rdcast`.
- **Counts**
  - `getNumAtoms()` (617).
  - `getNumAtoms(bool onlyExplicit)` (620): when `false`, adds `getTotalNumHs()` per atom (ROMol.cpp:180-190).
  - `getNumHeavyAtoms` (622).
  - `getNumBonds(bool onlyHeavy = 1)` (645; ROMol.cpp:303-314).
  - `getAtomDegree(const Atom*)` (638).
- **Neighbors**
  - `getAtomNeighbors` returns an `ADJ_ITER_PAIR` (793); `getAtomBonds` returns an `OBOND_ITER_PAIR` (824).
  - Both PRECONDITION that the atom is owned by this mol (ROMol.cpp:356-368).
- **BGL exposure**
  - `getTopology()` returns `const MolGraph&` (876).
  - `getVertices()` / `getEdges()`, const and non-const (841-861).
  - `operator[](vertex_descriptor)` and `operator[](edge_descriptor)` (979-985).
  - BGL typedefs (354-364).
- **C++20 ranges** (421-516)
  - `atoms()`, `bonds()`, `atomNeighbors(at)`, `atomBonds(at)`, each const and non-const.
  - `checkedAtoms()`, `checkedBonds()`, `checkedAtomNeighbors(at)`, `checkedAtomBonds(at)`.
  - `CXXAtomIterator` is random-access (116-245; `static_assert` at 247-250). `CXXBondIterator` is bidirectional and its `size()` is O(n) (253-343).
  - The checked variants compare vertex/edge counts against counts recorded when `begin()` builds the element iterator (136-158, 271-285), and throw `std::runtime_error("molecule modified during iteration")`.
- **Legacy iterators**
  - 14 typedefs (390-406).
  - begin/end pairs (883-947) for Atoms, Bonds, AromaticAtoms, Heteros, QueryAtoms(query), MatchingAtoms(fnptr), Conformers.
  - `hasQuery()` (918; ROMol.cpp:538-550).
- **Stability**
  - `Atom*` and `Bond*` stay valid until that element is deleted.
  - Indices (`Atom::d_index`, `Bond::d_index`, vertex descriptors) shift on every removal (RWMol.cpp:391-394, 466, 619-623).
  - `getBondWithIdx` relies on edge-list order matching `d_index`; it never reads `d_index`.

### Python

| Name | Binding | Returns / lifetime |
|---|---|---|
| `Mol.GetAtoms()` | nb/Mol.cpp:461-465 | `_AtomSeqHolder1` (342-354) over `checkedAtoms()`; `keep_alive<0,1>` keeps the Mol |
| `Mol.GetBonds()` | nb/Mol.cpp:466-470 | `_BondSeqHolder1` (368-380); `__len__` and `__getitem__` are O(n) each |
| `Atom.GetNeighbors()` | nb/Atom.cpp:272-280 | `_AtomSeqHolder2` (nb/Mol.cpp:355-367) over `checkedAtomNeighbors`; keeps the Atom |
| `Atom.GetBonds()` | nb/Atom.cpp:282-289 | `_BondSeqHolder2` (nb/Mol.cpp:381-393) |
| `Mol.GetConformers()` | nb/Mol.cpp:529-533 | `_ROConformerSeq` (394-406): snapshot `vector<Conformer*>` (seqholders.hpp:182-206) |
| `Mol.GetAromaticAtoms()`, `GetAtomsMatchingQuery(qa)` | nb/Mol.cpp:960-983 | `_ROQAtomSeq` (407-419; seqholders.hpp:208-252) |
| `Mol.GetAtomWithIdx(idx)`, `GetBondWithIdx(idx)`, `GetBondBetweenAtoms(idx1, idx2)` | nb/Mol.cpp:482-554 | raw pointer with `reference_internal` (element keeps the Mol alive); `None` if no bond |
| `Mol.GetNumAtoms()`, `GetNumAtoms(onlyExplicit=True)`, `GetNumHeavyAtoms()`, `GetNumBonds(onlyHeavy=True)` | nb/Mol.cpp:471-498 | |
| `Atom.GetIdx`, `GetDegree`, `GetTotalDegree` | nb/Atom.cpp:195-208 | |
| `Bond.GetIdx`, `GetBeginAtomIdx`, `GetEndAtomIdx`, `GetOtherAtomIdx`, `GetBeginAtom`, `GetEndAtom`, `GetOtherAtom` | nb/Bond.cpp:154-172 | atoms returned with `reference_internal` |

Holder details (nb/seqholders.hpp):

- `__iter__` goes through `nb::make_iterator` with `keep_alive<0,1>`.
- `__getitem__` takes an `int` only, so there is no slicing.
- The atom and bond holders accept negative indices (36-46, 58-68). `QueryAtomIterSeq` and `ConformerIterSeq` reject them (239-251, 196-204).
- The `[vstart, vend)` range is captured when the holder is built (31-32, 53-54; ROMol.h:237-243).
- `ReadOnlySeq`, `AtomCountFunctor`, `ConformerCountFunctor` (seqholders.hpp:71-180) are unused in nbWrap and contain `std::cerr` debug prints.
- Python has no access to BGL, bookmarks, or the legacy iterators.

Boost differences:

- `Mol.GetAtoms` / `GetBonds` are Python monkeypatches that call `GetAtomWithIdx` / `GetBondWithIdx` per index (rdkit/Chem/__init__.py:47-111).
- `Atom.GetNeighbors` / `GetBonds` return tuples (bp/Atom.cpp:69-83, 290-293).

---

## 3. Editing

### C++ (RWMol unless noted)

| Function | Location | Returns | Behavior |
|---|---|---|---|
| `addAtom(bool updateLabel=true)` | RWMol.h:70; RWMol.cpp:217-233 | index | new empty Atom; sets bookmark `ci_RIGHTMOST_ATOM` if `updateLabel`; adds (0,0,0) to each conformer |
| `addAtom(Atom*, bool updateLabel=true, bool takeOwnership=false)` | RWMol.h:83-86 → ROMol.cpp:377-401 | index | copies the atom unless `takeOwnership` |
| `addBond(unsigned, unsigned, BondType=UNSPECIFIED)` | RWMol.h:129-130; RWMol.cpp:501-546 | **new bond count**, not index | `AROMATIC` marks the bond and both atoms aromatic (512-522); clears both endpoints' valence cache (531-532); in batch mode, schedules the new bond if an endpoint is scheduled (536-543) |
| `addBond(Atom*, Atom*, BondType)` | RWMol.h:132-133; RWMol.cpp:548-551 | count | |
| `addBond(Bond*, bool takeOwnership=false)` | RWMol.h:144-146 → ROMol.cpp:403-432 | count | no valence-cache clear, no aromatic propagation |
| `createPartialBond(idx, type)` / `finishPartialBond(idx, bookmark, type)` | RWMol.h:148-185; RWMol.cpp:632-653 | `Bond*` / count | the partial `Bond` is never inserted and never deleted |
| `removeAtom(unsigned)`, `removeAtom(Atom*)`, `removeAtom(Atom*, bool clearProps)` | RWMol.h:114-118; RWMol.cpp:344-499 | | deletes the `Atom` (496) |
| `removeBond(unsigned, unsigned)` | RWMol.h:188; RWMol.cpp:553-630 | | silent no-op if there is no bond; deletes the `Bond` (628) |
| `replaceAtom(idx, Atom*, bool updateLabel=false, bool preserveProps=false)` | RWMol.h:99-100; RWMol.cpp:235-275 | | copies the new atom, deletes the old one (249); `updateLabel` is an unnamed, unused parameter (235) |
| `replaceBond(idx, Bond*, bool preserveProps=false, bool keepSGroups=true)` | RWMol.h:198-199; RWMol.cpp:277-325 | | copies, deletes the old bond (310); lowers endpoint explicit-H counts when bond order rises (292-302) |
| `insertMol(const ROMol&)` | RWMol.h:57; RWMol.cpp:144-215 | | appends copies; offsets bond stereo atoms, `_ringStereoAtoms`/`_ringStereoOtherAtom`, conformer coords (only when this mol has no conformers or the same count), stereo groups (ABS groups merged), sgroups |
| `beginBatchEdit` / `rollbackBatchEdit` / `commitBatchEdit` | RWMol.h:211-216; RWMol.cpp:655-684 | | `begin` throws `ValueErrorException` on re-entry (656-658) |
| `getLastAtom`, `getActiveAtom`, `setActiveAtom` ×2 | RWMol.h:102-112; RWMol.cpp:327-342 | | bookmark-based |
| `clear()` | RWMol.h:204-209 | | |
| public raw setters: `Atom::setIdx`, `Bond::setIdx`, `setBeginAtomIdx`, `setEndAtomIdx`, `setBeginAtom`, `setEndAtom`, `Bond::setOwningMol` | Atom.h:170-175; Bond.h:196-263 | | do not touch the graph |

What each edit invalidates:

| Operation | RingInfo | computed props | endpoint valence cache | bookmarks | StereoGroups | SubstanceGroups | conformers | bond stereo atoms |
|---|---|---|---|---|---|---|---|---|
| `addAtom` | kept | kept | — | ci_RIGHTMOST_ATOM | — | — | extended with (0,0,0) | — |
| `addBond(i,j,t)` | **kept (stale)** | kept | cleared | — | — | — | — | — |
| `addBond(Bond*)` | **kept (stale)** | kept | **not cleared** | — | — | — | — | — |
| `removeBond` | reset (613) | kept | **not cleared** | pruned (571-584) | bond removed (616) | referencing groups removed (615) | — | neighbors' cleared; CIS/TRANS→NONE (586-610) |
| `removeAtom` (default `clearProps=true`) | reset via `clearComputedProps(true)` (486) | cleared on mol, atoms, bonds | not cleared | pruned (375-388) | atom removed (481) | referencing groups removed (478) | position removed (397-407) | reindexed, or cleared if they name the atom (467-475) |
| `replaceAtom` | kept | kept | new atom's cache copied in | repointed (253-259) | repointed (262-274) | — | — | — |
| `replaceBond` | kept | kept | not cleared | repointed (318-324) | **not repointed** | removed only if `!keepSGroups` (313-315) | — | — |
| `insertMol` | **kept (stale)** | kept | — | — | merged | merged | see row above | offset |
| `commitBatchEdit` | reset (678) | cleared (681) | not cleared | pruned | pruned | pruned | compacted (900-913) | remapped (878-898) |

### Python

| Name | Binding | Notes |
|---|---|---|
| `RWMol.AddAtom(atom)` | nb/Mol.cpp:1102-1103; 137-140 | copies; `updateLabel=true`; returns index |
| `RWMol.AddBond(beginAtomIdx, endAtomIdx, order=UNSPECIFIED)` | 1099-1101; 133-136 | returns new bond count |
| `RWMol.RemoveAtom(idx)`, `RemoveBond(idx1, idx2)` | 1094-1097 | |
| `RWMol.ReplaceAtom(index, newAtom, updateLabel=False, preserveProps=False)` | 1104-1109 | docstring says `updateLabel` works |
| `RWMol.ReplaceBond(index, newBond, preserveProps=False, keepSGroups=True)` | 1110-1115 | |
| `RWMol.InsertMol(mol)` | 1122-1123 | |
| `RWMol.BeginBatchEdit`, `RollbackBatchEdit`, `CommitBatchEdit` | 1125-1129 | |
| `with RWMol(...) as m:` | 1090-1093; 172-188 | exception → rollback; otherwise commit |
| `RWMol.SetStereoGroups(list)` | 1119-1120; 151-166 | validates atom ownership only, not bonds |
| `EditableMol.AddAtom`, `AddBond`, `RemoveAtom`, `RemoveBond`, `ReplaceAtom`, `ReplaceBond` (no `keepSGroups`), `Begin/Rollback/CommitBatchEdit`, `GetMol` | nb/EditableMol.cpp:105-135 | |

Not exposed in Python: `createPartialBond` / `finishPartialBond`, active-atom API, `clear`, `addBond(Bond*)`, bookmarks.

---

## 4. Atom state

| State | Storage (Atom.h) | C++ accessors | Python (nb/Atom.cpp) | Kind |
|---|---|---|---|---|
| atomic number | `uint8 d_atomicNum` 418 | `getAtomicNum` / `setAtomicNum(int)` 146-148 | `GetAtomicNum` / `SetAtomicNum` 187-190 | stored; setter takes `int` |
| symbol | — | `getSymbol` 151 (Atom.cpp:266-274; dummies read prop `dummyLabel`) | `GetSymbol` 192 | computed |
| formal charge | `int8 d_formalCharge` 417 | `getFormalCharge` / `setFormalCharge(int)` 212-214 | `Get/SetFormalCharge` 233-234 | stored |
| isotope | `uint16 d_isotope` 427 | `getIsotope` 238 / `setIsotope(unsigned)` 236 (Atom.cpp:642) | `Get/SetIsotope` 250-251 | stored |
| mass | — | `getMass` 233 (Atom.cpp:644-655) | `GetMass` 249 | computed from the periodic table |
| explicit H count | `uint8 d_numExplicitHs` 416 | `get/setNumExplicitHs` 223-225 | `Get/SetNumExplicitHs` 245-246 | stored |
| noImplicit | `bool df_noImplicit` 415 | `get/setNoImplicit` 218-220 | `Get/SetNoImplicit` 236-243 | stored |
| explicit valence | `int8 d_explicitValence` 422 (−1 = unset) | `getValence(EXPLICIT)` 196; `calcExplicitValence(strict)` 360; deprecated `getExplicitValence` 199 | `GetValence(which)` 226, `GetExplicitValence` (deprecated) 221 | cache; needs `updatePropertyCache` |
| implicit valence / implicit Hs | `int8 d_implicitValence` 422 | `getValence(IMPLICIT)`, `getNumImplicitHs` 193, `calcImplicitValence` 367, deprecated `getImplicitValence` 202 | `GetNumImplicitHs` 218, `GetImplicitValence` (deprecated) 224 | cache |
| total valence / total Hs | — | `getTotalValence` 190; `getTotalNumHs(includeNeighbors=false)` 186 | `GetTotalValence` 228, `GetTotalNumHs` 210 | needs the cache |
| degree / total degree | graph | `getDegree` 178 (returns 0 with no owner, Atom.cpp:276-278); `getTotalDegree` 182 | `GetDegree`, `GetTotalDegree` 198-208 | total needs the cache |
| valence violation | — | `hasValenceViolation` 205 (Atom.cpp:657-697; needs owner; recomputes) | `HasValenceViolation` 230 | computed |
| radical electrons | `uint8` 423 | `get/setNumRadicalElectrons` 208-209 | 252-253 | stored; written by `assignRadicals`, only for `noImplicit` atoms (MolOps.cpp:437-443) |
| hybridization | `uint8 d_hybrid` 425; enum 93-103 | `get/setHybridization` 250-254 | 260-264 | stored; written by `setHybridization` |
| aromatic | `bool df_isAromatic` 414 | `get/setIsAromatic` 228-230 | 247-248 | stored; written by `setAromaticity`, `Kekulize`, `RWMol::addBond(AROMATIC)` |
| chiral tag | `uint8 d_chiralTag` 424; enum 106-117 | `set/getChiralTag` 241-247; `invertChirality` 243 | `Set/GetChiralTag`, `InvertChirality` 256-258 | stored; the non-tetrahedral permutation lives in prop `_chiralPermutation` (Atom.cpp:836-861) |
| atom map number | prop `molAtomMapNumber` | `setAtomMapNum(int, bool strict=true)` 382-391; `getAtomMapNum` 394-398 | `Get/SetAtomMapNum(mapno, strict=False)` 549-553 | stored as a prop |
| R label, alias, value, supplemental SMILES label | props | free functions Atom.h:440-458 | module functions `Get/SetAtomRLabel`, `Get/SetAtomAlias`, `Get/SetAtomValue`, `Get/SetSupplementalSmilesLabel` 603-631 | stored as props |
| CIP code / rank / possible | props `_CIPCode`, `_CIPRank`, `_ChiralityPossible` (types.h:71-77) | `getProp` | `GetProp` | written by perception (§7, §10) |
| monomer info | `AtomMonomerInfo *dp_monomerInfo` 432 (owning raw pointer) | `getMonomerInfo` ×2 369-370; `setMonomerInfo` takes ownership 372 (Atom.cpp:633-636) | `GetMonomerInfo`, `GetPDBResidueInfo` (`reference_internal`) 539-544; `SetMonomerInfo`, `SetPDBResidueInfo` (copy) 545-548 | stored; classes MonomerInfo.h:26-120 |
| macro atom info | `unique_ptr<MacroAtomInfo>` 433 | 374-379 | not bound | stored |
| index | `atomindex_t d_index` 428 (`uint32`, details.h:14) | `getIdx` / `setIdx` ×2 163-175 | `GetIdx` 195 | stored; setter public |
| owner | `ROMol *dp_mol` 431 | `hasOwningMol` 154; `getOwningMol() const → ROMol&` 157-160; `setOwningMol` protected 410-412 | `HasOwningMol`, `GetOwningMol` 266-270 | back pointer |
| query-ness | virtual | `hasQuery` 269, `getQueryType` 271, `setQuery` / `expandQuery` (PRECONDITION(0)) 274-282, `getQuery` (`nullptr`) 277, `Match` 316 | `HasQuery`, `GetQueryType`, `DescribeQuery`, `Match`, `GetSmarts` 254-317; `QueryAtom.ExpandQuery` / `SetQuery` 595-601 | QueryAtom overrides these (QueryAtom.h:94-129) |
| flags | `uint64 d_flags` 429 | `setFlags`, `getFlags` ×2 404-406 | not bound | scratch; not serialized |
| perturbation order | graph | `getPerturbationOrder` 342 | not bound | computed |
| π electrons, aromatic-ish | — | `numPiElectrons` 465, `isAromaticAtom` 463, `isEarlyAtom` 461 | module `GetNumPiElectrons` 632 | computed |

Python-only atom extras:

- `IsInRing`, `IsInRingSize` run `findSSSR` lazily (nb/Atom.cpp:78-90, 299-306).
- `UpdatePropertyCache`, `NeedsUpdatePropertyCache`, `ClearPropertyCache` (527-537).

Calls that need a prior perception step (`updatePropertyCache` or sanitization):

- `getNumImplicitHs` PRECONDITIONs `d_implicitValence > -1` unless `noImplicit`, whether or not the atom has an owner (Atom.cpp:300-309).
- `getValence` returns 0 with no owner, and otherwise PRECONDITIONs the cache (Atom.cpp:319-335).
- Callers that reach the cache include `getTotalNumHs`, `getTotalDegree`, `getTotalValence`, `ROMol::getNumAtoms(false)`, `getNumBonds(false)`, `setConjugation` (ConjugHybrid.cpp:163) and `setHybridization` (`getTotalDegree`, ConjugHybrid.cpp:186-209).

---

## 5. Bond state

| State | Storage (Bond.h) | C++ | Python (nb/Bond.cpp) | Kind |
|---|---|---|---|---|
| type | `uint8 d_bondType` 403; enum 57-81 (22 values) | `get/setBondType` 163-165; `getBondTypeAsDouble` 168; free `getTwiceBondType` 434 | `Get/SetBondType`, `GetBondTypeAsDouble` 114-120 | stored |
| valence contribution | — | virtual `getValenceContrib(const Atom*)` 175 (Bond.cpp:195-210: dative bonds count only toward the end atom; needs owner) | `GetValenceContrib` 139-145 | computed |
| direction | `uint8 d_dirTag` 404; enum 84-93 | `set/getBondDir` 314-316 | `Get/SetBondDir` 122-125 | stored |
| stereo | `uint8 d_stereo` 405; enum 96-107 | `setStereo` 330-336 (PRECONDITION: stereo atoms required for CIS/TRANS); `getStereo` 338 | `Get/SetStereo` 127-130 | stored |
| stereo atoms | `INT_VECT *dp_stereoAtoms` 397 (allocated lazily) | `setStereoAtoms(bgn, end)` 352 (PRECONDITION neighbors, Bond.cpp:255-268); `getStereoAtoms()` const 355-360 / non-const 362-367 | `GetStereoAtoms` (copy) 131; `SetStereoAtoms` 134 | stored, as indices |
| aromatic | `bool df_isAromatic` 401 | 178-180 | `Get/SetIsAromatic` 147-148 | stored |
| conjugated | `bool df_isConjugated` 402 | 183-185 | `Get/SetIsConjugated` 150-152 | stored; written by `setConjugation` |
| begin / end | `atomindex_t d_beginAtomIdx, d_endAtomIdx` 400 | `getBegin/EndAtomIdx` 223-230; `getOtherAtomIdx` 237 (POSTCONDITION on a bad index, Bond.cpp:88-97); `setBegin/EndAtomIdx` 244-250 (range-checked only when owned, Bond.cpp:99-111); `set/getBeginAtom`, `EndAtom`, `getOtherAtom` 257-282 (need owner) | `GetBeginAtomIdx`, `GetEndAtomIdx`, `GetOtherAtomIdx`, `GetBeginAtom`, `GetEndAtom`, `GetOtherAtom` 157-172 | stored |
| index | `d_index` 399 | `getIdx` / `setIdx` 209-216 | `GetIdx` 154 | stored |
| owner | `ROMol *dp_mol` 396 | `hasOwningMol` 188; `getOwningMol() const → ROMol&` 191-194; **`setOwningMol` public** 196-198 | `HasOwningMol`, `GetOwningMol` 108-112 | back pointer |
| query-ness | virtual 290-311 | `hasQuery`; `setQuery` / `getQuery` / `expandQuery` all PRECONDITION(0) in the base (`getQuery` throws, Bond.cpp:221-224); `Match` 311 (Bond.cpp:226-235) | `HasQuery`, `DescribeQuery`, `Match`, `GetSmarts` 174-198; `QueryBond.ExpandQuery` / `SetQuery` 382-388 | QueryBond overrides these (QueryBond.h:60-114) |
| atropisomer chirality | in `d_stereo` | `invertChirality` 201 (Bond.cpp:326-339) | `InvertChirality` 137 | stored |
| macro bond info | `unique_ptr<MacroBondInfo>` 398 | 369-374 | not bound | stored |
| flags | `uint64 d_flags` 406 | 387-389 | not bound | scratch |
| `updatePropertyCache(strict)` | — | no-op 381 | not bound | — |

Free helpers: `isDative` ×2, `canSetDoubleBondStereo`, `canHaveDirection` (Bond.h:411-430).

Python-only bond extras: `IsInRing`, `IsInRingSize` run `findSSSR` lazily (nb/Bond.cpp:66-79, 182-188).

---

## 6. Properties

### C++: `RDProps` (base of ROMol, Atom, Bond, Conformer, SubstanceGroup)

- **Storage:** `mutable Dict d_props` (RDProps.h:19-22). The source comment calls the `mutable` "a quirk of history".
- **Public API (12 declarations)**
  - `clear()` 37
  - `getDict()` ×2 39-40
  - `getPropList(includePrivate=true, includeComputed=true)` 48-66
  - `setProp<T>(key, val, computed=false) const` 79-93
  - `getProp<T>(key, T&)` 113; `getProp<T>(key)` 119
  - `getPropIfPresent` 127
  - `hasProp` 132
  - `clearProp const` 142-152
  - `clearComputedProps() const` 155-165
  - `updateProps(source, preserveExisting=false)` 173-175
- **Computed-prop bookkeeping**
  - A `vector<string>` of names is stored inside the same dict under `__computedProps` (types.h:51).
  - `setProp`, `clearProp` and `clearComputedProps` each read, copy and write back that vector (RDProps.h:84-91, 143-150, 156-164).
  - "Private" means the name starts with `_` (RDProps.h:59). Reserved names: `common_properties` (types.h:71-159).
- **`Dict`**
  - `vector<Pair>` with linear search (Dict.h:166-185, 221-228, 298-305).
  - Empty key → `ValueErrorException` (294-296); missing key → `KeyErrorException` (227).
  - A `std::string` target converts any stored type through `rdvalue_tostring` (231-239).
  - Arithmetic targets parse string values through `lexical_cast` under a `LocaleSwitcher` (RDValue.h:268-292).
  - A type mismatch throws `std::bad_any_cast` or `boost::bad_lexical_cast` (RDProps.h:103-108).
- **ROMol-specific**
  - `setName(const std::string&) const` / `getName()` wrap `_Name` (ROMol.h:955-958; ROMol.cpp:591-599).
  - `clearComputedProps(bool includeRings=true) const` (ROMol.h:960; ROMol.cpp:601-616) resets RingInfo, then clears computed props on the mol, atoms and bonds. Conformers and SubstanceGroups are not touched.

### Python property families

| Method | Mol (nb/Mol.cpp) | Atom (nb/Atom.cpp) | Bond (nb/Bond.cpp) | Conformer (nb/Conformer.cpp) | SubstanceGroup (nb/SubstanceGroup.cpp) |
|---|---|---|---|---|---|
| `SetProp(key, val[, computed])` | 748, has `computed` | 320, **no `computed`** | 201, **no `computed`** | 147, has `computed` | 211, has `computed` |
| `Set{Double,Int,Unsigned,Bool}Prop` | 759-802 | 356-457 | 225-299 | 156-192 | 218-241 |
| `Set/GetExplicitBitVectProp` | — | 481-501 | — | — | — |
| `GetProp(key, autoConvert=False[, default])` | 818, 831 | 329, 342 | 207, 217 | 199, 209 | 248, 262 |
| `Get{Double,Int,Unsigned,Bool}Prop(key[, default])` | 844-931 | 374-479 | 237-314 | 217-280 | 276-315; also `GetUnsignedVectProp`, `GetStringVectProp` 316-326 |
| `HasProp` | 804 | 504 | 316 | 193 | 242 |
| `ClearProp` | 932 | 512 | 322 | 281 | 340 |
| `ClearComputedProps` | 986 (`includeRings=True`) | — | — | 286 | — |
| `GetPropNames(includePrivate=False, includeComputed=False)` | 939 | 519 | 328 | 288 | 327 |
| `GetPropsAsDict` defaults `(includePrivate, includeComputed)` | **False, False** (952) | **True, True** (523) | **True, True** (332) | **False, False** (300) | **True, True** (333) |
| `SetName` / `GetName` | 810-817 | — | — | — | — |

Shared helpers (nb/props.hpp):

- `HasProp` returns `int` (328-333; nb/Atom.cpp:65-69; nb/Bond.cpp:50-53).
- `GetProp<T>` maps `KeyErrorException` to `KeyError` and other errors to `ValueError` (170-183).
- `GetPyPropImpl` implements `autoConvert` (227-314).
- `GetPropsAsDict` skips `AnyTag` values and reaches `UNDER_CONSTRUCTION` on unknown tags (77-167).

---

## 7. Derived and cached state

| State | Where stored | Who computes | Who invalidates | Staleness detection |
|---|---|---|---|---|
| ring info (SSSR, symmetrized, fast) | `RingInfo *dp_ringInfo` ROMol.h:1020; `getRingInfo() const → RingInfo*` 773 | `findSSSR` (FindRings.cpp:958-), `symmetrizeSSSR` (1085-1114), `fastFindRings` (1116-1126), `findRingFamilies` (1128-), `RingInfo::addRing` | `RingInfo::reset` (RingInfo.cpp:303-316) via `clearComputedProps(true)`, `removeBond` (RWMol.cpp:613), `commitBatchEdit` (678), start of `sanitizeMol` (MolOps.cpp:612) | only `isInitialized()` and `getRingType()` (RingInfo.h:49-70). Per-index queries return 0 or empty for unknown indices (RingInfo.cpp:16-58, 129-145) |
| lazy ring filling | — | Python `Atom/Bond.IsInRing` (nb/Atom.cpp:78-90; nb/Bond.cpp:66-79); `legacyStereoPerception` → `fastFindRings` (Chirality.cpp:2289-2291); `Kekulize` (Kekulize.cpp:664-666); `setAromaticity` (Aromaticity.cpp:1146-1150); atropisomer cleanup `checkBond` (MolOps.cpp:547-549); canonicalization (Canon.cpp:1432-1434) | | consumers accept any initialized ring type |
| ring families (URF) | public member `dp_urfData` (RingInfo.h:335-336); `d_atomRingFamilies` | `findRingFamilies` | `reset(true)`, `resetRingFamilies` (RingInfo.cpp:293-297) | `areRingFamiliesInitialized()`; `findRingFamilies` returns early, ignoring `includeDativeBonds` / `includeHydrogenBonds` (FindRings.cpp:1130-1134) |
| fused-ring caches | `d_fusedRings`, `d_numFusedBonds` (RingInfo.h:332-333) | `initFusedRings` (RingInfo.cpp:244-269), `numFusedBonds` (211-223), lazily from non-const queries | **never cleared**, not even by `reset` (303-316) | compares sizes only (245, 214) |
| valence cache (explicit/implicit) | `Atom::d_explicitValence`, `d_implicitValence` (Atom.h:422) | `Atom::updatePropertyCache` / `calc*Valence` (Atom.cpp:618-631, 739-742); `ROMol::updatePropertyCache` (ROMol.cpp:618-625); `sanitizeMol` | `clearPropertyCache` (Atom.cpp:749-752; ROMol.cpp:637-641); `RWMol::addBond(i,j)` on both endpoints | `needsUpdatePropertyCache` only asks whether it was ever computed (Atom.cpp:744-747; ROMol.cpp:627-635). Setters do not invalidate. Atom copies carry the cache along |
| conjugation | `Bond::df_isConjugated` | `setConjugation(ROMol&)` (ConjugHybrid.cpp:151-174) | nothing | none |
| hybridization | `Atom::d_hybrid` | `setHybridization(ROMol&)` (ConjugHybrid.cpp:176-) | nothing | none (default `UNSPECIFIED`) |
| aromaticity | `Atom/Bond::df_isAromatic` | `setAromaticity` (Aromaticity.cpp:1139-1180); parsers; cleared by `Kekulize(markAtomsBonds)` | nothing automatic | none |
| radicals | `Atom::d_numRadicalElectrons` | `assignRadicals` (MolOps.cpp:437-) | nothing | none |
| legacy CIP labels | atom prop `_CIPCode`, **not marked computed** (Chirality.cpp:1799) | `assignStereochemistry` when the legacy switch is on (via `assignAtomChiralCodes`, Chirality.cpp:1724-) | `assignStereochemistry(cleanIt=true)` (2300-2325), `removeStereochemistry` (3837-3859) | none; survives `clearComputedProps` |
| new CIP labels | `_CIPCode` (e.g. CIPLabeler/configs/AtropisomerBond.cpp:52); mol `_CIPComputed` (CIPLabeler.cpp:376) | `CIPLabeler::assignCIPLabels` (CIPLabeler.h:68-69, 97-100) | `_CIPComputed` cleared at CIPLabeler.cpp:323 | flag only |
| CIP ranks | atom prop `_CIPRank`, computed (Chirality.cpp:1327) | `Chirality::assignAtomCIPRanks(const ROMol&, UINT_VECT&)` (1308-1329) | `clearComputedProps`, `removeStereochemistry` | none |
| stereo perceived | mol prop `_StereochemDone`, computed (Chirality.cpp:2922) | `MolOps::assignStereochemistry` (2907-2923) | `clearComputedProps`; explicit clears in `assignChiralTypesFrom3D` (3407-3408), `assignChiralTypesFromMolParity` (3526-3527), `removeStereochemistry` (3838-3839) | flag presence skips the run unless `force` (2909-2911). Set **non-computed** by MolInterchange/Parser.cpp:926, 990 and ChemTransforms/MolZip.cpp:654 |
| `_ChiralityPossible` | atom prop | `assignStereochemistry(flagPossibleStereoCenters=true)` | `cleanIt` | none |
| matrices | mol computed props `DistanceMatrix[BO][AtomWts]`, `…_Paths`, `AdjacencyMatrix[BO]`, `_3DDistanceMatrix_Conf<id>[_AtomWeights]` (Matrices.cpp:166-232, 287-331, 390-432) | `getDistanceMat`, `getAdjacencyMatrix`, `get3DDistanceMat`, all taking `const ROMol&` | `clearComputedProps`; `force=true` | adjacency key ignores `emptyVal` and `bondsToUse` (Matrices.cpp:297-304) |
| batch-edit pending deletions | `unique_ptr<dynamic_bitset> dp_delAtoms`, `dp_delBonds` (ROMol.h:1024-1025) | `beginBatchEdit`, `removeAtom`, `removeBond` | `commitBatchEdit`, `rollbackBatchEdit` | copied by the copy constructor (ROMol.cpp:121-130) |

---

## 8. Sanitization and perception

| C++ (MolOps.h) | Mutates | Impl | Python (nb/MolOps.cpp) | Python behavior |
|---|---|---|---|---|
| `sanitizeMol(RWMol&, unsigned &operationThatFailed, unsigned ops=SANITIZE_ALL)` 573; `sanitizeMol(RWMol&)` 577 | yes | MolOps.cpp:605-719 | `SanitizeMol(mol, sanitizeOps=ALL, catchErrors=False)` 1259 (helper 421-438) | in place, on a `Mol` via `static_cast<RWMol&>`; returns `SanitizeFlags`; `catchErrors` swallows every exception (432) |
| `SanitizeFlags` BETTER_ENUM, 14 values, `ALL = 0xFFFFFFF` | — | MolOps.h:520-535 | enum `SanitizeFlags` 1179-1196 | |
| `detectChemistryProblems(const ROMol&, ops)` 601-603 | no; works on a copy (MolOps.cpp:721-723) | | `DetectChemistryProblems` 3204 (helper 915-923) | tuple of exception objects |
| `Kekulize(RWMol&, markAtomsBonds=true, canonical=true, maxBackTracks=100)` 764-766 | yes | Kekulize.cpp | `Kekulize(mol, clearAromaticFlags=False, canonical=True)` 1824 (440-444) | in place; parameter renamed, default inverted |
| `KekulizeIfPossible(…) → bool` 794-797 | yes, only on success | | `KekulizeIfPossible` 1857 (445-449) | in place; **returns `None`** |
| `setAromaticity(RWMol&, model=DEFAULT, int(*func)(RWMol&)=nullptr) → int` 656-658; models 616-623 | yes | Aromaticity.cpp:1139-1180 | `SetAromaticity(mol, model=AROMATICITY_DEFAULT)` 1938 (461-464); enum 1916-1923 | in place; return value dropped; `AROMATICITY_CUSTOM` is exposed but no function can be passed |
| `setMMFFAromaticity(RWMol&)` 626 | yes | | not in rdmolops | |
| `cleanUp(RWMol&)` 680 | yes | MolOps.cpp:315- | `Cleanup` 1872 | in place |
| `cleanUpOrganometallics(RWMol&)` 697 | yes | MolOps.cpp:331- | `CleanupOrganometallics` 1913 | in place |
| `assignRadicals(RWMol&)` 700 | yes | MolOps.cpp:437- | `AssignRadicals` 1976 | in place |
| `adjustHs(RWMol&)` 717 | yes | MolOps.cpp:375-435 | not bound | |
| `setConjugation(ROMol&)` 800 | yes | ConjugHybrid.cpp:151 | `SetConjugation` 1952 | in place |
| `setHybridization(ROMol&)` 803 | yes | ConjugHybrid.cpp:176 | `SetHybridization` 1964 | in place |
| `cleanupChirality(RWMol&)` 1124 | yes | Chirality.cpp:3051 | `CleanupChirality` 1885 | in place |
| `cleanupAtropisomers(RWMol&[, Hybridizations&])` 1128, 1131 | yes | | `CleanupAtropisomers` 1898 | in place |
| `findSSSR(const ROMol&, VECT_INT_VECT &res, dative=false, hbonds=false)` 867; `findSSSR(const ROMol&, VECT_INT_VECT *res=nullptr, …)` 872 | yes (RingInfo, through `const`) | FindRings.cpp:958- | `GetSSSR(mol, includeDativeBonds=False, includeHydrogenBonds=False)` 1276 (323-328) | returns rings and fills RingInfo |
| `fastFindRings(const ROMol&)` 884 | yes, through `const` | FindRings.cpp:1116 | `FastFindRings` 1340 | in place |
| `findRingFamilies(const ROMol&, …)` 886 | yes, through `const` | FindRings.cpp:1128 | `FindRingFamilies` 1353 | in place |
| `symmetrizeSSSR`: 4 overloads 928-959, enum `SymmetrizeSSSRAlgorithm` 890-894 | yes | FindRings.cpp:1085-1114 | `GetSymmSSSR(mol, includeDativeBonds, includeHydrogenBonds, algorithm, recalcSSSR)` 1306 (501-509); enum 1279-1282 | returns rings and fills RingInfo; argument order differs from C++ |
| `set/getUseLegacyRingFinding` 815-818 | global | FindRings.cpp:936-946 | 1284-1288 | |

Actual `sanitizeMol` order (MolOps.cpp:609-719):

1. `clearComputedProps()`, which also resets ring info
2. `cleanUp`
3. `cleanUpOrganometallics`
4. `updatePropertyCache(strict = PROPERTIES bit)`
5. Kekulize, via the internal `kekulizeForSanitize`
6. `symmetrizeSSSR(recalc=false)`, whose algorithm comes from the legacy-ring env switch
7. `assignRadicals`
8. `setAromaticity`
9. `setConjugation`
10. `setHybridization`
11. `cleanupAtropisomers`
12. `cleanupChirality`
13. `adjustHs`
14. `updatePropertyCache(true)`

The header's documented sequence (MolOps.h:541-552) leaves out steps 1, 3 and 11.

---

## 9. Hydrogens

H-count state, recapped from §4:

- `d_numExplicitHs` and `df_noImplicit` are stored fields.
- Implicit Hs come from the implicit-valence cache.
- `getTotalNumHs(includeNeighbors)` adds explicit H neighbors in the graph only when `includeNeighbors` is true (Atom.cpp:289-298).

| C++ (MolOps.h) | In place / new | Python (nb/MolOps.cpp) |
|---|---|---|
| `struct AddHsParameters {explicitOnly, addCoords, addResidueInfo, skipQueries}`, all `false` | — | `AddHsParameters` 1359-1368 |
| `addHs(RWMol&, const AddHsParameters&, const UINT_VECT *onlyOnAtoms=nullptr)` 204 (AddHs.cpp:532) | in place | — |
| `inline ROMol *addHs(const ROMol&, explicitOnly=false, addCoords=false, onlyOnAtoms=nullptr, addResidueInfo=false)` 229-237 | new; caller deletes | `AddHs(mol, explicitOnly=False, addCoords=False, onlyOnAtoms=None, addResidueInfo=False)` 1425 → new |
| `inline void addHs(RWMol&, bool, bool, const UINT_VECT*, bool)` 240-245 | in place | `AddHs(mol, params, onlyOnAtoms=None)` 1392 → new (306-315) |
| `removeHs(const ROMol&, bool implicitOnly, bool updateExplicitCount=false, bool sanitize=true)` `[[deprecated]]` 292-295 (AddHs.cpp:1100) | new | `RemoveHs(mol, implicitOnly=False, updateExplicitCount=False, sanitize=True)` 1467-1471 → new; binds the deprecated overload with no warning |
| `removeHs(RWMol&, bool implicitOnly, …)` `[[deprecated]]` 298-301 (AddHs.cpp:1093-1099) | in place | — |
| `struct RemoveHsParameters`: 16 fields (MolOps.h:302-331); defaults `removeWithWedgedBond`, `removeMapped`, `removeInSGroups`, `showWarnings`, `removeNonimplicit` = true; the rest false | — | `RemoveHsParameters` 1478-1524, all 16 fields |
| `removeHs(RWMol&, const RemoveHsParameters& = {}, bool sanitize=true)` 335-337 (AddHs.cpp:1008-1081; calls `sanitizeMol` when Hs were removed, 1064-1066) | in place | — |
| `removeHs(const ROMol&, const RemoveHsParameters& = {}, bool sanitize=true)` 340-342 | new | `RemoveHs(mol, params, sanitize=True)` 1525-1532 → new |
| `removeAllHs(RWMol&, bool sanitize=true)` 345 / `removeAllHs(const ROMol&, bool)` 348-349 | in place / new | `RemoveAllHs(mol, sanitize=True)` 1533 → new |
| `mergeQueryHs(const ROMol&, mergeUnmappedOnly=false, mergeIsotopes=false)` 376-378 / `mergeQueryHs(RWMol&, …)` 381-383 | new / in place | `MergeQueryHs` 1538 → new |
| `hasQueryHs(const ROMol&) → pair<bool,bool>` 396 | — | `HasQueryHs` 1552 → tuple |
| `needsHs(const ROMol&)` 1282 | — | `NeedsHs` 3316 |
| `setTerminalAtomCoords(ROMol&, idx, otherIdx)` 259-260 | in place | `SetTerminalAtomCoords` 1327 |

---

## 10. Stereo

**Atom chiral tags.** `ChiralType`, 9 values (Atom.h:106-117).

- Tetrahedral sense is relative to bond order, with implicit H last (Atom.h:59-73).
- Non-tetrahedral permutation lives in prop `_chiralPermutation` (Atom.cpp:836-861).
- `getPerturbationOrder` (Atom.h:342), `invertChirality` (243).

**Bond stereo.** `BondStereo`, 8 values (Bond.h:96-107), plus stereo atoms (§5) and `BondDir`.

**StereoGroups** (StereoGroup.h)

- `enum class StereoGroupType {ABSOLUTE, OR, AND}` (31-35).
- `class StereoGroup` (44-90): holds `vector<Atom*>`, `vector<Bond*>`, readId, writeId. `operator==` ignores the ids (78-85).
- Free functions: `removeAtomFromGroups`, `removeBondFromGroups`, `removeGroupsWithAtom(s)`, `assignStereoGroupIds`, `forwardStereoGroupIds` (91-108).
- Mol access: `getStereoGroups() const&` (ROMol.h:992-994) and `setStereoGroups(vector)` (1003). The setter merges multiple ABS groups and moves the result to the end (ROMol.cpp:434-464).

**Perception entry points** (MolOps.h; implementations in Chirality.cpp)

| C++ | Line | Python (nb/MolOps.cpp unless noted) |
|---|---|---|
| `assignStereochemistry(ROMol&, cleanIt=false, force=false, flagPossibleStereoCenters=false)` | 1229-1231 (Chirality.cpp:2907-2923) | `AssignStereochemistry` 2287 |
| `assignStereochemistryFrom3D(ROMol&, confId=-1, replaceExistingTags=true)` | 1162 (Chirality.cpp:3768-3783: silently returns if the conformer is not 3D; forces perception) | `AssignStereochemistryFrom3D` 2320 |
| `assignChiralTypesFrom3D(ROMol&, confId=-1, replace=true)` | 1148 | `AssignAtomChiralTagsFromStructure` 2362 |
| `assignChiralTypesFromBondDirs(ROMol&, confId=-1, replace=true)` | 1173 | `AssignChiralTypesFromBondDirs` 2306 |
| `assignChiralTypesFromMolParity(ROMol&, replace=true)` | 1272 | `AssignAtomChiralTagsFromMolParity` 2377 |
| `detectBondStereochemistry(ROMol&, confId=-1)` (deprecated) | 1178 | `DetectBondStereochemistry` 1216 |
| `RDKit::DetectBondStereoChemistry(ROMol&, const Conformer*)` | FileParsers/MolFileStereochem.h:22 | `DetectBondStereoChemistry` 1209 |
| `setDoubleBondNeighborDirections(ROMol&, const Conformer* = nullptr)` | 1181 | `SetDoubleBondNeighborDirections` 1226 |
| `setBondStereoFromDirections(ROMol&)` | 1198 | `SetBondStereoFromDirections` 1237 |
| `clearSingleBondDirFlags`, `clearAllBondDirFlags`, `clearDirFlags` | 1185-1194 | not bound |
| `removeStereochemistry(ROMol&)` | 1238 | `RemoveStereochemistry` 2343 |
| `findPotentialStereoBonds(ROMol&, cleanIt=false)` | 1263 | `FindPotentialStereoBonds` 2336 |
| `cleanupChirality(RWMol&)` | 1124 | `CleanupChirality` 1885 |
| `Chirality::findPotentialStereo(ROMol&, bool cleanIt, bool flagPossible=true)` / `(const ROMol&)` | Chirality.h:147-151 | `FindPotentialStereo(mol, cleanIt=False, flagPossible=True)` nb/ChiralityOps.cpp:23-30 |
| `Chirality::cleanupStereoGroups(ROMol&)` | Chirality.h:154 | `CleanupStereoGroups` nb/ChiralityOps.cpp:31 |
| `Chirality::assignAtomCIPRanks(const ROMol&, UINT_VECT&)` | Chirality.h:70-71 | `ComputeAtomCIPRanks` 2291 |
| `Chirality::assignLegacyCIPLabels(ROMol&, flag=false)` | Chirality.h:157-158 | not bound |
| `Chirality::findMesoCenters` (425-427), `simplifyEnhancedStereo` (417-418), `addStereoAnnotations` (400-403), `wedgeMolBonds` (325-327), `wedgeBond` (328-329), `reapplyMolBlockWedging` (344-346), `removeNonExplicit3DChirality` (49) | Chirality.h, lines in first column | `FindMesoCenters` 2391, `SimplifyEnhancedStereo` 3301, `AddStereoAnnotations` 3276, `WedgeMolBonds` 2630, `WedgeBond` 2718, `ReapplyMolBlockWedging` 2647, `RemoveNonExplicit3DChirality` 2659 |
| `CIPLabeler::assignCIPLabels(ROMol&, maxIter=0)` and the atoms/bonds overload | CIPLabeler/CIPLabeler.h:68-69, 97-100 | `rdCIPLabeler.AssignCIPLabels` (CIPLabeler/nbWrap/rdCIPLabeler.cpp:67) |
| `StereoType`, `StereoDescriptor` (7 values), `StereoSpecified`, `struct StereoInfo` | Chirality.h:92-136 | nb/Chirality.cpp:21-63; `StereoDescriptor` binds only 5 values |

Python StereoGroup:

- `StereoGroup.GetGroupType`, `GetAtoms` (tuple), `GetBonds` (tuple), `GetReadId`, `GetWriteId`, `SetWriteId` (nb/StereoGroup.cpp:88-103). The class has no constructor.
- `CreateStereoGroup(stereoGroupType, mol, atomIds=[], bondIds=[], readId=0)`: `take_ownership`, `keep_alive<0,2>` (105-109).
- `ForwardStereoGroupIds` (111); `Mol.GetStereoGroups` (nb/Mol.cpp:955-958); `RWMol.SetStereoGroups` (1119).

### Legacy vs new perception switch

- Environment variable `RDK_USE_LEGACY_STEREO_PERCEPTION`, default **true** (Chirality.h:36, 39-41).
- `setUseLegacyStereoPerception` writes it with `setenv`; `getUseLegacyStereoPerception` reads it with `getenv` on every call (Chirality.cpp:848-858; utils.h:96-106).
- Read sites: `assignStereochemistry` (Chirality.cpp:2917-2921) and `canBeStereoBond` (2623-, read at 2656), which `shouldBeACrossedBond` calls.
- Only the legacy path writes `_CIPCode` during `assignStereochemistry`. The new path (`stereoPerception`, 2580-2817) only clears it. The documented behavior is at MolOps.h:1202-1208.
- Python `Set/GetUseLegacyStereoPerception` (nb/MolOps.cpp:3220-3226).
- A second switch, `RDK_ENABLE_NONTETRAHEDRAL_STEREO` (default true), is read in FindStereo.cpp:125, 327, 579 and Chirality.cpp:3411. Python `Set/GetAllowNontetrahedralChirality` (3213-3219).

---

## 11. Conformers

- **Ownership**
  - `ROMol` holds `std::list<boost::shared_ptr<Conformer>> d_confs` (ROMol.h:381, 1021).
  - `addConformer(Conformer*, bool assignId=false)` always takes ownership. It PRECONDITIONs that atom counts match, and `assignId` picks max+1 (ROMol.cpp:678-694).
- **Access**
  - `getConformer(int id=-1)`, const and non-const (ROMol.h:740-744): a negative id returns the first conformer; a missing id throws `ConformerException` (ROMol.cpp:643-667).
  - `removeConformer(unsigned)` is a silent no-op when the id is absent (669-676).
  - `clearConformers` (750), `getNumConformers` (764), and four iterators (939-947).
- **Conformer members** (Conformer.h)
  - `resize`, `reserve` (82-85).
  - `hasOwningMol` (88); `getOwningMol() const → ROMol&`, which PRECONDITIONs (91-94).
  - `getPositions` const / non-const (97-100). Only the const one checks `numAtoms == mol.numAtoms`, and its precondition message is empty (Conformer.cpp:22-29).
  - `getAtomPos` ×4 (103-116; URANGE_CHECK, Conformer.cpp:31-45).
  - `setAtomPos` ×2 (119-132): grows the position vector when `atomId >= size`.
  - `getId` / `setId` (134-137); `getNumAtoms` (140).
  - `is3D` / `set3D` (143-144); `df_is3D` defaults to **true** (154).
  - `setOwningMol` is protected (147-151).
  - Free `hasNonZeroZCoords` (166-174). `ConformerException` (26-38).
- **Python** (nb/Conformer.cpp)
  - `GetNumAtoms` 112; `HasOwningMol`, `GetOwningMol` (`reference_internal`) 115-118; `GetId`, `SetId` 120-122.
  - `GetAtomPosition` returns a `Point3D` copy (124-129).
  - `GetPositions` returns a numpy copy (130; 29-50).
  - `SetPositions` takes a C-contiguous float64 array, 2D or 3D (132; 52-82).
  - `SetAtomPosition(aid, loc)` from a sequence, and `SetAtomPosition(atomId, Point3D)` (134-140).
  - `Set3D`, `Is3D` (142-144), plus the property family (§6).
  - `Mol.AddConformer(conf, assignId=False)` **copies** the conformer (nb/Mol.cpp:513-522). `Mol.GetConformer(id=-1)` returns `reference_internal` (524-527). `RemoveConformer`, `RemoveAllConformers` (535-539).

---

## 12. Substance groups and other attachments

**SubstanceGroup** (SubstanceGroup.h)

- Class (54-238) holds indices: `d_atoms`, `d_patoms`, `d_bonds` (231-233), plus brackets, CStates and attach points (235-237).
- `SubstanceGroup(ROMol*, type)` sets only the group's side of ownership (93-96).
- `getOwningMol` (133-136); `getIsValid` / `setIsValid` (140-143); `getIndexInMol` searches the owner's vector by address (147; SubstanceGroup.cpp:43-60).
- `add*WithIdx` / `add*WithBookmark` (150-155); `remove*WithIdx`, documented as able to leave invalid state (157-161).
- `adjustToRemovedAtom` / `Bond`, `includesAtom` / `Bond` (195-211); `setOwningMol` is public (217).
- Free functions:
  - `getSubstanceGroups(ROMol&)` returns a mutable reference to the vector; it is a friend of ROMol (274-277; ROMol.h:1027-1031).
  - `addSubstanceGroup(ROMol&, SubstanceGroup)` takes a copy and `push_back`s it (284-285; SubstanceGroup.cpp:363-372).
  - `removeSubstanceGroupsReferencingAtom` / `Bond` (292-300).
  - `SubstanceGroupChecks` (240-269).
- `SubstanceGroupException : std::runtime_error` (35-42).
- Python:
  - Class methods at nb/SubstanceGroup.cpp:155-345.
  - Module functions `GetMolSubstanceGroups` (copy, 349-351), `GetMolSubstanceGroupWithIdx` (352-355), `ClearMolSubstanceGroups` (356), `CreateMolSubstanceGroup` (358-361), `CreateMolDataSubstanceGroup` (362-366), `AddMolSubstanceGroup` (367-370). The last four return pointers into the vector with `rv_policy::reference` + `keep_alive<0,1>`.

**Bookmarks** (ROMol.h:683-731)

- Atom and bond bookmarks are `map<int, list<Atom*|Bond*>>`.
- `get*Bookmarks()` returns a raw pointer to the internal map (708, 731).
- `getUnique*WithBookmark` does not check uniqueness (ROMol.cpp:236-240, 258-262).
- Python: not bound, although the `Mol` docstring advertises them (nb/Mol.cpp:197-199).

**Other attachments:** see §4 for `MonomerInfo` / `MacroAtomInfo` and §5 for `MacroBondInfo`. Ring info (§7), props (§6) and batch-edit bitsets (§7) also hang off the mol.

### Reference models in one mol

- `StereoGroup` and bookmarks hold **pointers**.
- `SubstanceGroup`, `Bond` stereo atoms and `RingInfo` hold **indices**.
- Haptic-bond endpoints are a **string property** `_MolFileBondEndPts` that `removeAtom` reparses (RWMol.cpp:418-456).

---

## 13. Fragments and whole-mol transforms

| C++ | Location | Result | Python (nb/MolOps.cpp) |
|---|---|---|---|
| `getMolFrags(const ROMol&, vector<int> &mapping)` | MolOps.h:72-73 | count + out-param | `GetMolFrags(mol, asMols=False, sanitizeFrags=True, frags=None, fragsMolAtomMapping=None)` 2203-2205 (597-644): tuple of index tuples, or tuple of new Mols; `frags` and `fragsMolAtomMapping` are lists the binding **appends to**, honored only when `asMols=True` |
| `getMolFrags(const ROMol&, vector<vector<int>>&)` | 85-86 | count + out-param | ″ |
| `getMolFrags(const ROMol&, vector<unique_ptr<ROMol>>&, sanitizeFrags=true, vector<int>* =nullptr, vector<vector<int>>* =nullptr, copyConformers=true)` | 108-112 | count + out-params | ″ |
| `getMolFrags(const ROMol&, sanitizeFrags=true, …) → vector<shared_ptr<ROMol>>` | 132-136 | new mols | ″ |
| `getMolFragsWithQuery<T>` ×2 | 153-180 | map, or count + out-param | `SplitMolByPDBResidues` 2219, `SplitMolByPDBChainId` 2234 (226-265) |
| `renumberAtoms(const ROMol&, const vector<unsigned>&) → ROMol*` | 511-512 (Renumber.cpp:19) | new; caller owns | `RenumberAtoms(mol, newOrder)` 3075 (196-206) |
| `combineMols(const ROMol&, const ROMol&, Point3D offset) → ROMol*` | ChemTransforms/ChemTransforms.h:205-207 | new | `CombineMols(mol1, mol2, offset=None)` 1598-1606 |
| `deleteSubstructs(const ROMol&, const ROMol&, onlyFrags=false, useChirality=false) → ROMol*` | ChemTransforms.h:38-41 | new | `DeleteSubstructs` 1590 |
| `replaceSubstructs`, `replaceSidechains`, `replaceCore` ×2, `MurckoDecompose` | ChemTransforms.h:81, 100, 136, 174, 190 | new | `ReplaceSubstructs` 1652, `ReplaceSidechains` 2750, `ReplaceCore` 2797 / 2882, `MurckoDecompose` 1594 |
| `MolFragmenter::fragmentOnBonds` ×2 | ChemTransforms/MolFragmenter.h:54, 62 | new | `FragmentOnBonds` 2915, `FragmentOnSomeBonds` 2920, `FragmentOnBRICSBonds` 2891 |
| `copyMolSubset(…) → unique_ptr<RWMol>` | Subset.h:92 | new | `CopyMolSubset` ×4 3370-3392 |
| `adjustQueryProperties(const ROMol&, const AdjustQueryParameters* = nullptr) → ROMol*` / `(RWMol&, …)` | MolOps.h:489-494 | new / in place | `AdjustQueryProperties` 3193 → new |
| `hapticBondsToDative(const ROMol&) → ROMol*` / `(RWMol&)` | 1296-1299 | new / in place | `HapticBondsToDative` 1992 → new |
| `dativeBondsToHaptic(const ROMol&) → ROMol*` / `(RWMol&)` | 1312-1315 | new / in place | `DativeBondsToHaptic` 2008 → new |
| `expandAttachmentPoints(RWMol&, …)`, `collapseAttachmentPoints(RWMol&, …)` | 1380-1399 | in place | `ExpandAttachmentPoints` 3248-3249, `CollapseAttachmentPoints` 3261-3262: in place via `static_cast<RWMol&>` (1097-1105) |
| `getAttachmentPointLabelNumber`, `isMarkedAttachmentPoint` | 1414-1427 | — | Boost only (bp/MolOps.cpp:3418, 3432) |
| `getFormalCharge`, `getShortestPath`, `getAvgMolWt`, `getExactMolWt`, `getMolFormula`, `getNumAtomsWithDistinctProperty`, `countAtomElec`, `atomHasConjugatedBond` | 56, 1084, 1326-1354, 1278, 53, 59 | — | `GetFormalCharge` 2246, `GetShortestPath` 2259, `CountAtomElec` 3319, `AtomHasConjugatedBond` 3321 |
| `getAdjacencyMatrix`, `getDistanceMat` ×2, `get3DDistanceMat` | 984-1066 | `double*` (ownership varies, H-D7) | `GetAdjacencyMatrix` 1787, `GetDistanceMatrix` 1734, `Get3DDistanceMatrix` 1761: numpy copies |

---

## 14. Error model

### C++ exception types

- `MolSanitizeException : std::exception` (SanitException.h:28-40), with:
  - `AtomSanitizeException` (45-60) → `AtomValenceException` (65-79), `AtomKekulizeException` (81-95)
  - `KekulizeException` (97-113)
  - virtual `getType()` and `copy()` on each.
- `ValueErrorException`, `IndexErrorException`, `KeyErrorException`, all `: std::runtime_error` (RDGeneral/Exceptions.h:21, 41, 57).
- `Invar::Invariant : std::runtime_error` (Invariant.h:50). It is thrown by `PRECONDITION`, `POSTCONDITION`, `CHECK_INVARIANT`, `URANGE_CHECK` and `UNDER_CONSTRUCTION`, each of which **logs to `rdErrorLog` before throwing** (100-160).
  - Build-dependent: under `INVARIANT_ASSERT_METHOD` these become `assert`; under `INVARIANT_SILENT_METHOD` they compile to nothing (34-46, 161-180).
- Others:
  - `ConformerException : std::exception` (Conformer.h:26-38)
  - `SubstanceGroupException : std::runtime_error` (SubstanceGroup.h:35-42)
  - plain `std::runtime_error` from checked iterators (ROMol.h:139, 144, 274)
  - `std::bad_any_cast` / `boost::bad_lexical_cast` from props
  - `boost::numeric::bad_numeric_cast` from the `rdcast` templates (unverified for each call site)
  - `FileParseException` from `MolOps::Hybridizations`' default and copy constructors (MolOps.h:1096-1102)

### What Python sees (nanobind)

| C++ | Python | Where |
|---|---|---|
| sanitize family | `rdchem.MolSanitizeException`, `AtomSanitizeException`, `AtomValenceException`, `AtomKekulizeException`, `KekulizeException`. All derive from `rdchem._BaseSanitException(ValueError)`. Translators attach `_type`, `_atom_idx`, `_atom_indices` | nb/rdchem.cpp:111-211 |
| `Invar::Invariant` | `RuntimeError(e.toUserString())`, so an out-of-range `GetAtomWithIdx` raises `RuntimeError`, not `IndexError`, after logging the failed expression **(observed)** | Code/RDBoost/nbWrap/RDBase.cpp:298-306 |
| `ValueErrorException` | `rdBase.ValueErrorException` (subclass of `ValueError`) | RDBase.cpp:288-289 |
| `IndexErrorException` | `rdBase.IndexErrorException` (subclass of `IndexError`) | RDBase.cpp:286-287 |
| `KeyErrorException` | `KeyError(key)` | RDBase.cpp:290-297 |
| `ConformerException` | `rdchem.ConformerException` (`ValueError`) | nb/Mol.cpp:220-221 |
| `std::runtime_error`, `SubstanceGroupException` | `RuntimeError` through nanobind's default translator (unverified) | — |
| `bad_any_cast` in `GetXProp` | `ValueError` | nb/props.hpp:177-181 |

Boost differences:

- `IndexError` and `ValueError` are raised as the builtins themselves (Code/RDBoost/Wrap.cpp:37-47).
- Invariant → `RuntimeError` (Wrap.cpp:63-72, under `#ifdef INVARIANT_EXCEPTION_METHOD`).
- The sanitize translator sets `cause` on the exception **class** object and raises with the message string (bp/rdchem.cpp:92-97).
- `ConformerException` → `ValueError("Bad Conformer Id")` (bp/rdchem.cpp:57-59; bp/Mol.cpp:278).

---

## Quirks and hazards

Tags: unmarked = read directly in source; **(inferred)** / **(unverified)** as defined at the top.

### A. `const` that mutates; raw pointers and ownership

- **H-A1** These `const` members mutate state or return mutable handles:
  - `ROMol::setName` (ROMol.cpp:591-593)
  - `ROMol::clearComputedProps` resets RingInfo (601-605)
  - `RDProps::setProp` / `clearProp` / `clearComputedProps`, through the `mutable` dict (RDProps.h:19, 80, 142, 155)
  - `Bond::getStereoAtoms() const` allocates through `const_cast` (Bond.h:355-360)
  - `ROMol::getRingInfo() const` returns a non-const `RingInfo*` (ROMol.h:773)
  - `Atom::getOwningMol() const`, `Bond::getOwningMol() const`, `Conformer::getOwningMol() const` and `SubstanceGroup::getOwningMol() const` return a non-const `ROMol&` (Atom.h:157; Bond.h:191; Conformer.h:91; SubstanceGroup.h:133)
  - `Bond::getBeginAtom/getEndAtom/getOtherAtom() const` return a non-const `Atom*` (Bond.h:270-282)
- **H-A2** Free functions that take `const ROMol&` and write to it:
  - `findSSSR`, `fastFindRings`, `findRingFamilies` fill RingInfo (MolOps.h:867-888)
  - `getAdjacencyMatrix`, `getDistanceMat`, `get3DDistanceMat` cache props (Matrices.cpp:230-232, 330, 429)
  - `Chirality::assignAtomCIPRanks` sets `_CIPRank` (Chirality.cpp:1308-1328)
- **H-A3** Python getters that perceive on demand: `Atom.IsInRing`, `IsInRingSize`, `Bond.IsInRing`, `IsInRingSize` call `findSSSR` whenever ring info is below SSSR. This replaces FAST ring info (nb/Atom.cpp:78-90; nb/Bond.cpp:66-79).
- **H-A4** Ownership conventions are mixed:
  - `addAtom` / `addBond` copy unless `takeOwnership` (ROMol.cpp:384-388, 418-422).
  - `addConformer` always takes ownership (ROMol.cpp:690-692).
  - `setMonomerInfo` takes ownership and deletes the previous value (Atom.cpp:633-636).
  - `Atom::copy` / `Bond::copy` return objects the caller owns.
  - Copy-returning functions differ: `addHs(const ROMol&)`, `removeHs`, `mergeQueryHs`, `adjustQueryProperties`, `renumberAtoms`, `hapticBondsToDative` and `combineMols` return raw `ROMol*` the caller owns; `copyMolSubset` returns `unique_ptr<RWMol>` (Subset.h:92); `getMolFrags` returns `shared_ptr` or `unique_ptr` (MolOps.h:108-136); v2 `MolFromSmiles` returns `unique_ptr`, v1 `SmilesToMol` a raw pointer (SmilesParse.h:69, 118).
- **H-A5** `createPartialBond` returns a heap `Bond*` that `finishPartialBond` neither inserts nor deletes (RWMol.cpp:632-653). The caller must delete it, and the header does not say so (RWMol.h:148-175). **(inferred)** leak.
- **H-A6** Public members expose internals:
  - `getAtomBookmarks()` / `getBondBookmarks()` return pointers to internal maps (ROMol.h:708, 731).
  - `RingInfo::dp_urfData` is a public data member (RingInfo.h:335-336).
  - `getSubstanceGroups(ROMol&)` returns a mutable reference to the vector (SubstanceGroup.h:274-275).
  - `Bond::setOwningMol` is public while `Atom::setOwningMol` is protected (Bond.h:196; Atom.h:410).
  - `setIdx`, `setBeginAtomIdx` and friends are public and do not touch the graph.

### B. Constructor, conversion and copy traps

- **H-B1** Python `RWMol(pkl, propertyFlags)` / `RWMol(bytes, propertyFlags)` call `new ReadWriteMol(std::string, unsigned)` (nb/Mol.cpp:1067-1081), but `ReadWriteMol` only has `(const ROMol&, bool quickCopy, int)` (123-127). **(observed)**: `RWMol(mol.ToBinary(AllProps), AllProps)` has no properties, while `Mol(...)` with the same arguments keeps them. The path:
  1. The pickle string converts implicitly through the non-explicit `ROMol(const std::string&)` (ROMol.h:538), which unpickles with AllProps.
  2. `propertyFlags` then lands in `quickCopy`, so any nonzero flag value drops props, conformers, sgroups and bookmarks.
  - Boost has the same constructor set (bp/Mol.cpp:182-186, 977).
- **H-B2** `ROMol(const std::string&)` and `RWMol(const ROMol&, bool=false, int=-1)` are non-explicit, so strings convert to molecules and `ROMol`s convert to full `RWMol` copies silently. Example: `insertSubstanceGroups(…, const RWMol &other, …)` (RWMol.cpp:85) is called with a `const ROMol&` (214), which deep-copies the whole source mol on every `insertMol`. **(inferred)**
- **H-B3** Copying a mol in the middle of a batch edit copies the pending-deletion bitsets, even into an `ROMol` (ROMol.cpp:121-130). `RWMol::clear()` leaves them in place (RWMol.h:204-209).
- **H-B4** A moved-from `ROMol` has `dp_ringInfo == nullptr` (ROMol.h:565, 576). Any later `clearComputedProps()` dereferences it (ROMol.cpp:604). **(inferred)**
- **H-B5** Python `Mol.__setstate__` destroys and placement-news the live C++ object (Wrap_nb.h:170-171). Any outstanding Atom, Bond, Conformer or RingInfo wrappers then dangle. **(inferred)**
- **H-B6** `ROMol` copy-assignment is deleted while `RWMol` copy-assignment exists (ROMol.h:607; RWMol.h:49).
- **H-B7** Name hiding: `RDProps::clear()` (RDProps.h:37) is public and non-virtual, and only resets props. `RWMol::clear()` (RWMol.h:204) wipes the whole molecule. Through an `ROMol&` that refers to an `RWMol`, `clear()` resolves to the props-only version.

### C. Special-member bugs and owner back-pointers

- **H-C1** `QueryAtom::operator=(QueryAtom&&)` calls `QueryAtom::operator=(std::move(other))`, which is itself (QueryAtom.h:79-86). `QueryBond` has the same shape (QueryBond.h:48-55). **(inferred)** infinite recursion.
- **H-C2** The `QueryBond` move constructor and move assignment use `std::move` on the raw `dp_query` pointer (QueryBond.h:45-47, 53, 117), which leaves the source still owning it. The destructor deletes it (QueryBond.cpp:23-26). **(inferred)** double delete.
- **H-C3** `Bond::operator=(const Bond&)` does not copy `d_stereo` and keeps the owner pointer, while the copy constructor copies `d_stereo` and drops the owner (Bond.cpp:25-44 vs 48-71).
- **H-C4** `Atom(Atom&&) = default` copies the owning raw `dp_monomerInfo` (Atom.h:134, 432), and `~Atom` deletes it (Atom.cpp:252). **(inferred)** double delete. Separately, `Atom::operator=` overwrites `dp_monomerInfo` without deleting the old one (Atom.cpp:215-219, 227-233). **(inferred)** leak.
- **H-C5** `EditableMol` (Python helper class) owns a raw `RWMol*` and has an implicitly generated copy constructor (nb/EditableMol.cpp:23-78). **(inferred)** double delete if copied in C++; not reachable from Python.
- **H-C6** Owner back-pointers: `Atom::dp_mol`, `Bond::dp_mol`, `Conformer::dp_mol`, `SubstanceGroup::dp_mol` (Atom.h:431; Bond.h:396; Conformer.h:156; SubstanceGroup.h:227). Many accessors PRECONDITION on them (`getOwningMol`, `getBeginAtom`, `hasValenceViolation`, `setStereoAtoms`, `getPerturbationOrder`).
- **H-C7** Copying preserves the owner pointer inconsistently:
  - Atom copy and assignment drop it (Atom.cpp:201).
  - Bond copy drops it; Bond assignment keeps it (Bond.cpp:27, 52).
  - Conformer copy and move keep it (Conformer.h:58-76).
  - So `Conformer(mol.getConformer())`, or Python `Chem.Conformer(conf)`, reports the original mol as its owner without being in it.

### D. Wrong answers, staleness and silent results

- **H-D1** Ring info goes stale silently:
  - Not reset by `addAtom`, `addBond` (both overloads), `insertMol`, `replaceAtom` or `replaceBond` (RWMol.cpp:144-325, 501-546; ROMol.cpp:377-432).
  - Per-index queries return 0 or empty for unknown indices (RingInfo.cpp:52-58, 129-135), and `isInitialized()` stays true.
  - **(observed)** Closing hexane into a ring with `RWMol.AddBond(0, 5)` leaves `NumRings()` at 0 and `IsInRing()` false.
- **H-D2** RingInfo fused-ring caches are never cleared, even by `reset()`, and are validated only by vector size (RingInfo.cpp:244-247, 211-223, 303-316). A re-perceived ring set of the same size reuses stale fused data. **(inferred)**
- **H-D3** `findRingFamilies` returns early when families already exist and ignores `includeDativeBonds` / `includeHydrogenBonds` (FindRings.cpp:1128-1135).
- **H-D4** `setAromaticity` and `Kekulize` use whatever ring info is initialized, including FAST rings, which the header says are only reliable for `numAtomRings > 0` (Aromaticity.cpp:1146-1150; Kekulize.cpp:664-668; MolOps.h:877-883).
- **H-D5** The valence cache is stale-prone:
  - Setters never invalidate it (Atom.h:148, 214, 218, 223).
  - `needsUpdatePropertyCache` only asks whether it was ever computed (Atom.cpp:744-747).
  - `removeBond`, `removeAtom` and `addBond(Bond*)` leave endpoint caches alone; only `addBond(i,j)` clears them (RWMol.cpp:531-532 vs 553-630; ROMol.cpp:403-432).
  - Atom copies carry the cache along (Atom.cpp:213-214).
- **H-D6** `getAdjacencyMatrix` has two problems:
  - It fills the array with `memset(res, emptyVal, …)` on `double`s (Matrices.cpp:308). **(inferred)** only `emptyVal=0` gives the intended value.
  - Its cache key ignores `emptyVal` and `bondsToUse` (297-304).
- **H-D7** Matrix ownership varies by overload:
  - `getAdjacencyMatrix` and the first `getDistanceMat` are cached in mol props; do not delete the result (MolOps.h:978-982, 1004-1008).
  - The `getDistanceMat(activeAtoms, bonds)` overload is not cached; the caller deletes (1031-1034).
  - `get3DDistanceMat`'s ownership depends on `propNamePrefix`. A `nullptr` prefix caches under `"_…"` (Matrices.cpp:395-399).
  - The nb wrapper deletes the result when the prefix is null or empty (nb/MolOps.cpp:538-540). The Python default is `""` (1761-1763), so that path is correct. Passing `None` would reach the double-free path only if nanobind maps `None` to `nullptr` (unverified).
  - The cache names `DistanceMatrix…` and `AdjacencyMatrix…` are not underscore-private (Matrices.cpp:175, 297), so they show up in `GetPropNames(includeComputed=True)`.
- **H-D8** CIP labels and the stereo flag outlive edits:
  - `_CIPCode` is stored non-computed by both legacy and new code (Chirality.cpp:1799; CIPLabeler/configs/AtropisomerBond.cpp:52), so it survives `clearComputedProps` and sanitization.
  - `_StereochemDone` is computed in `assignStereochemistry` (Chirality.cpp:2922) but set non-computed by MolInterchange/Parser.cpp:926, 990 and MolZip.cpp:654.
- **H-D9** `replaceBond` deletes the old `Bond*` without repointing StereoGroups, which store `Bond*` (RWMol.cpp:277-325; StereoGroup.h:48). `replaceAtom` does repoint atoms (262-274). **(inferred)** dangling pointer.
- **H-D10** `replaceAtom` ignores `updateLabel` (unnamed parameter, RWMol.cpp:235), though RWMol.h:94-95 and the Python docstrings say otherwise (nb/Mol.cpp:1104-1109; nb/EditableMol.cpp:117-121).
- **H-D11** Edits with surprising side effects:
  - `replaceBond` lowers endpoint explicit-H counts when bond order rises (RWMol.cpp:292-302).
  - `addBond(AROMATIC)` marks both atoms aromatic (512-522).
  - `setStereoGroups` merges ABS groups and moves them last (ROMol.cpp:434-464).
  - `removeHs(…, sanitize=true)` runs a full `sanitizeMol` (AddHs.cpp:1064-1066).
- **H-D12** `ROMol::getConformer` builds its error message with `mesg += id` where `id` is an `int`, which appends one `char` instead of the number (ROMol.cpp:659-661).
- **H-D13** `getUniqueAtomWithBookmark` / `getUniqueBondWithBookmark` do not check uniqueness, contrary to the header (ROMol.h:693-695, 716-718). They call `front()` on a possibly empty list (ROMol.cpp:236-240, 258-262). **(inferred)** undefined behavior.
- **H-D14** `getNumBonds(onlyHeavy=true)` returns every graph bond, including bonds to explicit H atoms; `false` adds per-atom H counts (ROMol.cpp:303-314). The Python docstring says "include only bonds to heavy atoms" (nb/Mol.cpp:493-498).
- **H-D15** Unowned atoms answer inconsistently:
  - `getDegree()` → 0 and `getValence()` → 0 (Atom.cpp:276-278, 320-322).
  - `getNumImplicitHs`, `getTotalNumHs` and `getTotalDegree` throw (300-309).
  - `hasValenceViolation`, `getPerturbationOrder` and Python `GetNeighbors` throw on the missing owner.
- **H-D16** Conformer corner cases:
  - `setAtomPos` grows the position vector silently, so Python `SetAtomPosition(aid >= n, …)` desynchronizes the conformer from the mol (Conformer.h:119-127).
  - `addConformer(…, assignId=false)` accepts duplicate ids (ROMol.cpp:678-694), and `getConformer(id)` returns the first match.
  - `removeConformer` with an unknown id is a silent no-op (669-676).
  - `df_is3D` defaults to true (Conformer.h:154).
  - `assignStereochemistryFrom3D` silently does nothing on non-3D conformers (Chirality.cpp:3770-3772).
- **H-D17** `Atom::Match` and `Bond::Match` are asymmetric by design (Atom.h:284-316; Atom.cpp:709-738). `Atom::getQuery()` returns `nullptr` while `Bond::getQuery()` throws a precondition (Atom.cpp:703; Bond.cpp:221-224).

### E. Preconditions that fail without prior perception

- **H-E1** Anything that reads the valence cache before `updatePropertyCache` hits a PRECONDITION, which Python sees as `RuntimeError`:
  - `getNumImplicitHs`, `getTotalNumHs`, `getTotalDegree`, `getTotalValence`, `getValence` (Atom.cpp:305-329)
  - `ROMol::getNumAtoms(false)`, `getNumBonds(false)` (ROMol.cpp:185-187, 309-311)
  - `setConjugation` and `setHybridization` (ConjugHybrid.cpp:163, 186-209)
- **H-E2** RingInfo queries PRECONDITION `df_init` (RingInfo.cpp:17-171). The inline `atomRings()`, `bondRings()`, `atomRingFamilies()` and `bondRingFamilies()` do not check, and return empty (RingInfo.h:122, 199, 304, 312). So Python `NumRings()` throws while `AtomRings()` returns `()` (nb/RingInfo.cpp:140, 145).
- **H-E3** `numRelevantCycles` and `atomRelevantCycles` also PRECONDITION that ring families exist (RingInfo.cpp:276-279, 322-324).
- **H-E4** `Bond::setStereo(STEREOCIS|STEREOTRANS)` PRECONDITIONs 2 stereo atoms (Bond.h:330-336). `setStereoAtoms` PRECONDITIONs that the atoms are neighbors (Bond.cpp:255-262).
- **H-E5** Every failed PRECONDITION logs to `rdErrorLog` before throwing, even when the caller catches it (Invariant.h:108-114).

### F. Defaults, out-parameters, boolean parameters, overloads

- **H-F1** Defaults that change answers:
  - `getNumBonds(onlyHeavy=1)`; `getConformer(id=-1)` returns the first conformer; `getTotalNumHs(includeNeighbors=false)`.
  - `addAtom(updateLabel=true)` sets a bookmark.
  - `getPropList(true, true)` in C++ vs Python `GetPropNames(False, False)`.
  - `assignStereochemistry(force=false)` no-ops once `_StereochemDone` is set.
  - `setAromaticity(DEFAULT)` = RDKit model; `symmetrizeSSSR(DEFAULT)` reads an env var.
  - `Kekulize(canonical=true)`: the header doc says `false` is the default (MolOps.h:729 vs 765; 777-780 vs 796).
- **H-F2** Out-parameters:
  - `sanitizeMol(…, unsigned &operationThatFailed, …)` (MolOps.h:573-575)
  - `getMolFrags` ×3 (72-112)
  - `getMolFragsWithQuery(…, map&)` (177-180)
  - `findSSSR(…, res&)` / `(…, res*)` (867-875)
  - `symmetrizeSSSR(…, res&)` (928-951)
  - `assignAtomCIPRanks(…, ranks&)` (Chirality.h:70-71)
  - Python `GetMolFrags(frags=[], fragsMolAtomMapping=[])` appends into caller-supplied lists (nb/MolOps.cpp:623-638).
- **H-F3** Boolean-parameter clusters:
  - `ROMol(other, bool, int)`
  - `replaceAtom(idx, atom, bool, bool)`, `replaceBond(idx, bond, bool, bool)`
  - `addHs(mol, bool, bool, ptr, bool)`; `removeHs(mol, bool, bool, bool)` (deprecated)
  - `getMolFrags(mol, bool, ptr, ptr, bool)`
  - `getDistanceMat(mol, bool, bool, bool, char*)`
  - `assignStereochemistry(mol, bool, bool, bool)`
  - `Kekulize(mol, bool, bool, unsigned)`
  - `symmetrizeSSSR(mol, res, enum, bool, bool, bool)`
- **H-F4** Overload pairs where the argument's type silently picks between in-place and copy:
  - `addHs`, `removeHs`, `removeAllHs`, `mergeQueryHs`, `adjustQueryProperties`, `hapticBondsToDative` and `dativeBondsToHaptic` mutate when given an `RWMol&` and return a new `ROMol*` when given a `const ROMol&` (MolOps.h:204-245, 292-349, 376-383, 489-494, 1296-1315).
  - `symmetrizeSSSR` has `bool` vs `SymmetrizeSSSRAlgorithm` overloads (934-959).
  - `findSSSR` has `res&` vs `res* = nullptr` (867-875).
  - `getNumAtoms()` vs `getNumAtoms(bool)` (ROMol.h:617-620).
- **H-F5** Return conventions differ between siblings:
  - `addAtom` returns the new index; `addBond` returns the bond count (RWMol.cpp:232, 545). The same holds in Python `AddAtom` / `AddBond`.
  - `KekulizeIfPossible` and `setAromaticity` return values are dropped in Python (nb/MolOps.cpp:445-449, 461-464).

### G. Python: in-place vs new, naming, defaults, lifetimes

- **H-G1** `Chem.Mol` is nominally read-only, yet module functions mutate it in place:
  - `SanitizeMol`, `Kekulize`, `KekulizeIfPossible`, `SetAromaticity`, `SetConjugation`, `SetHybridization`, `AssignRadicals`, `Cleanup`, `CleanupChirality`, `CleanupAtropisomers`, `CleanupOrganometallics`, `ExpandAttachmentPoints`, `CollapseAttachmentPoints`, all through `static_cast<RWMol&>(romol)` (nb/MolOps.cpp:417-498, 1097-1105). **(inferred)** formally undefined behavior when the dynamic type is `ROMol`; `RWMol` declares no data members (RWMol.h:33-221).
  - Also in place: `AssignStereochemistry`, the `Assign*` functions, `RemoveStereochemistry`, `FastFindRings`, `GetSSSR`, `GetSymmSSSR`.
  - These return new mols instead: `AddHs`, `RemoveHs`, `RemoveAllHs`, `MergeQueryHs`, `RenumberAtoms`, `CombineMols`, `DeleteSubstructs`, `AdjustQueryProperties`, `HapticBondsToDative`, `DativeBondsToHaptic`.
  - `removeStereochemistry` also `static_cast`s to `RWMol&` to call an `ROMol` member (Chirality.cpp:3863).
- **H-G2** Python defaults or names that differ from C++:
  - `Kekulize(clearAromaticFlags=False)` vs C++ `markAtomsBonds=true` (nb/MolOps.cpp:1824; MolOps.h:764).
  - `SetAtomMapNum(strict=False)` vs C++ `strict=true` (nb/Atom.cpp:551-553; Atom.h:382).
  - `FindPotentialStereo(cleanIt=False)`: C++ has no default (nb/ChiralityOps.cpp:23-30; Chirality.h:147-148).
  - `GetSymmSSSR` reorders arguments relative to `symmetrizeSSSR` (nb/MolOps.cpp:1306-1309; MolOps.h:928-932).
  - `GetPropsAsDict` defaults differ per class (§6).
  - `SetProp` on Atom and Bond has no `computed` argument (nb/Atom.cpp:59-63; nb/Bond.cpp:55-58).
  - `ClearComputedProps` exists only on Mol and Conformer.
  - `HasProp` returns `int` (nb/props.hpp:329-333).
  - `EditableMol.ReplaceBond` has no `keepSGroups` (nb/EditableMol.cpp:52-56, 122-125).
- **H-G3** Naming:
  - Python uses CamelCase with `Get`/`Set` prefixes, but `molzip` and `molzipFragments` are lowercase (nb/MolOps.cpp:3002-3038).
  - Boolean getters come in three shapes: `GetIsAromatic`, `HasQuery`, `IsInRing`.
  - `Is3D` / `Set3D`.
  - Two functions differ only in a capital letter: `DetectBondStereoChemistry(mol, conformer)` and `DetectBondStereochemistry(mol, confId=-1)`, bound to different C++ functions (nb/MolOps.cpp:1209, 1216).
  - Python ↔ C++ renames: `Cleanup` ↔ `cleanUp`; `CleanupOrganometallics` ↔ `cleanUpOrganometallics`; `AssignAtomChiralTagsFromStructure` ↔ `assignChiralTypesFrom3D`; `AssignAtomChiralTagsFromMolParity` ↔ `assignChiralTypesFromMolParity`; `GetSSSR` ↔ `findSSSR`; `GetSymmSSSR` ↔ `symmetrizeSSSR`; `ComputeAtomCIPRanks` ↔ `assignAtomCIPRanks`; `GetDistanceMatrix` ↔ `getDistanceMat`; `RemoveAllConformers` ↔ `clearConformers`; `GetAtomPosition` ↔ `getAtomPos`; Python `Mol` ↔ C++ `ROMol`.
- **H-G4** Wrong keyword name: `SetNoImplicit`'s argument is written `"what"_a "Sets a marker …"` (nb/Atom.cpp:236-240). **(inferred)** string-literal concatenation makes the keyword name the whole docstring, and the method gets no doc.
- **H-G5** Doc and binding nits:
  - `GetImplicitValence`'s runtime deprecation warning points at a `getExplicit=` keyword that does not exist (nb/Atom.cpp:147); `GetExplicitValence`'s warning correctly says `which=` (143). `GetValence` takes `which` (226). The two docstrings spell it `GetValence(Chem.ValenceType,EXPLICIT)` and `getValence(Chem.ValenceType,IMPLICIT)` (222, 225).
  - The C++ deprecation text says `getValence(true)` for an enum parameter (Atom.h:199-202).
  - `SetAtomAlias` and `SetAtomValue` name their second parameter `rlabel`; `GetAtomValue` / `SetAtomValue` are documented as "alias" (nb/Atom.cpp:610-619).
  - `Bond.GetEndAtomIdx` is documented as "first atom" (nb/Bond.cpp:159-160).
  - `GetRingInfo`: "Returns the number of molecule's RingInfo object" (nb/Mol.cpp:1021-1025).
  - The `Mol` docstring advertises bookmark methods that are not bound (nb/Mol.cpp:197-199).
  - `RingInfo::numRings` doc says it returns ring families under URF; the code returns `d_atomRings.size()` (RingInfo.h:188-189; RingInfo.cpp:170-174).
  - `sanitizeMol`'s documented step list is incomplete (§8).
- **H-G6** Python sequence holders capture `[vstart, vend)` at creation. Checked iterators only compare counts recorded at `begin()` (nb/seqholders.hpp:31-32, 53; ROMol.h:151-158, 237-243). **(inferred)** If the mol is edited between `GetAtoms()` / `GetNeighbors()` and iteration or indexing, `__len__` is stale and iteration can read past the graph or through invalidated adjacency iterators, undetected.
- **H-G7** `reference_internal` keeps the owner alive, not the element. These Python wrappers dangle after the C++ object is deleted **(inferred)**:
  - `Atom` after `RemoveAtom`, `ReplaceAtom` or `CommitBatchEdit` (RWMol.cpp:496, 249, 868)
  - `Bond` after `RemoveBond` or `ReplaceBond` (628, 310)
  - `Conformer` after `RemoveConformer` or `RemoveAllConformers` (ROMol.cpp:669-676; ROMol.h:750)
  - `MonomerInfo` after `SetMonomerInfo` (Atom.cpp:633-636)
  - `SubstanceGroup` objects returned by `CreateMolSubstanceGroup`, `AddMolSubstanceGroup` or `GetMolSubstanceGroupWithIdx` after any later `push_back` reallocates the vector (nb/SubstanceGroup.cpp:33-69, 352-370; SubstanceGroup.cpp:363-372)
- **H-G8** `RingInfo.AddRing(None, None)`: both lists become `nullptr`, the size check passes (0 == 0), and `atomIds->begin()` dereferences null (nb/RingInfo.cpp:98-113; Wrap_nb.h:101-111). **(inferred)**
- **H-G9** `GetProp(key, autoConvert=True)` on an existing `AnyTag` property hits `break`, falls through, and raises `KeyError` (nb/props.hpp:267-269, 310-313). **(inferred)** `GetPropsAsDict` silently skips such properties (124-126).
- **H-G10** `_BaseSanitException.GetType()` returns `self._msg`, the message, while the translators store the type in `_type`, which no method reads (nb/rdchem.cpp:117-118, 145, 160, 175, 189). **(observed)** for an `AtomValenceException`. In Boost, `e.cause.GetType()` returns the C++ type name and `cause` is a class attribute that every raise overwrites (bp/rdchem.cpp:92-97).
- **H-G11** Python enum `StereoDescriptor` lacks `Bond_AtropCW` and `Bond_AtropCCW` (nb/Chirality.cpp:41-48 vs Chirality.h:103-111). What happens when one of these values is returned to Python is **(unverified)**.
- **H-G12** `AROMATICITY_CUSTOM` is exposed (nb/MolOps.cpp:1922), but `SetAromaticity` cannot take a function. **(inferred)** It always fails the PRECONDITION (Aromaticity.cpp:1167-1171).
- **H-G13** Leftover debug output:
  - `QueryAtom.ExpandQuery` prints to `std::cerr` on every call (nb/Atom.cpp:41, 46-47).
  - The unused `ReadOnlySeq` prints in its constructor and destructor (nb/seqholders.hpp:102-118).
  - `ReadWriteMol` has an unused `std::shared_ptr<RWMol> dp_mol` (nb/Mol.cpp:190-191).
- **H-G14** `RWMol.SetStereoGroups` validates that atoms belong to the mol but not bonds (nb/Mol.cpp:151-166).
- **H-G15** Names and keywords present in only one binding:
  - Boost only: `GetAttachmentPointLabelNumber`, `IsMarkedAttachmentPoint` (bp/MolOps.cpp:3418, 3432), `Mol.AddConformersFromTrajectory` (bp/Mol.cpp:484).
  - Boost `GetNumAtoms` also takes a deprecated `onlyHeavy` keyword (bp/Mol.cpp:433-435).
  - nanobind only: `__getstate__` / `__setstate__`, native `Mol.GetAtoms` / `GetBonds`, `StereoInfo.NOATOM`.
  - Boost `Atom.GetNeighbors` / `GetBonds` return tuples, which slice; nanobind returns holders indexed by `int` only.
- **H-G16** Process-global monkeypatching from optional modules: `IPythonConsole` replaces `Mol.GetSubstructMatch` (rdkit/Chem/Draw/IPythonConsole.py:484-485); `PandasTools` replaces `Mol.__ge__` (rdkit/Chem/PandasTools.py:656-657).

### H. Global and build-time switches

- **H-H1** Three runtime switches live in the process **environment**. Setters call `setenv`; getters call `getenv` on every read:
  - `RDK_USE_LEGACY_STEREO_PERCEPTION` (default true)
  - `RDK_ENABLE_NONTETRAHEDRAL_STEREO` (default true)
  - `RDK_USE_LEGACY_RING_FINDING` (default false)
  - Sources: Chirality.h:35-47; Chirality.cpp:836-858; MolOps.h:810-818; FindRings.cpp:936-946.
  - Any value other than `"0"`, including `"false"`, reads as true (utils.h:96-106).
  - Read by `assignStereochemistry`, `findPotentialStereo` / `FindStereo.cpp`, and `symmetrizeSSSR`, which `sanitizeMol` calls.
- **H-H2** `extern bool Chirality::useLegacyStereoPerception` is declared twice (Chirality.h:51-57) and is never defined under `Code/`. A grep of `Code/` found no definition.
- **H-H3** `MolPickler` default pickle properties form a process-global `NoProps` (MolPickler.cpp:78, 270-283), so `pickle` / `ToBinary()` drop props unless the global is changed.
- **H-H4** Build-time switches:
  - `INVARIANT_EXCEPTION_METHOD` / `ASSERT_METHOD` / `SILENT_METHOD` change whether PRECONDITIONs throw, assert or vanish (Invariant.h:34-46, 161-180). `URANGE_CHECK` takes 3 arguments in ASSERT mode and 2 otherwise (141 vs 168).
  - `USE_NEW_STEREOCHEMISTRY` switches the CIP ranking algorithm (Chirality.cpp:1315-1321).
  - On WIN32, `ROMol`'s `initMol`, `destroy`, `addAtom`, `addBond` and `initFromOther` are protected rather than private (ROMol.h:1035-1038).
  - `RDK_BUILD_THREADSAFE_SSS` gates the pickler mutex (MolPickler.cpp:271-280).

### I. Data-model issues

- **H-I1** The class docs say properties "have no impact upon chemistry", which is the justification for `const` setters (ROMol.h:97-99; Atom.h:53-55; Bond.h:43-45). Yet chemistry lives in props:
  - atom map number (Atom.h:382-398)
  - non-tetrahedral `_chiralPermutation` (Atom.cpp:837-861)
  - dummy symbol `dummyLabel` (Atom.cpp:269-271)
  - haptic endpoints `_MolFileBondEndPts` (RWMol.cpp:418-456; MolOps.h:1363-1365)
  - `_CIPCode` and `_StereochemDone`
- **H-I2** Atom setters narrow silently: `setAtomicNum(int)` into `uint8`, `setFormalCharge(int)` into `int8`, `setNumExplicitHs` / `setNumRadicalElectrons(unsigned)` into `uint8`, `setIsotope(unsigned)` into `uint16` (Atom.h:148, 209, 214, 223; Atom.cpp:642; storage Atom.h:416-427). **(inferred)** out-of-range values truncate.
- **H-I3** Four reference models coexist (§12): pointers for StereoGroups and bookmarks; indices for SubstanceGroups, stereo atoms and RingInfo; strings for haptic endpoints. Each edit path patches each model by hand (RWMol.cpp:253-274, 375-481, 571-616).
- **H-I4** Hot paths are linear: `getBondWithIdx` is O(idx) (ROMol.cpp:324-329); `CXXBondIterator::size()` is O(n) (ROMol.h:331-339); Python bond holders are O(n) per `len` or index; `Dict` lookup is linear (Dict.h:166-185).
- **H-I5** `MolOps::Hybridizations` has default and copy constructors that throw `FileParseException`, a misuse guard with an unrelated exception type (MolOps.h:1094-1120).

---

## Size

### C++ public member functions

Counts are declarations, overloads included; constructors, destructor and assignment are excluded. Distinct names are in parentheses.

| Class | Own | Plus inherited `RDProps` | Location |
|---|---|---|---|
| ROMol | 108 (68) | 12 (10); plus 2 conditional `save/load` | ROMol.h:411-1003 |
| RWMol | 22 (16) on top of ROMol → 130 total | — | RWMol.h:57-216 |
| Atom | 59 (55) | + RDProps | Atom.h:143-406 |
| Bond | 45 (41) | + RDProps | Bond.h:160-389 |
| Conformer | 17 (12) | + RDProps | Conformer.h:82-144 |
| RingInfo | 39 (39) | public data member `dp_urfData` | RingInfo.h:49-323 |

`MolOps.h`:

- 83 function declarations (62 distinct names) in `namespace MolOps`.
- Another 4 in `MolOps::details`.
- 4 types: `AddHsParameters`, `RemoveHsParameters`, `AdjustQueryParameters`, `Hybridizations`.
- 4 enums: `SanitizeFlags`, `AromaticityModel`, `AdjustQueryWhichFlags`, `SymmetrizeSSSRAlgorithm`.

### Python (nanobind)

Counts are `.def` registrations; overloads count separately. Distinct names in parentheses.

| Class | Methods | Ctors | Location |
|---|---|---|---|
| `Mol` | 64 (48: 44 + 4 dunders) | 6 | nb/Mol.cpp:426-1029 |
| `RWMol` | +18 (18: 12 + 6 dunders); the effective API is 62 distinct names | 6 | nb/Mol.cpp:1054-1133 |
| `Atom` | 73 (68), plus static `NOATOM` | 3 | nb/Atom.cpp:172-554 |
| `QueryAtom` | +2 | 0 | nb/Atom.cpp:595-601 |
| `Bond` | 49 (44) | 0 | nb/Bond.cpp:104-336 |
| `QueryBond` | +2 | 0 | nb/Bond.cpp:382-388 |
| `Conformer` | 32 (26) | 3 | nb/Conformer.cpp:106-302 |
| `RingInfo` | 26 (26) | 0 | nb/RingInfo.cpp:123-159 |
| `EditableMol` | 10 (10) | 1 | nb/EditableMol.cpp:103-135 |
| `StereoGroup` | 6 (6) | 0 | nb/StereoGroup.cpp:88-103 |
| `SubstanceGroup` | 44 (39) | 0 | nb/SubstanceGroup.cpp:155-345 |
| module `rdmolops` | 109 registrations (100 distinct names) | — | nb/MolOps.cpp + nb/ChiralityOps.cpp |
