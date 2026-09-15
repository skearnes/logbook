# 2026-09-12 — Next-generation RDKit: scope, shape, and mechanics

- **Date:** 2026-09-12
- **Author:** Steven Kearnes
- **Acknowledgments:** Prepared with [Claude Code](https://claude.com/claude-code) (Claude Opus 5)
- **Status:** draft (working)
- **Tags:** rdkit, next-generation, api-design, cpp20, cpp-modules, nanobind,
  packaging, conda-forge, proposal
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

## Summary

Working notes toward a proposal for a next-generation RDKit API that ships
beside the existing one in the `rdkit` package. The proposal settles scope,
shape, and mechanics; the new API itself and the chemistry under it come later.
Decisions so far:

- The new API owns the types whose rules it changes and re-exports others by
  recorded decision; it inherits nothing implicitly.
- Its types wrap existing code, and implementations are replaced one type at a
  time.
- The new API lives in its own directory behind a CI-enforced boundary.
- It is written as C++20 headers that can gain module interfaces later.
- In Python it is `rdkit.v3`, bound with nanobind and built only with RDKit's
  nanobind wrappers. Molecules convert between the APIs by copying.
- conda-forge ships it as new outputs of the existing recipe.

Which types below the molecule it re-exports, and whether that code becomes a
layer both APIs share, is open.

Follows [2026-07-18 — Reimplementing RDKit in Rust](../2026-07-18-rdkit-rust-reimplementation/README.md)
and the [2026-07-15 maintenance plans](../2026-07-15-rdkit-maintenance-plans/README.md).
Measured against RDKit `master` at `23378a7` (2026-09-12) and conda-forge
`rdkit-feedstock` at `6918f4e` (rdkit 2026.03.6).

## Constraints

From feedback on the July entries:

- Existing code keeps working; API changes are otherwise welcome.
- C++20 or newer, no Rust. C++ modules are welcome once support matures, but
  not required.
- One package: installing `rdkit` provides both generations, the way pydantic 2
  ships `pydantic.v1`.
- Little duplicated code.
- Scope is shape and mechanics, not the new API or its routines (the hydrogen
  model, for example).
- Splitting RDKit into several projects is an option.

## Evidence

### Where the code sits

| Layer | Contents | Lines |
| --- | --- | --- |
| Below the molecule type | Libraries that don't link `GraphMol`, even indirectly: `RDGeneral`, `DataStructs`, `Geometry`, `Numerics/Alignment`, `Numerics/EigenSolvers`, `ML/InfoTheory`, `ML/Cluster/Murtagh`, `SimDivPickers`, `Catalogs`, `ChemicalFeatures`, `RDStreams` | ~20k |
| Molecule core | Files directly in `GraphMol/`, plus `SmilesParse`, `Substruct`, `CIPLabeler` | ~62k |
| Built on the molecule type | Everything else | ~220k |
| Bindings | `Wrap/` (Boost.Python) ~41k, `nbWrap/` (nanobind) ~38k, `JavaWrappers` and `MinimalLib` ~10k | ~88k |

Counts are lines in `.cpp`, `.h`, `.hpp`, `.c`, `.cc`, `.yy`, and `.ll` files
under `Code/`. Non-binding rows skip files named `*test*`, `*catch*`, or
`*bench*` and the `Demos`, `Fuzz`, `Bench`, and `RDBoost` directories. Link
dependencies come from each `rdkit_library()` call's `LINK_LIBRARIES`.

- Most code is written against `ROMol` and `RWMol`; only ~20k lines sit below
  the molecule type.
- The boundary is blurry. The generic force-field core and `Numerics/Optimizer`
  reach `GraphMol` only through `Trajectory/Snapshot.h`, which needs only
  `Geometry`; the MMFF and UFF terms include molecule headers directly.

### RDKit's `v2` C++ API

- Since 2024.03.1, nine parser and supplier headers declare a `v2` namespace
  beside an `inline namespace v1` whose functions call v2
  (`FileParsers.h:131-139`, `SmilesParse.h:132`). The release notes keep v1
  supported "for the forseeable future".
- Nothing is duplicated because both versions share `RWMol`. v2 changes
  ownership and parameter passing but keeps v1's defaults. A namespace could
  change behavior too, just not a shared type's invariants
  ([Decision 1](#1-the-new-api-owns-the-types-whose-rules-it-changes)).

### What existing code can do to a molecule

- `const` doesn't mean unchanged. `RDProps` keeps properties in a
  `mutable Dict`, "a quirk of history" (`RDProps.h:19-20`), so `setProp`,
  `clearProp`, and `clearComputedProps` are `const`.
- `getRingInfo() const` returns a non-const `RingInfo *` (`ROMol.h:774`), which
  `findSSSR`, `fastFindRings`, and `findRingFamilies` fill from a
  `const ROMol &` (`MolOps.h:867-886`).
- `MolToSmiles` takes a `const ROMol &` and stores `_smilesAtomOutputOrder` and
  `_smilesBondOutputOrder` on it (`SmilesWrite.cpp:747-750`).
- A stale cache fails when read: `getNumImplicitHs()` checks at runtime that
  `calcImplicitValence()` has run (`Atom.cpp:305-307`).
- In Python, `Chem.SanitizeMol`, `Chem.Kekulize`, `Chem.SetAromaticity`,
  `Chem.AssignStereochemistry`, and `AllChem.EmbedMolecule` change their
  argument in place.
- Behavior changes have shipped as process-wide switches. 2022.09 added a
  global flag for new stereo perception
  ([rdkit#5309](https://github.com/rdkit/rdkit/pull/5309)), and 2023.03 removed
  the per-call `SmilesParserParams.useLegacyStereo` in its favor.
- The flag is an environment variable, `RDK_USE_LEGACY_STEREO_PERCEPTION`: the
  setter calls `setenv`, and every query reads it back
  (`Chirality.cpp:848-858`). Legacy perception is still the default
  (`Chirality.h:37-39`), and non-tetrahedral stereo is switched the same way.

### RDKit's nanobind port

- [rdkit#9030](https://github.com/rdkit/rdkit/pull/9030), merged 2026-08-11,
  put an `nbWrap/` beside each of the 59 `Wrap/` directories, `rdBase`
  included (`Code/RDBoost/nbWrap/RDBase.cpp`). Both builds run the same Python
  tests.
- A build uses one: Boost.Python is on by default, nanobind is off, and enabling
  both is a CMake error.
- Modules use `nanobind_add_module(... NB_SHARED ...)` without `NB_DOMAIN`
  (`RDKitUtils.cmake:228`), so extension modules see each other's bound types.
  nanobind is pinned to 2.x ([rdkit#9535](https://github.com/rdkit/rdkit/pull/9535)).
- With shared types, each C++ type has one Python class. Binding a type again
  warns and returns the class registered first before setting any module
  attribute, and the second binding's `def` calls add methods to that class
  (nanobind 2.15.0 `src/nb_type.cpp:1323-1339,1665`,
  `include/nanobind/nb_class.h:666-688`). `ROMol` is bound as `Mol`
  (`GraphMol/nbWrap/Mol.cpp:426`).
- The discussion that started the port ([rdkit#9031](https://github.com/rdkit/rdkit/discussions/9031))
  reported two to five times lower wrapper overhead than Boost.Python in early
  `Point2D` benchmarks.
- conda-forge builds the Boost.Python wrappers (`libboost-python-devel` in
  `host`).

### How TensorFlow and pydantic ship two generations

TensorFlow is the closer analog, a C++ core under a Python API:

- TensorFlow 2 gave `tf` to the new API and kept 1.x as `tf.compat.v1`. Release
  1.14 added `tf.compat.v2` so libraries could publish code that runs on both.
- The `tf_export` decorator records each symbol's v1 and v2 names, and both
  namespaces are generated from those records, so unchanged symbols are written
  once.
- `tf_upgrade_v2` rewrites imports to `tf.compat.v1`. Removed modules such as
  `tf.contrib` had no path.
- Mixing is limited. The generations differ in process-wide runtime behaviors,
  eager execution can be switched only once per program, and some
  `tf.compat.v1` symbols don't work under TF2 behaviors.

pydantic 2.13.4 bundles a frozen 13,166-line copy of pydantic 1.10.26 as
`pydantic.v1`, and 1.10.17 added the same namespace to the 1.x line as a bridge.
Its engine ships separately as `pydantic-core`, and its migration guide says V1
and V2 models can't be mixed.

What carries over is one package and one implementation behind both surfaces.
What doesn't:

- giving existing top-level names to the new API, which the constraint on
  existing code rules out;
- a bridge release, which RDKit doesn't need because existing names never move;
- process-wide switches, which Decision 1 avoids;
- a frozen copy, which for RDKit's C++ would duplicate everything and still have
  to track compilers, Boost, Python, and NumPy.

### C++ modules today

- CMake supports them from 3.28 with the Ninja (1.11+) or Visual Studio
  generators, on GCC 14+, Clang 16+, or MSVC 14.34+
  ([cmake-cxxmodules](https://cmake.org/cmake/help/latest/manual/cmake-cxxmodules.7.html)).
  AppleClang isn't supported: Apple's clang 21 ships no `clang-scan-deps` and
  rejects `export module`.
  RDKit requires CMake 3.18 and builds as C++20.
- CMake doesn't support header units. `import std;` is experimental, and Clang
  documents mixing it with `#include`d standard headers as poorly supported
  ([Clang modules documentation](https://clang.llvm.org/docs/StandardCPlusPlusModules.html)).
- Compiled module interfaces are compiler- and flag-specific, so libraries ship
  interface sources.
- Existing headers can be wrapped. Included in a module's global module fragment
  (`module;` before `export module`), their declarations stay in the global
  module, and `export using` publishes chosen names. Included after
  `export module`, they would belong to the new module and stop matching
  `librdkit`. The `std` modules of libc++, libstdc++, and MSVC's STL all start
  with `module;` (`std.cppm.in:14`, `std.cc.in:24`, `std.ixx:5`, fetched
  2026-09-15).
- Macros don't cross `import`. Of 520 headers, 53 use `PRECONDITION`, 20
  `CHECK_INVARIANT`, and 7 `BOOST_LOG`.
- A module interface can't expose an internal-linkage entity, such as a
  namespace-scope `static` function, from inline code or declarations. Reading
  an internal constant's value is exempt, and `static` members take their
  class's linkage ([basic.link](https://eel.is/c++draft/basic.link)). GCC
  accepts exposures of internal entities from the global module fragment with a
  warning, `-Wexpose-global-module-tu-local`
  ([GCC options](https://gcc.gnu.org/onlinedocs/gcc/C_002b_002b-Dialect-Options.html)).
- A probe built here with Homebrew Clang 23.1.1 and GCC 16.2.0, against
  conda-forge's `librdkit-dev` 2026.03.6 and Boost 1.90
  ([run.sh](assets/module_probe/run.sh)):
  - A module wrapping `ROMol.h` and `SmilesParse.h` in a global module fragment
    compiles on both. An importer parses a SMILES through it without naming
    `RWMol`, and on Clang links against `librdkit` and runs. Naming
    `RDKit::ROMol` in the importer fails on both.
  - GCC warns 30 times that the wrapped headers expose internal constants, 27
    times from Boost and 3 from RDKit (`RDTypeTag::AnyTag`, `zero_tolerance`).
    Clang reports none.
  - A function attached to a named module gets a module-specific symbol
    (`_ZW8attached6answerv`), and a wrapped one keeps its ordinary symbol.
  - Wrapped, GCC warns on an inline function that uses an internal function or
    a constant's address and rejects one that uses an anonymous-namespace type;
    Clang accepts all three. Attached, GCC rejects all three and Clang warns
    only on the function. Reading a constant's value passes everywhere.
- Requiring CMake 3.28 at the root enables policy CMP0155, which scans every
  C++20 source and excludes scanned sources from unity builds, the main lever in the July [compile-time plan](../2026-07-15-rdkit-maintenance-plans/assets/COMPILE_TIME_PLAN.md).
  Policies are set per directory and recorded on each target, and
  `CXX_SCAN_FOR_MODULES` overrides them per target or source
  ([cmake_policy](https://cmake.org/cmake/help/latest/command/cmake_policy.html),
  [UNITY_BUILD](https://cmake.org/cmake/help/latest/prop_tgt/UNITY_BUILD.html)).

### The conda-forge recipe

- One CMake build feeds five outputs. `librdkit` (`Unspecified`, `base`,
  `data`, `runtime`), `librdkit-dev` (`dev`), and `rdkit` (`extras`, `python`,
  the `stubs` target, and a `pip install` for `.dist-info`) install by CMake
  component. `rdkit-dev` is an alias, and `rdkit-postgresql` runs the
  cartridge's install script.
- `librdkit` is pinned exactly through `run_exports`, and every dependent pins
  it exactly.
- The generators are `make` and `NMake Makefiles JOM`, and neither scans
  modules.
- The compilers are `gxx` 15, conda-forge's `clangxx` 21 (not AppleClang), and
  `vs2022`. `linux-aarch64`, `linux-ppc64le`, and `osx-arm64` are
  cross-compiled.
- `rdkit` builds the Boost.Python wrappers and generates stubs with
  `pybind11-stubgen`. `pyproject.patch` packages `rdkit`, `rdkit.*`, and
  `rdkit-stubs`.

## Decisions

Each decision ends with a strict alternative: how it would work, what it would
enable, and why it isn't the choice now. Strictness that a future capability
such as C++ modules needs is worth adopting early when it is cheap.

### 1. The new API owns the types whose rules it changes

- A type whose rules, behavior, or representation the new API changes is its
  own, holding the existing object privately as a wrapper rather than a `using`
  alias. The molecule, atoms, and bonds are the expected case, along with any
  type whose public members reach them, such as `Conformer`
  (`getOwningMol()`, `Conformer.h:91`) and `StereoGroup` (`Atom *` members,
  `StereoGroup.h:47`).
- A type it uses unchanged, such as a point or a bit vector, can be
  re-exported. Each re-export, like anything else it inherits (a function, a
  default), is an explicit, recorded decision.
- Existing objects enter owned types only through explicit conversion.

Why owned types:

- **A shared type can't gain rules.** Existing code keeps working, so a shared
  molecule must keep allowing every state existing code builds and relies on. A
  rule the new API wants, such as the July hydrogen model's rule that hydrogens
  are always graph nodes, needs a type existing code can't reach.
- **Nothing guards a shared molecule.** Namespaces change free functions, not
  member functions, and `const` functions change molecules too: writing a
  SMILES adds properties, and ring finding fills ring info
  ([evidence](#what-existing-code-can-do-to-a-molecule)). An `RWMol` can pass
  through existing functions and return in a state the new API forbids, so the
  new API would recheck its rules on every entry. An owned type checks once, at
  conversion, and after that only its own methods change it.
- **Global switches don't mix.** RDKit shipped new stereo perception behind an
  environment variable in 2022.09, removing the per-call parser option, and it
  is still off by default. A switch applies to every library in the process,
  the same limit on mixing TensorFlow's generations. An owned type can carry a
  new behavior as its own, so both behaviors run in one process and meet at
  conversion, though code it wraps still reads existing switches
  ([Open questions](#open-questions)).
- **Later changes stay possible.** Decision 2 replaces implementations behind
  types, but through an alias `ROMol`'s members become API, down to
  `getTopology()` returning a Boost Graph Library graph (`ROMol.h:877`). If
  C++20 can make "forgot to sanitize" a compile error, that state lives in the
  type, and a shared `RWMol` has no room for it.

Owning costs Decision 4's non-inlined call per access and Decision 6's copy at
conversion.

A re-export avoids both, and old and new code share its objects. It is recorded
with what it gives up:

- its header joins the public headers (Decision 3);
- argument-dependent lookup finds its namespace's free functions for new-API
  callers;
- its implementation can't be replaced (Decision 2) without breaking new-API
  code;
- in Python it keeps its existing class and names unless `rdkit.v3` wraps it in
  a Python class (Decision 6).

Strict alternative: own every type.

- **How.** Wrap every type the API exposes, points and bit vectors included, so
  public headers include no existing header and need no approved exceptions.
- **Enables.** Public headers and a future module interface never pull in
  existing headers, their macros, or Boost. Every type stays replaceable,
  Python gets PEP 8 names without wrapper classes, and argument-dependent
  lookup never reaches existing namespaces.
- **Not now.** Each value type would need a wrapper and a copy at every
  boundary, against the constraint on duplicated code. Modules don't require
  it, since interface files can wrap existing headers (Decision 4).
- **Adopt when.** A re-exported type needs a new representation, or a module
  interface must not reach existing headers.

### 2. Implemented over existing code, replaced type by type

- A new type starts with existing code behind it. Its implementation is replaced
  when the existing representation can't hold what the API needs (an unknown
  coordinate, say) or measurement shows it is too slow.
- Duplication is limited to what each replacement forces.
- If differential tests show a replacement reproduces existing behavior exactly,
  the existing API may become an adapter over it, as v1 is over v2. This is
  optional.

Strict alternative: wrapped routines read no process-wide switch.

- **How.** Before the new API wraps an existing routine, each switch the
  routine reads becomes a parameter that defaults to the switch for existing
  callers. Stereo perception, for example, reads its environment variable at
  `Chirality.cpp:2658,2919`.
- **Enables.** New-API behavior stops depending on environment variables and
  other libraries' settings, which closes the open question on switches, and
  differential tests can pin both behaviors in one process.
- **Not now.** It changes existing code routine by routine, and which routines
  the new API wraps isn't decided.
- **Adopt when.** The first wrapped routine reads a switch.

### 3. A new directory with an enforced boundary

- The new API gets a top-level directory beside `Code/`, as its own CMake
  project built from the root behind an option that is on by default from the
  start.
- CI enforces the boundary. Public headers include only new-API headers, the
  standard library, and existing headers explicitly approved under Decision 1;
  implementation files may include any existing header.
- A directory comes first for three reasons:
  - changes to the existing code the new API wraps need both test suites in one
    CI run;
  - the new API links the existing implementation;
  - one tarball means one feedstock, with no cross-feedstock ABI pins while the
    API is changing.
- It becomes a separate project if any of these happens:
  - it needs its own release schedule;
  - a different group develops it;
  - RDKit splits first and the code the new API wraps becomes its own project;
  - it no longer links `librdkit`'s molecule code.

Strict alternative: build against an installed RDKit.

- **How.** A CI job installs RDKit, then configures the new directory on its
  own with `find_package(rdkit)`, so it sees only installed headers and
  exported targets.
- **Enables.** The directory can become a separate project whenever a trigger
  above fires, and any use of build-tree headers or targets fails in CI.
- **Not now.** It adds a second configure and build to CI while no split is
  planned.
- **Adopt when.** A split trigger becomes likely.

### 4. Headers now, modules later

CI checks from the start keep headers ready for modules:

- The public API works without macros. Checks and logging are functions built
  on `std::source_location`.
- Every header compiles on its own.
- Public headers declare no anonymous namespaces. GCC rejects an inline
  function that uses an anonymous-namespace type even from a global module
  fragment ([evidence](#c-modules-today)).
- Components have no dependency cycles, since modules can't import in a cycle,
  and each maps to one future module.
- The include boundary from Decision 3 holds.

Other internal linkage isn't linted; the second strict alternative below says
when it would be.

Consequences:

- Without modules, a wrapped type can't be usable yet unnameable. Wrappers
  therefore use pimpl and pay a non-inlined call per access. The July prototype
  measured de-inlining `ROMol`'s hot accessors at under 1%, with mixed sign ([BGL_DECOUPLING_FINDINGS.md](../2026-07-15-rdkit-maintenance-plans/assets/BGL_DECOUPLING_FINDINGS.md)),
  and `getAtomWithIdx()` is already out of line. With modules, wrappers can
  inline.
- Ninja and CMake 3.28 aren't needed before adopting modules.
- Adding modules means interface files that wrap the same headers, shipped
  beside them. The new directory sets CMP0155 in its own scope while the root
  keeps its CMake 3.18 policies, so scanning covers only the new targets and
  existing unity builds survive.
- Adopt modules when all of these hold:
  - conda-forge's `vs2022` toolset is 14.34+ (`gxx` 15 and `clangxx` 21
    already qualify);
  - the compiler behind the macOS PyPI wheels can scan modules;
  - a wrapper module over `ROMol.h` builds on GCC, Clang, and MSVC (it builds
    with GCC 16 and Clang 23; the feedstock's compilers and MSVC are
    unchecked);
  - the feedstock builds with Ninja.
- `import std;` waits until CMake's support for it isn't experimental.

Strict alternative: modules first.

- **How.** Write the new API as named module interfaces from the start,
  attached to `rdkit.v3.*` modules, with existing headers only in implementation
  units or global module fragments. Every build uses CMake 3.28+ with Ninja or
  Visual Studio, and the internal-linkage lint below applies, since GCC rejects
  those exposures inside a named module.
- **Enables.** Wrapped types stay reachable but unnameable, so wrappers inline
  without pimpl, and neither macros nor existing names reach importers.
- **Not now.** conda-forge builds with `make` and `NMake Makefiles JOM`, and
  Apple's clang 21 rejects `export module` (checked with `-std=c++20`). An
  entity attached to a named module can't also be declared in a header
  ([basic.link](https://eel.is/c++draft/basic.link)) and gets a different
  symbol, so SWIG, MinimalLib, and the cartridge would need a separate
  interface.
- **Adopt when.** The adoption conditions above hold and every consumer of the
  new API can import modules.

Strict alternative: lint internal linkage.

- **How.** Public headers declare no namespace-scope `static` function or
  variable, and their namespace-scope constants are `inline constexpr`.
- **Enables.** The same headers compile attached to a named module, as modules
  first needs: there GCC rejects an inline function that uses an internal
  function or a constant's address ([evidence](#c-modules-today)).
- **Not now.** Wrapped headers don't need it. GCC only warns on those
  exposures, Clang accepts them, and MSVC is unchecked.
- **Adopt when.** With modules first, or when a compiler rejects a wrapped
  exposure. New headers make it nearly free to adopt early.

### 5. Packaging

- The new API installs through its own CMake components into two new conda-forge
  outputs, a runtime library and a `-dev` package, each pinned exactly. `rdkit`
  depends on the runtime, so `conda install rdkit` provides both generations.
  `librdkit` and `rdkit-postgresql` are unchanged.
- `nanobind` joins the host requirements. `libboost-python-devel` leaves once
  the recipe builds the nanobind wrappers.
- As a separate project, the new API would get its own feedstock, as
  `pydantic-core` has.

Strict alternative: a stable C++ ABI.

- **How.** Version the new API's shared library by ABI, check each release with
  an ABI checker such as libabigail, and pin dependents to a compatible range
  instead of an exact build.
- **Enables.** C++ packages built against the new API survive RDKit releases
  without rebuilds, and the new API could later release on its own schedule.
- **Not now.** The API is still changing, and exact pins keep one feedstock
  simple (Decision 3).
- **Adopt when.** The API is declared stable, or a split trigger makes separate
  releases likely.

### 6. `rdkit.v3`

- The Python subpackage is `rdkit.v3`, with headers under `rdkit/v3/` and any
  future C++ modules named `rdkit.v3.*`.
- At the start the new API is reached only by importing `rdkit.v3` explicitly.
  Exposing it from the top level comes later ([Deferred](#deferred)).
- Why `v3`:
  - it stays accurate once it is the recommended API;
  - it follows the parsers' `RDKit::v1` and `RDKit::v2` namespaces (2024.03);
  - it matches `pydantic.v1` and `tf.compat.v1`;
  - it clashes with no subpackage, even ignoring case, and can't be confused
    with RDKit's year.month releases.
- The docs need a sentence on why Python has no `rdkit.v2`.
- **Proposed C++ namespace: `rdkit::v3`.**
  - Unqualified names resolve through enclosing namespaces, and no language
    feature prevents that, even in out-of-line definitions. Under `rdkit::v3`
    the enclosing `rdkit` declares nothing, so the compiler rejects existing
    names used without `RDKit::`.
  - Later generations become `rdkit::v4` and so on, and sibling namespaces
    aren't searched.
  - The only lint needed bans `using namespace RDKit` in the new directory.
  - The cost is two top-level namespaces that differ only by case. The tree has
    no lowercase `rdkit` namespace.
- **Alternative: `RDKit::v3` with a linter.** It keeps one top-level namespace,
  but unqualified names fall back to `RDKit` and silently reach existing
  functions. Blocking that takes a custom Clang check (clang-query or
  clang-tidy) covering names in expressions, types, and templates, maintained
  across Clang versions.
- In either namespace, argument-dependent lookup finds `RDKit` functions for
  `RDKit` arguments, which implementation code that wraps existing objects
  relies on.
- Python names follow PEP 8: snake_case functions, methods, and modules, and
  CapWords classes. The existing API keeps its names and modules.
- The bindings use nanobind only. They are their own set of wrappers but share
  the existing wrappers' nanobind domain, so conversion can pass existing
  objects. They bind no C++ type the existing wrappers bind: the second binding
  gets the class registered first, sets no module attribute, and adds its
  methods to that class, so the result depends on import order
  ([evidence](#rdkits-nanobind-port)).
- Molecules convert by copying in both directions, because a shared molecule
  would let existing code break the new type's rules. The conversions and
  re-exports are the only places `rdkit.v3` names existing types. The existing
  API never depends on `rdkit.v3`.

Strict alternative: a separate nanobind domain.

- **How.** Build `rdkit.v3` with `NB_DOMAIN`, which gives it its own nanobind
  library and type registry (nanobind 2.15.0
  `cmake/nanobind-config.cmake:425-437`), and convert molecules through a copy
  that passes no bound object, such as RDKit's binary pickle.
- **Enables.** `rdkit.v3` can bind any C++ type, re-exported ones included,
  under PEP 8 names with no import-order effects, and can change nanobind
  versions on its own.
- **Not now.** Conversion couldn't take an existing `Mol` directly, and a
  re-exported type would get a second Python class that existing functions
  reject.
- **Adopt when.** `rdkit.v3` ships apart from the existing wrappers or needs its
  own nanobind version.

## Top-level imports

Importing `rdkit` or any subpackage, `rdkit.v3` included, runs
`rdkit/__init__.py`. That file:

- imports the compiled `rdBase`, whose comment says exceptions leak memory
  otherwise, and takes `__version__` from it; loading `rdBase` also sets up C++
  logging (`Code/RDBoost/nbWrap/RDBase.cpp:285`);
- under Jupyter or Colab, imports `rdkit.Chem.Draw.IPythonConsole` (and with it
  `rdkit.Chem`, `rdchem`, `rdChemReactions`, and `rdMolDraw2D`) and sends RDKit
  logs to stderr;
- configures the `rdkit` logger;
- patches iteration on `rdBase` vector types.

### Keeping Boost.Python out of `rdkit.v3`

`rdkit.v3` is built only with the nanobind wrappers:

- In that build `rdBase` is a nanobind module, so `rdkit/__init__.py` loads no
  Boost.Python as it stands.
- Conversion goes through nanobind. An `rdkit.v3` function can take and return
  an existing `Mol`, because the extension modules share bound types.
- `rdkit.v3` ships once the nanobind wrappers are the build RDKit ships.

### Deferred

- **`from rdkit import Mol`.** Re-exporting only adds names, so this stays
  possible, with two limits:
  - a new top-level name can't reuse a subpackage's name, since importing
    `rdkit.Geometry` would overwrite a class `Geometry`;
  - a new top-level module can't be a lowercase version of an existing one,
    since `rdkit/geometry` and `rdkit/Geometry` are one directory on default
    macOS and Windows filesystems.
- **Import cost.** Until `rdkit/__init__.py` changes, importing `rdkit.v3` also
  loads `rdBase` and, in notebooks, the existing drawing stack.

Revisiting the top-level imports must preserve:

- `rdkit.__version__`;
- `rdkit.rdBase` without a separate import;
- `rdBase` loading before other extension modules, which also sets up C++
  logging;
- the `rdkit` logger setup;
- notebook rendering of existing molecules;
- the vector iteration patch.

## Revised positions

- **Rust is out.** Its July appeal was typestate. Whether C++20 can make "forgot
  to sanitize" a compile error is an API question left open.
- **No clean-room rewrite of the core.** Replacing implementations type by type,
  behind owned types, keeps that option without paying for it up front. The [hydrogen model spec](../2026-07-18-rdkit-rust-reimplementation/assets/HYDROGEN_MODEL_SPEC.md) remains input to later API work.
- **Shared types, not namespaces, limit behavior changes** (Decision 1).

## Open questions

- Which types below the molecule the new API re-exports, and whether that code
  becomes a layer both APIs share. Answering it needs a survey of those ~20k
  lines and of the force-field code that reaches `GraphMol` only through
  `Snapshot`.
- When the nanobind build becomes the default upstream and on conda-forge, since
  `rdkit.v3` waits on it.
- How `rdkit.v3` keeps process-wide switches read by the code it wraps from
  changing its behavior. Stereo perception reads an environment variable on
  every call, and the SMILES parser's per-call option was removed in its favor
  ([evidence](#what-existing-code-can-do-to-a-molecule)).
- Stubs and pickling for `rdkit.v3`: nanobind's stub generator or the recipe's
  `pybind11-stubgen`.
- Support for the existing API: fixes only or features too, and for how long.
- Whether and when Java (SWIG), JavaScript (MinimalLib), and the PostgreSQL
  cartridge get the new API.
- How PyPI wheels ship `rdkit.v3`, and what that does to wheel size.
- Which RDKit splits beyond the new directory are worth making, and when.
- Which accessors need inlining before modules arrive.
- How behavior that differs between the two APIs is documented and tested.

## Next steps

- Survey the code below the molecule type: what it contains, which types the
  new API would re-export, which macros would need function forms, and what
  moving `Snapshot` would free.
- Ask RDKit's maintainers when the nanobind wrappers become the default build.
- In a nanobind build, prototype copying an existing `Mol` into a stand-in
  `rdkit.v3` type and back.
- Specify the CI checks: include boundary, macro-free public API,
  self-contained headers, no anonymous namespaces in public headers, no
  component cycles, and the namespace rule from Decision 6.
- Repeat the [module probe](assets/module_probe/run.sh), which passes with
  Homebrew GCC 16 and Clang 23, with conda-forge's `gxx` 15, `clangxx` 21, and
  `vs2022`. This checks one of Decision 4's adoption conditions.
