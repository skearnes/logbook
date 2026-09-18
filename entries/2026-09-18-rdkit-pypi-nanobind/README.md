# 2026-09-18 — RDKit on PyPI: what the nanobind wrappers change

- **Date:** 2026-09-18
- **Author:** Steven Kearnes
- **Acknowledgments:** Prepared with [Claude Code](https://claude.com/claude-code) (Claude Opus 5)
- **Status:** draft
- **Tags:** rdkit, nanobind, packaging, pypi, conda-forge, wheels, boost
- **License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

## Summary

RDKit's documented reason for treating conda-forge, not PyPI, as its Python
distribution is a 2019 blog post whose central claim is about Boost.Python. The
nanobind wrappers remove that claim: nanobind extensions share bound types
through the interpreter, not through a shared binding library, so the constraint
the post calls unavoidable no longer holds.

What remains after that is smaller and different in kind:

- RDKit's own 65 shared libraries still have to be vendored into a wheel, so the
  wheel is still not a C++ distribution.
- The wheel is still built and owned outside `rdkit/rdkit`, by a different
  author, and RDKit's CI does not build or test it.

The first is a real limit worth stating rather than fixing. The second is an
ownership gap, and closing it is the cheaper and larger win. Adopting nanobind's
stable ABI would collapse the wheel matrix from 22 files to 4.

Follows the
[2026-09-12 next-generation entry](../2026-09-12-rdkit-next-generation/README.md),
which left "how PyPI wheels ship `rdkit.v3`" open. Measured against the
`nanobind-narrow-objects` branch of `rdkit/rdkit` at `50a0242`, nanobind 2.15.0,
and `rdkit` 2026.3.6 on PyPI.

## Evidence

### The 2019 argument

Greg Landrum,
[Why the RDKit isn't available on PyPi](http://rdkit.blogspot.com/2019/11/why-rdkit-isnt-available-on-pypi.html),
2019-11-15, makes three claims:

1. **Not a pure Python package.** "The core problem is that the RDKit is not a
   pure python package; it's a mix of python and some compiled extension modules
   (shared libraries)."
2. **Boost.Python must be dynamically linked.** "But Boost::Python needs to be
   dynamically linked if you want to share types across extension modules."
3. **The alternatives cost too much.** Static linking "would increase the size of
   the packages pretty dramatically (the statically linked Windows conda packages
   are about twice as big as those for Linux or the Mac)," and one monolithic
   extension module would be large and slow to import.

The post recommends conda: "conda was designed to solve exactly the
cross-platform dependency management problem described above. And it does it
pretty well." It also names its own escape hatch: "I would be really happy if it
turns out that I'm wrong about this and it is, in fact, possible to get things
set up so that the RDKit can be made 'pip installable.'"

Claim 2 is the load-bearing one. Claim 1 describes a difficulty that vendoring
tools already handle for many projects; claim 3 is an argument against the
alternative that claim 2 forces. Without claim 2, claim 3 has nothing to rule
out.

### Why claim 2 held for Boost.Python

Boost.Python keeps its to- and from-Python converter registry as global state in
`libboost_python`. Statically linking it gives each extension module a private
registry, so a molecule created in one module is an unregistered type in the
next. The shared library is what makes the registry one registry.

That has a second cost visible in the build. The Boost.Python component is named
per Python minor version
([`CMakeLists.txt:464`](https://github.com/rdkit/rdkit/blob/master/CMakeLists.txt#L464)):

```cmake
find_package(Boost ${RDK_BOOST_VERSION} COMPONENTS "python${Python_VERSION_MAJOR}${Python_VERSION_MINOR}" "numpy${Python_VERSION_MAJOR}${Python_VERSION_MINOR}" REQUIRED CONFIG)
```

So every new CPython release needs a new Boost binary before RDKit can be built
against it, on conda-forge as much as on PyPI.

### How nanobind shares types

nanobind 2.15.0 does not use a shared binding library for this. `nb_module_exec`
looks up a capsule in the interpreter state dictionary, keyed by ABI tag and
domain, and creates it if absent (`src/nb_internals.cpp:533-556,643-649`):

```cpp
PyObject *dict = PyInterpreterState_GetDict(PyInterpreterState_Get());
PyObject *key = PyUnicode_FromFormat("__nb_internals_%s_%s__",
                                     abi_tag(), name ? name : "");
```

The header states the consequence directly (`src/nb_internals.h:370-374`):
`nb_internals` "storing information related to function/type bindings and
instances. Separate nanobind extensions within the same NB_DOMAIN furthermore
share `nb_internals` to communicate with each other."

Two independently built nanobind extensions with the same ABI tag and domain
therefore share one type registry whether or not they share a library. Claim 2
does not survive this.

RDKit passes `NB_SHARED`
([`RDKitUtils.cmake:228`](https://github.com/rdkit/rdkit/blob/master/Code/cmake/Modules/RDKitUtils.cmake#L228)),
which does install a `libnanobind`. That is a size choice — one copy of a small
runtime instead of 55 — not a correctness requirement, and it is RDKit's own
build artifact rather than an ABI-pinned dependency shared with the rest of the
distribution. `NB_STATIC` would work too.

### What the stable ABI would be worth

nanobind can build limited-API extensions. `STABLE_ABI` compiles the runtime
with `Py_LIMITED_API=0x030C0000` (`cmake/nanobind-config.cmake:440`) and is
silently dropped below Python 3.12, on non-CPython interpreters, without
`Python::SABIModule`, and on free-threaded builds (`:389-404`). Boost.Python
cannot use the limited API at all.

RDKit requires Python 3.9
([`CMakeLists.txt:372`](https://github.com/rdkit/rdkit/blob/master/CMakeLists.txt#L372))
and does not pass `STABLE_ABI`.

The current matrix, from the PyPI JSON API for `rdkit` 2026.3.6:

| | value |
| --- | --- |
| wheels in the release | 22 |
| CPython tags | 6 (`cp310`–`cp315`) |
| platform tags | 4 (`macosx_11_0_arm64`, `manylinux_2_28_x86_64`, `manylinux_2_28_aarch64`, `win_amd64`) |
| wheel size | 25–38 MB |
| total per release | 717 MB |

One abi3 wheel per platform makes that 4 files, and makes a new CPython release
stop being a release-blocking event.

### How much of the port is left

- 51 `nbWrap/` directories beside 51 `Wrap/` directories.
- 55 `NB_MODULE` against 69 `BOOST_PYTHON_MODULE` in non-test sources. Almost all
  of the difference is demo code under `Code/Demos/boost/` — `moduleA`,
  `moduleB`, `moduleC`, `moda`, `modb`, `overloads`, `python_objs`, `linalg` —
  which is not shipped.
- Exactly two shipped modules have no nanobind build:
  - `cQuantize`, declared with `rdkit_python_extension`
    (`Code/ML/Data/CMakeLists.txt:1`), which is a no-op when the Boost wrappers
    are off (`RDKitUtils.cmake:194,217`), so a nanobind build simply lacks it;
  - `rdTrajectory`, commented out at
    `Code/GraphMol/nbWrap/CMakeLists.txt:69-72`.

Wrapper behavior parity is tracked separately in
[rdkit#9601](https://github.com/rdkit/rdkit/issues/9601).

### Compiled Boost that remains

`serialization` and `iostreams` are still required
([`CMakeLists.txt:584-595`](https://github.com/rdkit/rdkit/blob/master/CMakeLists.txt#L584-L595)).
21 files include `boost/archive` or `boost/serialization`; 9 include
`boost/iostreams`.

Both are used only from C++, so a single static copy inside RDKit's own
libraries has no cross-module registry to split — the failure mode that forced
Boost.Python to be shared does not apply. Removing both would leave no Boost
shared library in the wheel at all. `iostreams` is the smaller job and mostly
gzip stream wrapping; `serialization` defines the pickle and binary formats, so
it carries a compatibility constraint rather than only an implementation choice.

### Who builds the wheel

- PyPI `rdkit` 2026.3.6 records author `Christopher Kuenneth` and homepage
  `https://github.com/kuelumbus/rdkit-pypi`.
- `rdkit/rdkit` contains no `pyproject.toml` and no `setup.py`; the root
  `setup.cfg` configures yapf.
- Neither `azure-pipelines.yml` nor `.github/workflows/` builds a wheel.
- The same repository does publish MinimalLib to npm on each release, with npm
  Trusted Publishing (`.github/workflows/minimallib-publish.yml`).

So the project already runs a publish pipeline from its own CI for a non-conda
artifact — for JavaScript, not for Python.

## Assessment

The 2019 post is a correct argument about Boost.Python that has been read since
as an argument about wheels. Once the nanobind wrappers are the shipped build,
its central claim no longer describes RDKit, and the remaining objections are
weaker than the ones it makes:

| 2019 claim | status under nanobind |
| --- | --- |
| Not a pure Python package | Holds. 65 `rdkit_library()` targets still get vendored and RPATH-rewritten. |
| Boost.Python must be dynamically linked | Retired. Type sharing goes through the interpreter, not a library. |
| Static linking or one big module costs too much | Moot. It was the cost of the alternative claim 2 forced. |

What is left is not a reason the wheel cannot work — it demonstrably does, at 22
files a release — but two narrower facts: a wheel is not a C++ distribution, and
this wheel is not built by the project that it ships.

## Options

Ranked, with a recommendation on each.

### 1. Finish the port, then move the floor to 3.12 and build `STABLE_ABI`

**Do this.** Two shipped modules remain, and the stable ABI is the only change
here that alters the economics rather than the tidiness: 4 wheels instead of 22,
and no scramble on each CPython release, for conda-forge as much as for PyPI.
The costs are real and should be measured first — abi3 forgoes type-slot
features nanobind otherwise uses, and dropping Python 3.9–3.11 is a user-facing
decision, not a build one. Note the stable ABI is silently disabled rather than
diagnosed when its preconditions fail (`nanobind-config.cmake:389-404`), so a
build that thinks it is producing abi3 wheels needs to assert on the tag.

### 2. Move the wheel build into `rdkit/rdkit`

**Do this, and it is independent of 1.** A `pyproject.toml` with scikit-build-core
and a `cibuildwheel` job in the existing CI turns an artifact the project does
not control into one it tests. This is what answers the maintenance half of the
2019 position — bug reports against the wheel currently land on a tracker whose
maintainers did not build it. The MinimalLib workflow is the precedent for the
release plumbing, and the work is mostly moving what `kuelumbus/rdkit-pypi`
already does, which argues for doing it with its author rather than around them.

### 3. Drop the remaining compiled Boost

**Worth doing, lower priority.** `iostreams` first: 9 files, and zlib covers what
it is used for. `serialization` is 21 files and touches the on-disk formats, so
it needs its own plan. The payoff is a wheel with no vendored Boost, which
shrinks the surface but does not unblock anything.

### 4. Ship a C++-usable wheel

**Don't.** Headers, a CMake config, and sonames left unmangled by `auditwheel`
and `delocate` could be added the way pyarrow and torch do it, but the
PostgreSQL cartridge, the SWIG Java wrappers, and Qt cannot ship this way
regardless. conda-forge should stay the complete
distribution, and the wheel should be documented as Python bindings only. The
2019 post is right about this part.

## Open questions

- Whether abi3's cost to call overhead is acceptable for RDKit's hot paths, which
  needs the same kind of measurement the July prototype used for de-inlining.
- Whether dropping Python 3.9–3.11 to get abi3 is acceptable, and on what
  schedule relative to those versions' upstream end of life.
- Whether `cQuantize` and `rdTrajectory` are ported or dropped; both have been
  left aside long enough that the answer may be "dropped."
- Whether conda-forge's `rdkit` output moves to the nanobind wrappers at the same
  release the wheels do, or before.
- What a nanobind build does to wheel size, given claim 3's numbers were about
  Boost.
- How `rdkit.v3` ships on PyPI, still open from the
  [2026-09-12 entry](../2026-09-12-rdkit-next-generation/README.md#open-questions).

## Next steps

- Build the nanobind wrappers with `STABLE_ABI` against Python 3.12 and check
  that the produced tag really is `abi3`, then benchmark against the non-abi3
  build.
- Port or retire `cQuantize` and `rdTrajectory`, so the nanobind build is a
  superset of what the Boost build ships.
- Ask on [rdkit#9601](https://github.com/rdkit/rdkit/issues/9601) when the
  nanobind wrappers become the default, since both options above wait on it.
- Sketch the `pyproject.toml` and `cibuildwheel` matrix against what
  `kuelumbus/rdkit-pypi` already builds, to size option 2 before proposing it.
- Measure what removing `boost::iostreams` would take, starting from the 9 files
  that include it.
