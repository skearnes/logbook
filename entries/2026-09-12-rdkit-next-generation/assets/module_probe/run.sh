#!/usr/bin/env bash
# Copyright 2026 Steven Kearnes
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Probe how Clang and GCC build C++20 modules that wrap RDKit headers.
#
# Usage: run.sh PREFIX CLANG GCC
#   PREFIX  conda prefix with librdkit-dev and libboost-headers installed
#   CLANG   a clang++ that builds named modules
#   GCC     a g++ that builds named modules with -fmodules
#
# Written for macOS: Clang compiles against the SDK from xcrun, and only Clang
# links and runs the importer, because conda-forge builds librdkit with libc++.

set -u

prefix=$1
clang=$2
gcc=$3
here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
mkdir "$work/clang" "$work/gcc"

sysroot=(-isysroot "$(xcrun --show-sdk-path)")
includes=(-I"$prefix/include/rdkit" -I"$prefix/include")

strip_paths() { sed -E 's#^[^ :]*/##'; }

echo "== compilers"
"$clang" --version | head -n 1
"$gcc" --version | head -n 1

echo
echo "== symbols: attached to a named module, or wrapped after module;"
for unit in attached wrapped; do
  cd "$work/clang"
  "$clang" -std=c++20 "${sysroot[@]}" -x c++-module "$here/$unit.cppm" \
    --precompile -o "$unit.pcm"
  "$clang" -std=c++20 "${sysroot[@]}" -c "$unit.pcm" -o "$unit.o"
  echo "clang $unit: $(nm "$unit.o" | grep -i answer)"
  cd "$work/gcc"
  "$gcc" -std=c++20 -fmodules -x c++ "$here/$unit.cppm" -c -o "$unit.o"
  echo "gcc $unit: $(nm "$unit.o" | grep -i answer)"
done

echo
echo "== internal-linkage exposures, wrapped (fragment) or attached (purview)"
for unit in fragment purview; do
  cat > "$work/use_$unit.cpp" <<EOF
import $unit;
int main() {
  return uses_static_function() + uses_anonymous_type() +
         uses_constant_value() + *uses_constant_address();
}
EOF
  cd "$work/clang"
  echo "clang $unit:"
  "$clang" -std=c++20 "${sysroot[@]}" -x c++-module "$here/$unit.cppm" \
    --precompile -o "$unit.pcm" 2>&1 | grep -E 'error|warning' | strip_paths
  "$clang" -std=c++20 "${sysroot[@]}" -fmodule-file="$unit=$unit.pcm" \
    -c "$work/use_$unit.cpp" -o "use_$unit.o" 2>/dev/null &&
    echo "  importer compiled"
  cd "$work/gcc"
  echo "gcc $unit:"
  "$gcc" -std=c++20 -fmodules -x c++ "$here/$unit.cppm" -c -o "$unit.o" \
    2>&1 | grep -E ': (error|warning):' | strip_paths
  "$gcc" -std=c++20 -fmodules -c "$work/use_$unit.cpp" -o "use_$unit.o" \
    2>/dev/null && echo "  importer compiled"
done

echo
echo "== a module wrapping ROMol.h and SmilesParse.h"
cd "$work/clang"
"$clang" -std=c++20 "${sysroot[@]}" "${includes[@]}" -x c++-module \
  "$here/wrapped_mol.cppm" --precompile -o rdkit.v3.wrapped_mol.pcm \
  2> clang_wrapped_mol.log
echo "clang module: $(grep -c ': error:' clang_wrapped_mol.log) errors," \
  "$(grep -c ': warning:' clang_wrapped_mol.log) warnings"
"$clang" -std=c++20 "${sysroot[@]}" -c rdkit.v3.wrapped_mol.pcm -o wrapped_mol.o
"$clang" -std=c++20 "${sysroot[@]}" "${includes[@]}" \
  -fmodule-file=rdkit.v3.wrapped_mol=rdkit.v3.wrapped_mol.pcm \
  -c "$here/import_mol.cpp" -o import_mol.o
"$clang" "${sysroot[@]}" import_mol.o wrapped_mol.o -nostdlib++ \
  -L"$prefix/lib" -Wl,-rpath,"$prefix/lib" -lc++ \
  -lRDKitSmilesParse -lRDKitGraphMol -lRDKitRDGeneral -o import_mol
echo "clang importer runs: $(./import_mol)"
"$clang" -std=c++20 "${sysroot[@]}" \
  -fmodule-file=rdkit.v3.wrapped_mol=rdkit.v3.wrapped_mol.pcm \
  -c "$here/names_existing.cpp" -o names_existing.o 2> clang_names.log
echo "clang importer naming RDKit::ROMol: $(grep -m 1 error clang_names.log | strip_paths)"

cd "$work/gcc"
"$gcc" -std=c++20 -fmodules "${includes[@]}" -x c++ "$here/wrapped_mol.cppm" \
  -c -o wrapped_mol.o 2> gcc_wrapped_mol.log
echo "gcc module: $(grep -c ': error:' gcc_wrapped_mol.log) errors," \
  "$(grep -c ': warning:' gcc_wrapped_mol.log) warnings"
grep -o "exposes TU-local entity '[^']*'" gcc_wrapped_mol.log | sort | uniq -c
"$gcc" -std=c++20 -fmodules "${includes[@]}" -c "$here/import_mol.cpp" \
  -o import_mol.o && echo "gcc importer compiled"
"$gcc" -std=c++20 -fmodules -c "$here/names_existing.cpp" \
  -o names_existing.o 2> gcc_names.log
echo "gcc importer naming RDKit::ROMol: $(grep -m 1 error gcc_names.log | strip_paths)"
