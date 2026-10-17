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

# Check that the mol.h sketch and the calls in example.cpp compile as C++20.
#
# Usage: check.sh COMPILER...
#   COMPILER  a C++20 compiler such as clang++ or g++-16

set -eu

here=$(cd "$(dirname "$0")" && pwd)
for compiler in "$@"; do
  flags=(-std=c++20 -fsyntax-only -Wall -Wextra -Wpedantic -Werror)
  echo '#include "mol.h"' | "$compiler" "${flags[@]}" -I"$here" -x c++ -
  "$compiler" "${flags[@]}" "$here/example.cpp"
  echo "$compiler: ok"
done
