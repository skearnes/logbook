// Copyright 2026 Steven Kearnes
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

static int helper() { return 1; }

namespace {
struct Hidden {
  int value = 2;
};
}  // namespace

static constexpr int limit = 3;

inline int uses_static_function() { return helper(); }
inline int uses_anonymous_type() { return Hidden{}.value; }
inline int uses_constant_value() { return limit; }
inline const int *uses_constant_address() { return &limit; }
