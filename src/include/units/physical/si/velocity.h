// The MIT License (MIT)
//
// Copyright (c) 2018 Mateusz Pusz
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <units/physical/dimensions.h>
#include <units/physical/si/length.h>
#include <units/physical/si/time.h>
#include <units/quantity.h>

namespace units::si {

struct metre_per_second : unit<metre_per_second> {};
struct dim_velocity : physical::dim_velocity<dim_velocity, metre_per_second, dim_length, dim_time> {};

struct kilometre_per_hour : deduced_unit<kilometre_per_hour, dim_velocity, kilometre, hour> {};

template<Unit U, Scalar Rep = double>
using velocity = quantity<dim_velocity, U, Rep>;

inline namespace literals {

// mps
constexpr auto operator"" mps(unsigned long long l) { return velocity<metre_per_second, std::int64_t>(l); }
constexpr auto operator"" mps(long double l) { return velocity<metre_per_second, long double>(l); }

// kmph
constexpr auto operator"" kmph(unsigned long long l) { return velocity<kilometre_per_hour, std::int64_t>(l); }
constexpr auto operator"" kmph(long double l) { return velocity<kilometre_per_hour, long double>(l); }

}  // namespace literals

}  // namespace units::si
