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

#include <mp-units/quantity_spec.h>
#include <mp-units/systems/isq/base_quantities.h>
#include <mp-units/systems/isq/mechanics.h>
#include <mp-units/systems/isq/space_and_time.h>

namespace mp_units::isq {

// TODO Add all the remaining ISQ definitions
QUANTITY_SPEC(activity, inverse(duration));
QUANTITY_SPEC(absorbed_dose, energy / mass);
QUANTITY_SPEC(ionizing_radiation_quality_factor, dimensionless);
QUANTITY_SPEC(dose_equivalent, absorbed_dose* ionizing_radiation_quality_factor);

}  // namespace mp_units::isq
