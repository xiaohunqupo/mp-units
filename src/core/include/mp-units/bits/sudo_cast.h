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

#include <mp-units/bits/hacks.h>
#include <mp-units/ext/type_traits.h>
#include <mp-units/framework/quantity_concepts.h>
#include <mp-units/framework/reference_concepts.h>
#include <mp-units/framework/unit.h>
#include <mp-units/framework/unit_magnitude.h>

namespace mp_units::detail {

template<typename... Ts>
constexpr bool has_common_type_v = requires { typename std::common_type_t<Ts...>; };

template<typename T, typename Other>
using maybe_common_type =
  std::conditional_t<has_common_type_v<T, Other>, std::common_type<T, Other>, std::type_identity<T>>::type;

/**
 * @brief Magnitude-only details about a unit conversion factor
 *
 * Factored out of `conversion_type_traits` and `conversion_value_traits` so that the
 * magnitude decomposition (which depends only on `M`, not on any representation type)
 * is computed once per unique magnitude and shared across all representation combinations.
 *
 * @note This is a low-level facility.
 *
 * @tparam M common magnitude between the two quantities
 */
template<UnitMagnitude auto M>
struct magnitude_traits {
  using c_mag_type = common_magnitude_type<M>;
  static constexpr UnitMagnitude auto num = numerator(M);
  static constexpr UnitMagnitude auto den = denominator(M);
  static constexpr UnitMagnitude auto irr = M * (den / num);
};

/**
 * @brief Type-related details about the conversion from one quantity to another
 *
 * This trait helps to determine what representations to use at which step in the conversion process,
 * in order to avoid overflow and underflow while not causing excessive computations.
 *
 * @note This is a low-level facility.
 *
 * @tparam M common magnitude between the two quantities
 * @tparam Rep1 first quantity representation type
 * @tparam Rep2 second quantity representation type
 */
template<UnitMagnitude auto M, typename Rep1, typename Rep2>
struct conversion_type_traits {
  using c_rep_type = maybe_common_type<Rep1, Rep2>;
  using c_mag_type = magnitude_traits<M>::c_mag_type;
  using multiplier_type = conditional<
    treat_as_floating_point<c_rep_type>,
    // ensure that the multiplier is also floating-point
    conditional<std::is_arithmetic_v<value_type_t<c_rep_type>>,
                // reuse user's type if possible
                std::common_type_t<c_mag_type, value_type_t<c_rep_type>>, std::common_type_t<c_mag_type, double>>,
    c_mag_type>;
  using c_type = maybe_common_type<c_rep_type, multiplier_type>;
};

/**
 * @brief Value-related details about the conversion from one quantity to another
 *
 * This trait provide ingredients to calculate the conversion factor that needs to be applied
 * to a number, in order to convert from one quantity to another.
 *
 * @note This is a low-level facility.
 *
 * @tparam M common magnitude between the two quantities
 * @tparam T common multiplier representation type
 */
template<UnitMagnitude auto M, typename T>
struct conversion_value_traits {
  using mag = magnitude_traits<M>;
  static constexpr T num_mult = get_value<T>(mag::num);
  static constexpr T den_mult = get_value<T>(mag::den);
  static constexpr T irr_mult = get_value<T>(mag::irr);
  static constexpr T ratio = num_mult / den_mult * irr_mult;
};


/**
 * @brief Single point of intentional narrowing/truncation with compiler diagnostics disabled
 *
 * Every `static_cast` that intentionally converts to a lower-precision type (e.g. `long double`
 * intermediate → `double` result) must go through this helper so that the diagnostic suppression
 * macros appear in exactly one place and are easy to audit.
 *
 * @tparam To target type
 * @tparam From source type (deduced)
 */
template<typename To, typename From>
[[nodiscard]] constexpr To silent_cast(From value) noexcept
{
  MP_UNITS_DIAGNOSTIC_PUSH
  MP_UNITS_DIAGNOSTIC_IGNORE_FLOAT_CONVERSION
  return static_cast<To>(value);
  MP_UNITS_DIAGNOSTIC_POP
}

/**
 * @brief Numerical scaling of a value between two units
 *
 * Contains all the scaling logic that depends only on the source/target unit and representation
 * types. By factoring this out of `sudo_cast` (which is parametrised on the full `Quantity` type),
 * the expensive instantiation is shared across all quantity types that happen to have the same
 * unit and representation — e.g. `quantity<m, double>`, `quantity<isq::length[m], double>`, and
 * `quantity<isq::radius[m], double>` all reuse the same `sudo_cast_value<m, double, km, double>`
 * instantiation.
 *
 * @note This is a low-level facility.
 *
 * @tparam FromUnit source unit
 * @tparam FromRep  source representation type
 * @tparam ToUnit   target unit
 * @tparam ToRep    target representation type
 */
template<Unit auto FromUnit, typename FromRep, Unit auto ToUnit, typename ToRep>
[[nodiscard]] constexpr ToRep sudo_cast_value(FromRep value)
{
  constexpr UnitMagnitude auto c_mag =
    mp_units::get_canonical_unit(FromUnit).mag / mp_units::get_canonical_unit(ToUnit).mag;
  using type_traits = conversion_type_traits<c_mag, FromRep, ToRep>;
  using multiplier_type = typename type_traits::multiplier_type;
  // Cast arg to the intermediate computation type; this cast is widening (never truncating).
  const auto arg = static_cast<typename type_traits::c_type>(value);
  // We need to return a representation type used by the user's quantity type. It might have
  // a lower precision than we get as a result of the intermediate scaling calculations.
  // For example, when converting between degree and radian we need to multiply/divide by `pi`
  // which is implemented in terms of `long double`. If the user's quantity type has `double`
  // or `float` representation, this will cause a warning on conversion from `long double`
  // (even with the `static_cast` usage). However, the value truncation is exactly what we want
  // in this case, so we need to suppress the warning here. All such casts go through
  // `silent_cast` which is the single point of intentional truncation in this file.
  if constexpr (is_integral(c_mag))
    return silent_cast<ToRep>(arg * get_value<multiplier_type>(numerator(c_mag)));
  else if constexpr (is_integral(pow<-1>(c_mag)))
    return silent_cast<ToRep>(arg / get_value<multiplier_type>(denominator(c_mag)));
  else {
    using value_traits = conversion_value_traits<c_mag, multiplier_type>;
    if constexpr (std::is_floating_point_v<multiplier_type>)
      // this results in great assembly
      return silent_cast<ToRep>(arg * value_traits::ratio);
    else
      // this is slower but allows conversions like 2000 m -> 2 km without loosing data
      return silent_cast<ToRep>(arg * value_traits::num_mult / value_traits::den_mult * value_traits::irr_mult);
  }
}


/**
 * @brief Explicit cast between different quantity types
 *
 * @note This is a low-level facility and is too powerful to be used by the users directly. They should either use
 * `value_cast` or `quantity_cast`.
 *
 * @tparam To a target quantity type to cast to
 */
template<Quantity To, typename FwdFrom, Quantity From = std::remove_cvref_t<FwdFrom>>
  requires(mp_units::castable(From::quantity_spec, To::quantity_spec)) &&
          (((equivalent(From::unit, To::unit)) && std::constructible_from<typename To::rep, typename From::rep>) ||
           (!equivalent(From::unit, To::unit)))  // && scalable_with_<typename To::rep>))
// TODO how to constrain the second part here?
[[nodiscard]] constexpr To sudo_cast(FwdFrom&& q)
{
  if constexpr (equivalent(From::unit, To::unit)) {
    // no scaling of the number needed
    return {silent_cast<typename To::rep>(std::forward<FwdFrom>(q).numerical_value_is_an_implementation_detail_), To::reference};
  } else {
    return {sudo_cast_value<From::unit, typename From::rep, To::unit, typename To::rep>(
              std::forward<FwdFrom>(q).numerical_value_is_an_implementation_detail_),
            To::reference};
  }
}


/**
 * @brief Explicit cast between different quantity_point types
 *
 * @note This is a low-level facility and is too powerful to be used by the users directly. They should either use
 * `value_cast` or `quantity_cast`.
 *
 * @tparam ToQP a target quantity point type to which to cast to
 */
template<QuantityPoint ToQP, typename FwdFromQP, QuantityPoint FromQP = std::remove_cvref_t<FwdFromQP>>
  requires(mp_units::castable(FromQP::quantity_spec, ToQP::quantity_spec)) &&
          (detail::same_absolute_point_origins(ToQP::point_origin, FromQP::point_origin)) &&
          (((equivalent(FromQP::unit, ToQP::unit)) &&
            std::constructible_from<typename ToQP::rep, typename FromQP::rep>) ||
           (!equivalent(FromQP::unit, ToQP::unit)))
[[nodiscard]] constexpr QuantityPoint auto sudo_cast(FwdFromQP&& qp)
{
  if constexpr (is_same_v<MP_UNITS_NONCONST_TYPE(ToQP::point_origin), MP_UNITS_NONCONST_TYPE(FromQP::point_origin)>) {
    return quantity_point{
      sudo_cast<typename ToQP::quantity_type>(std::forward<FwdFromQP>(qp).quantity_from(FromQP::point_origin)),
      FromQP::point_origin};
  } else {
    // it's unclear how hard we should try to avoid truncation here. For now, the only corner case we cater for,
    // is when the range of the quantity type of at most one of QP or ToQP doesn't cover the offset between the
    // point origins. In that case, we need to be careful to ensure we use the quantity type with the larger range
    // of the two to perform the point_origin conversion.
    // Numerically, we'll potentially need to do three things:
    //  (a) cast the representation type
    //  (b) scale the numerical value
    //  (c) add/subtract the origin difference
    // In the following, we carefully select the order of these three operations: each of (a) and (b) is scheduled
    // either before or after (c), such that (c) acts on the largest range possible among all combination of source
    // and target unit and representation.
    constexpr UnitMagnitude auto c_mag =
      mp_units::get_canonical_unit(FromQP::unit).mag / mp_units::get_canonical_unit(ToQP::unit).mag;
    using type_traits = conversion_type_traits<c_mag, typename FromQP::rep, typename ToQP::rep>;
    using value_traits = conversion_value_traits<c_mag, typename type_traits::multiplier_type>;
    using c_rep_type = typename type_traits::c_rep_type;
    if constexpr (value_traits::num_mult * value_traits::irr_mult > value_traits::den_mult) {
      // original unit had a larger unit magnitude; if we first convert to the common representation but retain the
      // unit, we obtain the largest possible range while not causing truncation of fractional values. This is optimal
      // for the offset computation.
      return sudo_cast<ToQP>(
        sudo_cast<quantity_point<FromQP::reference, FromQP::point_origin, c_rep_type>>(std::forward<FwdFromQP>(qp))
          .point_for(ToQP::point_origin));
    } else {
      // new unit may have a larger unit magnitude; we first need to convert to the new unit (potentially causing
      // truncation, but no more than if we did the conversion later), but make sure we keep the larger of the two
      // representation types. Then, we can perform the offset computation.
      return sudo_cast<ToQP>(
        sudo_cast<quantity_point<make_reference(FromQP::quantity_spec, ToQP::unit), FromQP::point_origin, c_rep_type>>(
          std::forward<FwdFromQP>(qp))
          .point_for(ToQP::point_origin));
    }
  }
}

}  // namespace mp_units::detail
