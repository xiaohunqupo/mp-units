
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

#include <mp-units/bits/dimension_concepts.h>
#include <mp-units/bits/quantity_concepts.h>
#include <mp-units/bits/quantity_spec_concepts.h>
#include <mp-units/bits/reference_concepts.h>
#include <mp-units/bits/representation_concepts.h>
#include <mp-units/bits/sudo_cast.h>
#include <mp-units/bits/unit_concepts.h>
#include <mp-units/customization_points.h>
#include <mp-units/reference.h>
#include <compare>
#include <utility>

// the below is not used in this header but should be exposed with it
#include <mp-units/bits/quantity_cast.h>
#include <mp-units/bits/value_cast.h>

namespace mp_units {

namespace detail {

template<auto UFrom, auto UTo>
concept IntegralConversionFactor = Unit<decltype(UFrom)> && Unit<decltype(UTo)> &&
                                   is_integral(get_canonical_unit(UFrom).mag / get_canonical_unit(UTo).mag);

template<typename QFrom, typename QTo>
concept QuantityConvertibleTo =
  Quantity<QFrom> && Quantity<QTo> && implicitly_convertible(QFrom::quantity_spec, QTo::quantity_spec) &&
  convertible(QFrom::unit, QTo::unit) &&
  (treat_as_floating_point<typename QTo::rep> ||
   (!treat_as_floating_point<typename QFrom::rep> && IntegralConversionFactor<QFrom::unit, QTo::unit>)) &&
  // TODO consider providing constraints of sudo_cast here rather than testing if it can be called (its return type is
  // deduced thus the function is evaluated here and may emit truncating conversion or other warnings)
  requires(QFrom q) { detail::sudo_cast<QTo>(q); };

template<quantity_character Ch, typename Func, typename T, typename U>
concept InvokeResultOf = std::regular_invocable<Func, T, U> && RepresentationOf<std::invoke_result_t<Func, T, U>, Ch>;

template<typename Func, typename Q1, typename Q2>
concept InvocableQuantities = Quantity<Q1> && Quantity<Q2> &&
                              InvokeResultOf<common_quantity_spec(Q1::quantity_spec, Q2::quantity_spec).character, Func,
                                             typename Q1::rep, typename Q2::rep> &&
                              requires { common_reference(Q1::reference, Q2::reference); };

template<typename Func, Quantity Q1, Quantity Q2>
  requires detail::InvocableQuantities<Func, Q1, Q2>
using common_quantity_for = quantity<common_reference(Q1::reference, Q2::reference),
                                     std::invoke_result_t<Func, typename Q1::rep, typename Q2::rep>>;

}  // namespace detail

/**
 * @brief A quantity
 *
 * Property of a phenomenon, body, or substance, where the property has a magnitude that can be
 * expressed by means of a number and a reference.
 *
 * @tparam R a reference of the quantity providing all information about quantity properties
 * @tparam Rep a type to be used to represent values of a quantity
 */
template<Reference auto R, RepresentationOf<get_quantity_spec(R).character> Rep = double>
class quantity {
public:
  Rep numerical_value_;  // needs to be public for a structural type

  // member types and values
  static constexpr Reference auto reference = R;
  static constexpr QuantitySpec auto quantity_spec = get_quantity_spec(reference);
  static constexpr Dimension auto dimension = quantity_spec.dimension;
  static constexpr Unit auto unit = get_unit(reference);
  using rep = Rep;

  // static member functions
  [[nodiscard]] static constexpr quantity zero() noexcept
    requires requires { quantity_values<rep>::zero(); }
  {
    return quantity(quantity_values<rep>::zero());
  }

  [[nodiscard]] static constexpr quantity one() noexcept
    requires requires { quantity_values<rep>::one(); }
  {
    return quantity(quantity_values<rep>::one());
  }

  [[nodiscard]] static constexpr quantity min() noexcept
    requires requires { quantity_values<rep>::min(); }
  {
    return quantity(quantity_values<rep>::min());
  }

  [[nodiscard]] static constexpr quantity max() noexcept
    requires requires { quantity_values<rep>::max(); }
  {
    return quantity(quantity_values<rep>::max());
  }

  // construction, assignment, destruction
  quantity() = default;
  quantity(const quantity&) = default;
  quantity(quantity&&) = default;

  template<detail::QuantityConvertibleTo<quantity> Q>
  constexpr explicit(!std::convertible_to<typename Q::rep, Rep>) quantity(const Q& q) :
      numerical_value_(detail::sudo_cast<quantity>(q).numerical_value_)
  {
  }

  template<QuantityLike Q>
    requires detail::QuantityConvertibleTo<
      quantity<quantity_like_traits<Q>::reference, typename quantity_like_traits<Q>::rep>, quantity>
  constexpr explicit(is_specialization_of<decltype(quantity_like_traits<Q>::to_numerical_value(std::declval<Q>())),
                                          convert_explicitly> ||
                     !std::convertible_to<typename Q::rep, Rep>) quantity(const Q& q) :
      quantity(make_quantity<quantity_like_traits<Q>::reference>(quantity_like_traits<Q>::to_numerical_value(q).value))
  {
  }

  quantity& operator=(const quantity&) = default;
  quantity& operator=(quantity&&) = default;

  // unit conversions
  template<UnitCompatibleWith<unit, quantity_spec> U>
    requires detail::QuantityConvertibleTo<quantity, quantity<detail::make_reference(quantity_spec, U{}), Rep>>
  [[nodiscard]] constexpr Quantity auto in(U) const
  {
    return quantity<detail::make_reference(quantity_spec, U{}), Rep>{*this};
  }

  template<UnitCompatibleWith<unit, quantity_spec> U>
    requires requires(quantity q) { value_cast<U{}>(q); }
  [[nodiscard]] constexpr Quantity auto force_in(U) const
  {
    return value_cast<U{}>(*this);
  }

  // data access
  template<Unit U>
    requires(U{} == unit)
  [[nodiscard]] constexpr rep& numerical_value_ref_in(U) & noexcept
  {
    return numerical_value_;
  }

  template<Unit U>
    requires(U{} == unit)
  [[nodiscard]] constexpr const rep& numerical_value_ref_in(U) const& noexcept
  {
    return numerical_value_;
  }

  template<Unit U>
    requires(U{} == unit)
  constexpr const rep&& numerical_value_ref_in(U) const&& noexcept = delete;

  template<UnitCompatibleWith<unit, quantity_spec> U>
    requires requires(quantity q) { q.in(U{}); }
  [[nodiscard]] constexpr rep numerical_value_in(U) const noexcept
  {
    return (*this).in(U{}).numerical_value_;
  }

  template<UnitCompatibleWith<unit, quantity_spec> U>
    requires requires(quantity q) { q.force_in(U{}); }
  [[nodiscard]] constexpr rep force_numerical_value_in(U) const noexcept
  {
    return (*this).force_in(U{}).numerical_value_;
  }

  // conversion operators
  template<typename Q_, QuantityLike Q = std::remove_cvref_t<Q_>>
    requires detail::QuantityConvertibleTo<
      quantity, quantity<quantity_like_traits<Q>::reference, typename quantity_like_traits<Q>::rep>>
  [[nodiscard]] explicit(is_specialization_of<decltype(quantity_like_traits<Q>::from_numerical_value(numerical_value_)),
                                              convert_explicitly> ||
                         !std::convertible_to<Rep, typename Q::rep>) constexpr
  operator Q_() const& noexcept(noexcept(quantity_like_traits<Q>::from_numerical_value(numerical_value_)) &&
                                std::is_nothrow_copy_constructible_v<rep>)
  {
    return quantity_like_traits<Q>::from_numerical_value(numerical_value_).value;
  }

  template<typename Q_, QuantityLike Q = std::remove_cvref_t<Q_>>
    requires detail::QuantityConvertibleTo<
      quantity, quantity<quantity_like_traits<Q>::reference, typename quantity_like_traits<Q>::rep>>
  [[nodiscard]] explicit(is_specialization_of<decltype(quantity_like_traits<Q>::from_numerical_value(numerical_value_)),
                                              convert_explicitly> ||
                         !std::convertible_to<Rep, typename Q::rep>) constexpr
  operator Q_() && noexcept(noexcept(quantity_like_traits<Q>::from_numerical_value(numerical_value_)) &&
                            std::is_nothrow_move_constructible_v<rep>)
  {
    return quantity_like_traits<Q>::from_numerical_value(std::move(numerical_value_)).value;
  }

  // member unary operators
  [[nodiscard]] constexpr Quantity auto operator+() const
    requires requires(rep v) {
      {
        +v
      } -> std::common_with<rep>;
    }
  {
    return make_quantity<reference>(+numerical_value_);
  }

  [[nodiscard]] constexpr Quantity auto operator-() const
    requires requires(rep v) {
      {
        -v
      } -> std::common_with<rep>;
    }
  {
    return make_quantity<reference>(-numerical_value_);
  }

  template<typename Q>
  friend constexpr decltype(auto) operator++(Q&& q)
    requires std::derived_from<std::remove_cvref_t<Q>, quantity> && requires(rep v) {
      {
        ++v
      } -> std::same_as<rep&>;
    }
  {
    ++q.numerical_value_;
    return std::forward<Q>(q);
  }

  [[nodiscard]] constexpr Quantity auto operator++(int)
    requires requires(rep v) {
      {
        v++
      } -> std::common_with<rep>;
    }
  {
    return make_quantity<reference>(numerical_value_++);
  }

  template<typename Q>
  friend constexpr decltype(auto) operator--(Q&& q)
    requires std::derived_from<std::remove_cvref_t<Q>, quantity> && requires(rep v) {
      {
        --v
      } -> std::same_as<rep&>;
    }
  {
    --q.numerical_value_;
    return std::forward<Q>(q);
  }

  [[nodiscard]] constexpr Quantity auto operator--(int)
    requires requires(rep v) {
      {
        v--
      } -> std::common_with<rep>;
    }
  {
    return make_quantity<reference>(numerical_value_--);
  }

  // compound assignment operators
  template<typename Q>
    requires std::derived_from<std::remove_cvref_t<Q>, quantity> && requires(rep a, rep b) {
      {
        a += b
      } -> std::same_as<rep&>;
    }
  friend constexpr decltype(auto) operator+=(Q&& lhs, const quantity& rhs)
  {
    lhs.numerical_value_ += rhs.numerical_value_;
    return std::forward<Q>(lhs);
  }

  template<typename Q>
    requires std::derived_from<std::remove_cvref_t<Q>, quantity> && requires(rep a, rep b) {
      {
        a -= b
      } -> std::same_as<rep&>;
    }
  friend constexpr decltype(auto) operator-=(Q&& lhs, const quantity& rhs)
  {
    lhs.numerical_value_ -= rhs.numerical_value_;
    return std::forward<Q>(lhs);
  }

  template<typename Q>
    requires std::derived_from<std::remove_cvref_t<Q>, quantity> && (!treat_as_floating_point<rep>) &&
             requires(rep a, rep b) {
               {
                 a %= b
               } -> std::same_as<rep&>;
             }
  friend constexpr decltype(auto) operator%=(Q&& lhs, const quantity& rhs)

  {
    gsl_ExpectsAudit(rhs != zero());
    lhs.numerical_value_ %= rhs.numerical_value_;
    return std::forward<Q>(lhs);
  }

  template<typename Q, typename Value>
    requires std::derived_from<std::remove_cvref_t<Q>, quantity> && (!Quantity<Value>) &&
             requires(rep a, const Value b) {
               {
                 a *= b
               } -> std::same_as<rep&>;
             }
  friend constexpr decltype(auto) operator*=(Q&& lhs, const Value& v)
  {
    lhs.numerical_value_ *= v;
    return std::forward<Q>(lhs);
  }

  template<typename Q1, QuantityOf<dimension_one> Q2>
    requires std::derived_from<std::remove_cvref_t<Q1>, quantity> && (Q2::unit == ::mp_units::one) &&
             requires(rep a, const typename Q2::rep b) {
               {
                 a *= b
               } -> std::same_as<rep&>;
             }
  friend constexpr decltype(auto) operator*=(Q1&& lhs, const Q2& rhs)
  {
    lhs.numerical_value_ *= rhs.numerical_value_;
    return std::forward<Q1>(lhs);
  }

  template<typename Q, typename Value>
    requires std::derived_from<std::remove_cvref_t<Q>, quantity> && (!Quantity<Value>) &&
             requires(rep a, const Value b) {
               {
                 a /= b
               } -> std::same_as<rep&>;
             }
  friend constexpr decltype(auto) operator/=(Q&& lhs, const Value& v)
  {
    gsl_ExpectsAudit(v != quantity_values<Value>::zero());
    lhs.numerical_value_ /= v;
    return std::forward<Q>(lhs);
  }

  template<typename Q1, QuantityOf<dimension_one> Q2>
    requires std::derived_from<std::remove_cvref_t<Q1>, quantity> && (Q2::unit == ::mp_units::one) &&
             requires(rep a, const typename Q2::rep b) {
               {
                 a /= b
               } -> std::same_as<rep&>;
             }
  friend constexpr decltype(auto) operator/=(Q1&& lhs, const Q2& rhs)
  {
    gsl_ExpectsAudit(rhs != rhs.zero());
    lhs.numerical_value_ /= rhs.numerical_value_;
    return std::forward<Q1>(lhs);
  }

private:
#if defined MP_UNITS_COMP_CLANG && MP_UNITS_COMP_CLANG < 17
  template<auto R2, typename Rep2>
#else
  template<Reference auto R2, typename Rep2>
#endif
  friend constexpr quantity<R2, std::remove_cvref_t<Rep2>> make_quantity(Rep2&&);

  template<typename Value>
    requires std::constructible_from<rep, Value&&>
  constexpr explicit quantity(Value&& v) : numerical_value_(std::forward<Value>(v))
  {
  }
};

// CTAD
template<QuantityLike Q>
explicit(
  is_specialization_of<decltype(quantity_like_traits<Q>::to_numerical_value(std::declval<Q>())), convert_explicitly>)
  quantity(Q) -> quantity<quantity_like_traits<Q>::reference, typename quantity_like_traits<Q>::rep>;

// binary operators on quantities
template<auto R1, typename Rep1, auto R2, typename Rep2>
  requires detail::InvocableQuantities<std::plus<>, quantity<R1, Rep1>, quantity<R2, Rep2>>
[[nodiscard]] constexpr Quantity auto operator+(const quantity<R1, Rep1>& lhs, const quantity<R2, Rep2>& rhs)
{
  using ret = detail::common_quantity_for<std::plus<>, quantity<R1, Rep1>, quantity<R2, Rep2>>;
  return make_quantity<ret::reference>(ret(lhs).numerical_value_ + ret(rhs).numerical_value_);
}

template<auto R1, typename Rep1, auto R2, typename Rep2>
  requires detail::InvocableQuantities<std::minus<>, quantity<R1, Rep1>, quantity<R2, Rep2>>
[[nodiscard]] constexpr Quantity auto operator-(const quantity<R1, Rep1>& lhs, const quantity<R2, Rep2>& rhs)
{
  using ret = detail::common_quantity_for<std::minus<>, quantity<R1, Rep1>, quantity<R2, Rep2>>;
  return make_quantity<ret::reference>(ret(lhs).numerical_value_ - ret(rhs).numerical_value_);
}

template<auto R1, typename Rep1, auto R2, typename Rep2>
  requires(!treat_as_floating_point<Rep1>) && (!treat_as_floating_point<Rep2>) &&
          detail::InvocableQuantities<std::modulus<>, quantity<R1, Rep1>, quantity<R2, Rep2>>
[[nodiscard]] constexpr Quantity auto operator%(const quantity<R1, Rep1>& lhs, const quantity<R2, Rep2>& rhs)
{
  gsl_ExpectsAudit(rhs != rhs.zero());
  using ret = detail::common_quantity_for<std::modulus<>, quantity<R1, Rep1>, quantity<R2, Rep2>>;
  return make_quantity<ret::reference>(ret(lhs).numerical_value_ % ret(rhs).numerical_value_);
}

template<auto R1, typename Rep1, auto R2, typename Rep2>
  requires detail::InvokeResultOf<(get_quantity_spec(R1) * get_quantity_spec(R2)).character, std::multiplies<>, Rep1,
                                  Rep2>
[[nodiscard]] constexpr Quantity auto operator*(const quantity<R1, Rep1>& lhs, const quantity<R2, Rep2>& rhs)
{
  return make_quantity<R1 * R2>(lhs.numerical_value_ * rhs.numerical_value_);
}

template<auto R, typename Rep, typename Value>
  requires(!Quantity<Value>) && (!Reference<Value>) &&
          detail::InvokeResultOf<get_quantity_spec(R).character, std::multiplies<>, Rep, const Value&>
[[nodiscard]] constexpr Quantity auto operator*(const quantity<R, Rep>& q, const Value& v)
{
  return make_quantity<R>(q.numerical_value_ * v);
}

template<typename Value, auto R, typename Rep>
  requires(!Quantity<Value>) && (!Reference<Value>) &&
          detail::InvokeResultOf<get_quantity_spec(R).character, std::multiplies<>, const Value&, Rep>
[[nodiscard]] constexpr Quantity auto operator*(const Value& v, const quantity<R, Rep>& q)
{
  return make_quantity<R>(v * q.numerical_value_);
}

template<auto R1, typename Rep1, auto R2, typename Rep2>
  requires detail::InvokeResultOf<(get_quantity_spec(R1) / get_quantity_spec(R2)).character, std::divides<>, Rep1, Rep2>
[[nodiscard]] constexpr Quantity auto operator/(const quantity<R1, Rep1>& lhs, const quantity<R2, Rep2>& rhs)
{
  gsl_ExpectsAudit(rhs != rhs.zero());
  return make_quantity<R1 / R2>(lhs.numerical_value_ / rhs.numerical_value_);
}

template<auto R, typename Rep, typename Value>
  requires(!Quantity<Value>) && (!Reference<Value>) &&
          detail::InvokeResultOf<get_quantity_spec(R).character, std::divides<>, Rep, const Value&>
[[nodiscard]] constexpr Quantity auto operator/(const quantity<R, Rep>& q, const Value& v)
{
  gsl_ExpectsAudit(v != quantity_values<Value>::zero());
  return make_quantity<R>(q.numerical_value_ / v);
}

template<typename Value, auto R, typename Rep>
  requires(!Quantity<Value>) && (!Reference<Value>) &&
          detail::InvokeResultOf<get_quantity_spec(R).character, std::divides<>, const Value&, Rep>
[[nodiscard]] constexpr Quantity auto operator/(const Value& v, const quantity<R, Rep>& q)
{
  return make_quantity<::mp_units::one / R>(v / q.numerical_value_);
}

template<auto R1, typename Rep1, auto R2, typename Rep2>
  requires requires { typename std::common_type_t<quantity<R1, Rep1>, quantity<R2, Rep2>>; } &&
           std::equality_comparable<typename std::common_type_t<quantity<R1, Rep1>, quantity<R2, Rep2>>::rep>
[[nodiscard]] constexpr bool operator==(const quantity<R1, Rep1>& lhs, const quantity<R2, Rep2>& rhs)
{
  using ct = std::common_type_t<quantity<R1, Rep1>, quantity<R2, Rep2>>;
  return ct(lhs).numerical_value_ == ct(rhs).numerical_value_;
}

template<auto R1, typename Rep1, auto R2, typename Rep2>
  requires requires { typename std::common_type_t<quantity<R1, Rep1>, quantity<R2, Rep2>>; } &&
           std::three_way_comparable<typename std::common_type_t<quantity<R1, Rep1>, quantity<R2, Rep2>>::rep>
[[nodiscard]] constexpr auto operator<=>(const quantity<R1, Rep1>& lhs, const quantity<R2, Rep2>& rhs)
{
  using ct = std::common_type_t<quantity<R1, Rep1>, quantity<R2, Rep2>>;
  return ct(lhs).numerical_value_ <=> ct(rhs).numerical_value_;
}

// make_quantity
#if defined MP_UNITS_COMP_CLANG && MP_UNITS_COMP_CLANG < 17
template<auto R, typename Rep>
#else
template<Reference auto R, typename Rep>
#endif
[[nodiscard]] constexpr quantity<R, std::remove_cvref_t<Rep>> make_quantity(Rep&& v)
{
  return quantity<R, std::remove_cvref_t<Rep>>(std::forward<Rep>(v));
}

}  // namespace mp_units

namespace std {

template<mp_units::Quantity Q1, mp_units::Quantity Q2>
  requires requires {
    {
      mp_units::common_reference(Q1::reference, Q2::reference)
    } -> mp_units::Reference;
    typename common_type_t<typename Q1::rep, typename Q2::rep>;
  }
struct common_type<Q1, Q2> {
public:
  using type = mp_units::quantity<mp_units::common_reference(Q1::reference, Q2::reference),
                                  common_type_t<typename Q1::rep, typename Q2::rep>>;
};

// dimensionless quantities support for interacting with raw values
template<mp_units::Quantity Q, typename Value>
  requires(!mp_units::Quantity<Value>) && (Q::dimension == mp_units::dimension_one) && (Q::unit == mp_units::one) &&
          requires { typename common_type_t<typename Q::rep, Value>; }
struct common_type<Q, Value> {
public:
  using type = mp_units::quantity<mp_units::one, common_type_t<typename Q::rep, Value>>;
};

template<mp_units::Quantity Q, typename Value>
  requires(!mp_units::Quantity<Value>) && (Q::dimension == mp_units::dimension_one) && (Q::unit == mp_units::one) &&
          requires { typename common_type_t<typename Q::rep, Value>; }
struct common_type<Value, Q> : common_type<Q, Value> {};

}  // namespace std
