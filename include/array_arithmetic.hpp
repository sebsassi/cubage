#pragma once

#include <array>

template <typename T>
concept ArithmeticAssignable = requires(T a, T b)
{
    a += b;
    a -= b;
    a *= b;
    a /= b;
};

template <typename T>
    requires std::same_as<T, std::array<typename T::value_type, std::tuple_size<T>::value>> &&
    ArithmeticAssignable<typename T::value_type>
constexpr T& operator+=(
    T& a, const T& b) noexcept
{
    for (std::size_t i = 0; i < std::tuple_size<T>::value; ++i)
        a[i] += b[i];
    return a;
}

template <typename T>
    requires std::same_as<T, std::array<typename T::value_type, std::tuple_size<T>::value>> &&
    ArithmeticAssignable<typename T::value_type>
[[nodiscard]] constexpr T operator+(
    const T& a, const T& b) noexcept
{
    T res = a;
    res += b;
    return res;
}

template <typename T>
    requires std::same_as<T, std::array<typename T::value_type, std::tuple_size<T>::value>> &&
    ArithmeticAssignable<typename T::value_type>
constexpr T& operator-=(
    T& a, const T& b) noexcept
{
    for (std::size_t i = 0; i < std::tuple_size<T>::value; ++i)
        a[i] -= b[i];
    return a;
}

template <typename T>
    requires std::same_as<T, std::array<typename T::value_type, std::tuple_size<T>::value>> &&
    ArithmeticAssignable<typename T::value_type>
[[nodiscard]] constexpr T operator-(
    const T& a, const T& b) noexcept
{
    T res = a;
    res -= b;
    return res;
}

template <typename T>
    requires std::same_as<T, std::array<typename T::value_type, std::tuple_size<T>::value>> &&
    ArithmeticAssignable<typename T::value_type>
constexpr T& operator*=(T& a, const typename T::value_type& b) noexcept
{
    for (std::size_t i = 0; i < std::tuple_size<T>::value; ++i)
        a[i] *= b;
    return a;
}

template <typename T>
    requires std::same_as<T, std::array<typename T::value_type, std::tuple_size<T>::value>> &&
    ArithmeticAssignable<typename T::value_type>
[[nodiscard]] constexpr T operator*(
    const T& a, const typename T::value_type& b) noexcept
{
    T res = a;
    res *= b;
    return res;
}

template <typename T>
    requires std::same_as<T, std::array<typename T::value_type, std::tuple_size<T>::value>> &&
    ArithmeticAssignable<typename T::value_type>
[[nodiscard]] constexpr T operator*(
    const typename T::value_type& b, const T& a) noexcept
{
    T res = a;
    res *= b;
    return res;
}

template <typename T>
    requires std::same_as<T, std::array<typename T::value_type, std::tuple_size<T>::value>> &&
    ArithmeticAssignable<typename T::value_type>
constexpr T& operator*=(
    T& a, const T& b) noexcept
{
    for (std::size_t i = 0; i < std::tuple_size<T>::value; ++i)
        a[i] *= b[i];
    return a;
}

template <typename T>
    requires std::same_as<T, std::array<typename T::value_type, std::tuple_size<T>::value>> &&
    ArithmeticAssignable<typename T::value_type>
[[nodiscard]] constexpr T operator*(
    const T& a, const T& b) noexcept
{
    T res = a;
    res *= b;
    return res;
}