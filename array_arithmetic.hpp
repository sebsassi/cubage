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

template <ArithmeticAssignable T, std::size_t N>
constexpr std::array<T, N>& operator+=(
    std::array<T, N>& a, const std::array<T, N>& b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] += b[i];
    return a;
}

template <ArithmeticAssignable T, std::size_t N>
[[nodiscard]] constexpr std::array<T, N> operator+(
    const std::array<T, N>& a, const std::array<T, N>& b)
{
    std::array<T, N> res = a;
    res += b;
    return res;
}

template <ArithmeticAssignable T, std::size_t N>
constexpr std::array<T, N>& operator-=(
    std::array<T, N>& a, const std::array<T, N>& b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] -= b[i];
    return a;
}

template <ArithmeticAssignable T, std::size_t N>
[[nodiscard]] constexpr std::array<T, N> operator-(
    const std::array<T, N>& a, const std::array<T, N>& b)
{
    std::array<T, N> res = a;
    res -= b;
    return res;
}

template <ArithmeticAssignable T, std::size_t N>
constexpr std::array<T, N>& operator*=(std::array<T, N>& a, const T& b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] *= b;
    return a;
}

template <ArithmeticAssignable T, std::size_t N>
[[nodiscard]] constexpr std::array<T, N> operator*(
    const std::array<T, N>& a, const T& b)
{
    std::array<T, N> res = a;
    res *= b;
    return res;
}

template <ArithmeticAssignable T, std::size_t N>
[[nodiscard]] constexpr std::array<T, N> operator*(
    const T& b, const std::array<T, N>& a)
{
    std::array<T, N> res = a;
    res *= b;
    return res;
}

template <ArithmeticAssignable T, std::size_t N>
constexpr std::array<T, N>& operator*=(
    std::array<T, N>& a, const std::array<T, N>& b)
{
    for (std::size_t i = 0; i < N; ++i)
        a[i] *= b[i];
    return a;
}

template <ArithmeticAssignable T, std::size_t N>
[[nodiscard]] constexpr std::array<T, N> operator*(
    const std::array<T, N>& a, const std::array<T, N>& b)
{
    std::array<T, N> res = a;
    res *= b;
    return res;
}