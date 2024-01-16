#pragma once

#include "array_arithmetic.hpp"

namespace cubage
{

template <typename T>
concept arithmetic = std::is_arithmetic<T>::value;

template <std::floating_point T, std::size_t Count, std::size_t Dim>
constexpr auto& orthonormalize(std::array<std::array<T, Dim>, Count>& basis)
{
    for (std::size_t i = 0; i < Count - 1; ++i)
    {
        normalize(basis[i]);
        for (std::size_t j = i + 1; j < Count; ++j)
            basis[j] -= project_along_unitv(basis[j], basis[i]);
    }
}

template <std::floating_point T, std::size_t N>
constexpr std::array<T, N> project_along_unitv(
    const std::array<T, N>& a, const std::array<T, N>& uvec)
{
    return dot(a, uvec)*uvec;
}

template <std::floating_point T, std::size_t N>
constexpr std::array<T, N>& normalize(std::array<T, N>& a)
{
    T inv_length = 1.0/length(a);
    a *= inv_length;
    return a;
}

template <arithmetic T, std::size_t N>
constexpr T dot(const std::array<T, N>& a, const std::array<T, N>& b)
{
    T sum = 0;
    for (std::size_t i = 0; i < N; ++i)
        sum += a[i]*b[i];
    return sum;
}

template <arithmetic T, std::size_t N>
constexpr T triple_dot(
    const std::array<T, N>& a, const std::array<T, N>& b,
    const std::array<T, N>& c)
{
    T sum = 0;
    for (std::size_t i = 0; i < N; ++i)
        sum += a[i]*b[i]*c[i];
    return sum;
}

template <std::floating_point T, std::size_t N>
constexpr T length(const std::array<T, N>& a)
{
    T largest = *std::ranges::max_element(a);
    if (largest == 0.0) return 0.0;
    T inv_largest = 1.0/largest;
    T sqr_sum = 0.0;
    for (const auto& element : a)
    {
        T scaled = element*inv_largest;
        sqr_sum += scaled*scaled;
    }
    return largest*std::sqrt(sqr_sum);
}

}