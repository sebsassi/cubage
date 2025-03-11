/*
Copyright (c) 2024 Sebastian Sassi

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do 
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
*/
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