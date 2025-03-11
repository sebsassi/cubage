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

namespace cubage
{

template <std::size_t Dimension, typename S, typename T, typename FuncType>
[[nodiscard]] static constexpr T symmetric_sum_axxxx(
    double base_val, double extra_val, FuncType f, const std::array<S, Dimension + 1>& vertices)
{
    std::array<double, Dimension + 1> std_point{};
    std_point.fill(base_val);

    T sum{};
    for (std::size_t i = 0; i < Dimension + 1; ++i)
    {
        std_point[i] = extra_val;
        sum += f(affine_transform_point(std_point, vertices));
        std_point[i] = base_val;
    }

    return sum;
}

template <std::size_t Dimension, typename S, typename T, typename FuncType>
[[nodiscard]] static constexpr T symmetric_sum_aaxxx(
    double base_val, double extra_val, FuncType f, const std::array<S, Dimension + 1>& vertices)
{
    std::array<double, Dimension + 1> std_point{};
    std_point.fill(base_val);

    T sum{};
    for (std::size_t i = 0; i < Dimension + 1; ++i)
    {
        std_point[i] = extra_val;
        for (std::size_t j = i + 1; j < Dimension + 1; ++j)
        {
            std_point[j] = extra_val;
            sum += f(affine_transform_point(std_point, vertices));
            std_point[j] = base_val;
        }
        std_point[i] = base_val;
    }

    return sum;
}

template <std::size_t Dimension, typename S, typename T, typename FuncType>
[[nodiscard]] static constexpr T symmetric_sum_abxxx(
    double base_val, const std::array<double, 2>& extra_vals, FuncType f, const std::array<S, Dimension + 1>& vertices)
{

    std::array<double, Dimension + 1> std_point{};
    std_point.fill(base_val);

    T sum{};
    for (std::size_t i = 0; i < Dimension + 1; ++i)
    {
        std_point[i] = extra_vals[0];
        for (std::size_t j = 0; j < i; ++j)
        {
            std_point[j] = extra_vals[1];
            sum += f(affine_transform_point(std_point, vertices));
            std_point[j] = base_val;
        }
        for (std::size_t j = i + 1; j < Dimension + 1; ++j)
        {
            std_point[j] = extra_vals[1];
            sum += f(affine_transform_point(std_point, vertices));
            std_point[j] = base_val;
        }
        std_point[i] = base_val;
    }
}

template <std::size_t Dimension, typename S, typename T, typename FuncType>
[[nodiscard]] static constexpr T symmetric_sum_aaaxx(
    double base_val, double extra_val, FuncType f, const std::array<S, Dimension + 1>& vertices)
{
    std::array<double, Dimension + 1> std_point{};
    std_point.fill(base_val);

    T sum{};
    for (std::size_t i = 0; i < Dimension + 1; ++i)
    {
        std_point[i] = extra_val;
        for (std::size_t j = i + 1; j < Dimension + 1; ++j)
        {
            std_point[j] = extra_val;
            for (std::size_t k = j + 1; k < Dimension + 1; ++k)
            {
                std_point[k] = extra_val;
                sum += f(affine_transform_point(std_point, vertices));
                std_point[k] = base_val;
            }
            std_point[j] = base_val;
        }
        std_point[i] = base_val;
    }

    return sum;
}



template <std::size_t Dimension, typename S>
[[nodiscard]] static constexpr inline S
affine_transform_point(
    const std::array<double, Dimension + 1>& std_point,
    const std::array<S, Dimension + 1>& vertices)
{
    S point{};
    for (std::size_t k = 0; k < Dimension + 1; ++k)
        point += std_point[k]*vertices[k];
    return point;
}

}