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
#include <numeric>

[[nodiscard]] constexpr std::size_t pt_index(std::size_t n, std::size_t k)
{
    return n*(n - 1)/2 + k - 1;
}

template <std::size_t N>
[[nodiscard]] constexpr auto generate_partition_triangle()
{
    std::array<std::size_t, N*(N + 1)/2> triangle{};
    triangle[0] = 1;
    for (std::size_t n = 2; n <= N; ++n)
    {
        triangle[pt_index(n, 1)] = triangle[pt_index(n - 1, 1)];
        for (std::size_t k = 2; k < n; ++k)
        {
            triangle[pt_index(n, k)] = triangle[pt_index(n - 1, k - 1)] + triangle[pt_index(n - k, k)];
        }
        triangle[pt_index(n, n)] = triangle[pt_index(n - 1, n - 1)];
    }
    return triangle;
}

[[nodiscard]] constexpr std::size_t
restricted_partition_number(std::size_t n, std::size_t k)
{
    constexpr auto tri = generate_partition_triangle<10>();
    const auto begin = tri.cbegin() + n*(n - 1)/2;
    const auto end = begin + k + 1;
    return std::accumulate(begin, end, 0);
}

template <std::size_t K>
[[nodiscard]] constexpr std::array<std::size_t, K> 
generate_restricted_partitions(std::size_t n)
{
    
}