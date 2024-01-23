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