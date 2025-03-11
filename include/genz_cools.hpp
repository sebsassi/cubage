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

#include <cmath>
#include <array>
#include <numeric>

#include "grundmann_moeller.hpp"
#include "mysovskikh.hpp"
#include "stroud.hpp"
#include "vector_math.hpp"
#include "simplex_region.hpp"

namespace cubage
{

constexpr std::size_t pt_index(std::size_t n, std::size_t k)
{
    return n*(n - 1)/2 + k - 1;
}

template <std::size_t N>
consteval auto generate_partition_triangle()
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

constexpr std::size_t restricted_partition_number(std::size_t n, std::size_t k)
{
    constexpr auto tri = generate_partition_triangle<10>();
    const auto begin = tri.cbegin() + n*(n - 1)/2;
    const auto end = begin + k + 1;
    return std::accumulate(begin, end, 0);
}

constexpr std::size_t partition_number(std::size_t n)
{
    // OEIS A000041
    constexpr std::array<std::size_t, 50> table = {
        1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 
        385, 490, 627, 792, 1002, 1255, 1575, 1958, 2436, 3010, 3718, 4565, 
        5604, 6842, 8349, 10143, 12310, 14883, 17977, 21637, 26015, 31185, 
        37338, 44583, 53174, 63261, 75175, 89134, 105558, 124754, 147273, 173525
    };
    return table[n];
}

template <std::size_t N>
constexpr std::size_t multinomial_coeff(
    const std::array<std::size_t, N>& counts)
{
    std::size_t res = 1;
    std::size_t sum = counts[0];
    for (std::size_t i = 1; i < N; ++i)
    {
        sum += counts[i];
        res *= binomial_coeff(sum, counts[i]);
    }
    return res;
}

constexpr std::size_t binomial_coeff(std::size_t n, std::size_t k)
{
    std::size_t res = 1;
    for (std::size_t i = 0; i < std::min(k, n - k); ++i)
    {
        res *= n - i;
        res /= i + 1;
    }
    return res;
}

template <typename T>
concept GenzCoolsIntegrable
    = ArrayLike<T>
    && std::tuple_size<T>::value > 1
    && std::tuple_size<T>::value <= 410
    && FloatingPointVectorOperable<T>;

template <std::size_t Order, GenzCoolsIntegrable DomainTypeParam, typename CodomainTypeParam>
    requires requires { 0 < Order; Order <= 4; }
struct GenzCools
{
public:
    using DomainType = DomainTypeParam;
    using CodomainType = CodomainTypeParam;
    using ReturnType = std::pair<IntegralResult<CodomainType>, std::size_t>;
    using Limits = Simplex<DomainType>;


    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] static constexpr ReturnType
    integrate(FuncType f, const Limits& limits)
    {
        const auto symmetric_sums = Rule::symmetric_sums(f);

        const std::array<CodomainType, Order + 1> sums{};
        std::size_t offset;
        for (std::size_t i = 0; i < weights.size(); ++i)
        {
            for (std::size_t j = 0; j < partition_number(i); ++j)
                sums[i] += symmetric_sums[offset + i];
            offset += partition_number(i);
        }
    }

private:
    constexpr std::size_t dimension = std::tuple_size<T>::value;
    constexpr std::size_t num_rules = 2*Order + 1;

    static constexpr std::size_t num_gm_weights()
    {
        std::size_t res = 0;
        for (std::size_t i = 0; i <= Order; ++i)
            res += partition_number(i);
        return res;
    }

    [[nodiscard]] static constexpr std::array<std::array<double, num_weights>, num_rules> 
    weights()
    {
        constexpr std::size_t gm_size = num_gm_weights();
        std::array<std::array<double, num_weights>, num_rules> rule_weights{};
        rule_weights[num_rules - 1][0] = 1.0;
        rule_weights[num_rules - 2][gm_size] = StroudGenerator<dimension, 1>::weights()[0];

        constexpr auto gm3_weights = GrundmannMoellerGenerator<dimension, 1>::weights();
        rule_weights[num_rules - 3][0] = gm3_weights[0];
        rule_weights[num_rules - 3][1] = gm3_weights[1];

        if constexpr (Order > 1)
        {
            constexpr auto stroud3_weights = StroudGenerator<dimension, 3>::weights();
            rule_weights[num_rules - 4][0] = stroud3_weights[0];
            rule_weights[num_rules - 4][gm_size] = stroud3_weights[1];
            rule_weights[num_rules - 4][gm_size + 1] = stroud3_weights[2];

            constexpr auto gm5_weights = GrundmannMoellerGenerator<dimension, 2>::weights();
            rule_weights[num_rules - 5][0] = gm5_weights[0];
            rule_weights[num_rules - 5][1] = gm5_weights[1];
            rule_weights[num_rules - 5][2] = gm5_weights[2];
            rule_weights[num_rules - 5][3] = gm5_weights[2];
        }

        if constexpr (Order > 2)
        {
            constexpr auto stroud5_weights = StroudGenerator<dimension, 5>::weights();
            rule_weights[num_rules - 6][0] = stroud5_weights[0];
            rule_weights[num_rules - 6][gm_size] = stroud5_weights[1];
            rule_weights[num_rules - 6][gm_size + 1] = stroud5_weights[2];
            rule_weights[num_rules - 6][gm_size + 2] = stroud5_weights[3];
            rule_weights[num_rules - 6][gm_size + 3] = stroud5_weights[4];

            constexpr auto gm7_weights = GrundmannMoellerGenerator<dimension, 3>::weights();
            rule_weights[num_rules - 7][0] = gm5_weights[0];
            rule_weights[num_rules - 7][1] = gm5_weights[1];
            rule_weights[num_rules - 7][2] = gm5_weights[2];
            rule_weights[num_rules - 7][3] = gm5_weights[2];
            rule_weights[num_rules - 7][4] = gm5_weights[3];
            rule_weights[num_rules - 7][5] = gm5_weights[3];
            rule_weights[num_rules - 7][6] = gm5_weights[3];
        }

        if constexpr (Order > 3)
        {
            constexpr auto mys_weights = MysovskikhD7Generator<Dimension>::weights();
            rule_weights[num_rules - 8][0] = mys_weights[0];
            rule_weights[num_rules - 8][3] = mys_weights[3];
            rule_weights[num_rules - 8][6] = mys_weights[6];
            rule_weights[num_rules - 8][gm_size + 4] = mys_weights[1];
            rule_weights[num_rules - 8][gm_size + 5] = mys_weights[2];
            rule_weights[num_rules - 8][gm_size + 6] = mys_weights[4];
            rule_weights[num_rules - 8][gm_size + 7] = mys_weights[5];
            rule_weights[num_rules - 8][gm_size + 8] = mys_weights[7];

            constexpr auto gm9_weights = GrundmannMoellerGenerator<dimension, 4>::weights();
            rule_weights[num_rules - 9][0] = gm5_weights[0];
            rule_weights[num_rules - 9][1] = gm5_weights[1];
            rule_weights[num_rules - 9][2] = gm5_weights[2];
            rule_weights[num_rules - 9][3] = gm5_weights[2];
            rule_weights[num_rules - 9][4] = gm5_weights[3];
            rule_weights[num_rules - 9][5] = gm5_weights[3];
            rule_weights[num_rules - 9][6] = gm5_weights[3];
            rule_weights[num_rules - 9][7] = gm5_weights[4];
            rule_weights[num_rules - 9][8] = gm5_weights[4];
            rule_weights[num_rules - 9][9] = gm5_weights[4];
            rule_weights[num_rules - 9][10] = gm5_weights[4];
            rule_weights[num_rules - 9][11] = gm5_weights[4];
        }

        make_null_rules(rule_weights);
        orthonormalize_null_rules(rule_weights, symmetry_factors());

        return rule_weights;
    }

    [[nodiscard]] static constexpr std::array<double, num_weights> 
    symmetry_factors()
    {

    }

    // Turn lower order rules into null rules by subtracting them from the
    // highest order rule.
    static constexpr inline auto& make_null_rules(auto& rules)
    {
        for (std::size_t i = 0; i < num_rules - 1; ++i)
        {
            rules[i] -= rules[num_rules - 1];
            rules[i] *= -1.0;
        }
        return rules;
    }
};

/*
Gram-Schmidt process with two modifications:
    1. Due to symmetry, multiple points have the same weight. For each symmetric collection of points, the weight is stored once in the weight vector. Therfore in taking the dot product each pair of weights is scaled by a symmetry factor, which is the number of points in the corresponding set of points.

    2. Each weight vector is scaled to have norm equal to the weight vector of the base (highest order) integration rule.
*/
template <std::floating_point T, std::size_t Count, std::size_t Dim>
constexpr auto& orthonormalize_null_rules(
    std::array<std::array<T, Dim>, Count>& rules, std::array<T, Dim>& symmetry_factors)
{
    T rule_norm_sqr = triple_dot(rules[0], symmetry_factors, rules[0]);
    T inv_rule_norm_sqr = 1.0/rule_norm_sqr;
    for (std::size_t j = 1; j < Count; ++j)
        rules[j] -= project_rule(
                rules[j], rules[0], symmetry_factors, inv_rule_norm_sqr);
    
    for (std::size_t i = 1; i < Count - 1; ++i)
    {
        normalize_null_rule(rules[i], symmetry_factors, rule_norm_sqr);
        for (std::size_t j = i + 1; j < Count; ++j)
            rules[j] -= project_rule(
                    rules[j], rules[i], symmetry_factors, inv_rule_norm_sqr);
    }
    return rules;
}

template <std::floating_point T, std::size_t N>
constexpr inline std::array<T, N>& normalize_null_rule(
    std::array<T, N>& null_rule, const std::array<T, N>& symmetry_factors,
    T base_rule_norm_sqr)
{
    T null_rule_norm_sqr = triple_dot(null_rule, symmetry_factors, null_rule);
    null_rule *= std::sqrt(base_rule_norm_sqr/null_rule_norm_sqr);
    return null_rule;
}

template <std::floating_point T, std::size_t N>
constexpr inline std::array<T, N> project_rule(
    const std::array<T, N>& a, const std::array<T, N>& b, const std::array<T, N>& symmetry_factors, T inv_bnorm_sqr)
{
    return (inv_bnorm_sqr*triple_dot(a, symmetry_factors, b))*b;
}

}