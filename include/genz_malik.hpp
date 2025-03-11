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

#include "box_region.hpp"

namespace cubage
{

template <typename FieldType>
[[nodiscard]] constexpr auto l1_norm(const FieldType& x)
{
    if constexpr (std::is_floating_point<FieldType>::value)
        return std::fabs(x);
    else
    {
        using value_type = typename FieldType::value_type;
        auto v = x | std::views::transform(
                static_cast<double(*)(double)>(std::fabs));
        return std::accumulate(v.begin(), v.end(), value_type{});
    }
}

template <typename FieldType>
concept GenzMalikIntegrable
    = ArrayLike<FieldType>
    && std::tuple_size<FieldType>::value > 1
    && std::tuple_size<FieldType>::value <= 32
    && FloatingPointVectorOperable<FieldType>;

/*
    Genz-Malik rule of degree 7 based on 

        A. C. Genz, A. A. Malik, "Remarks on algorithm 006 : An adaptive 
        algorithm for numerical integration Over an N-dimensional rectangular 
        region", J. Comput. Appl. Math. 6:295-302, 1980

        Jarle Berntsen, Terje O. Espelid, Alan Genz, "An Adaptive Algorithm for 
        the Approximate Calculation of Multiple Integrals", ACM Trans. Math. 
        Software 17:437-451, 1991
    
    This integration rule computes an approximation to an integral in a multi-dimensional hyperrectangular region. The error is estimated by comparing the computed value to a value from a lower order rule based on a subset of the evaluation points of the higher order rule.

    This rule is compatible with an adaptive subdivision strategy, and so 
    returns a suggestion for a subdivision axis. This axis is the axis with the 
    largest fourth difference as described in the second reference.

    The dimensionality of the domain of integration is limited to <= 32. If you 
    need to integrate over a 33 dimensional region, you should probably be 
    using Monte-Carlo.
*/
template <GenzMalikIntegrable DomainTypeParam, typename CodomainTypeParam>
    requires std::floating_point<CodomainTypeParam>
        || (FloatingPointVectorOperable<CodomainTypeParam>
            && ArrayLike<CodomainTypeParam>)
struct GenzMalikD7
{
    using DomainType = DomainTypeParam;
    using CodomainType = CodomainTypeParam;
    using ReturnType = std::pair<IntegralResult<CodomainType>, std::size_t>;
    using Limits = Box<DomainType>;
    using RegionType = SubdivisibleBox<DomainType>;

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] static constexpr ReturnType
    integrate(FuncType f, const Limits& limits)
    {
        constexpr std::size_t ndim = std::tuple_size<DomainType>::value;
        constexpr double v = double(1UL << ndim);
        constexpr std::array<double, 5> gm_weights_d7 = {
            ((400.0/19683.0)*double(ndim) + (-9120.0/19683.0))*double(ndim) + (12824.0/19683.0),
            980.0/6561.0,
            (-400.0/19683.0)*double(ndim) + (1820.0/19683.0),
            200.0/19683.0,
            (6859.0/19683.0)/v
        };

        constexpr std::array<double, 4> gm_weights_d5 = {
            ((50.0/729.0)*double(ndim) + (-950.0/729.0))*double(ndim) + 1.0,
            245.0/486.0,
            (-100.0/1458.0)*double(ndim) + (265.0/1458.0),
            25.0/729.0
        };

        const DomainType center = limits.center();
        const DomainType half_lengths = 0.5*limits.side_lengths();
        const double volume = limits.volume();

        const CodomainType central_value = f(center);
        const CodomainType& gm_sum_1 = central_value;

        constexpr double gm_point_0 = 0.358568582800318091990645153907937495454;
        constexpr double gm_point_1 = 0.948683298050513799599668063329815560116;
        const auto& [gm_sum_2, second_diff_2] = symmetric_sum_1_var(
                f, center, half_lengths, central_value, gm_point_0);
        const auto& [gm_sum_3, second_diff_3] = symmetric_sum_1_var(
                f, center, half_lengths, central_value, gm_point_1);
        const CodomainType gm_sum_4 = symmetric_sum_2_var(
                f, center, half_lengths);
        const CodomainType gm_sum_5 = symmetric_sum_n_var(
                f, center, half_lengths);

        const std::array<double, 5> volume_weights_d7 = {
            volume*gm_weights_d7[0], volume*gm_weights_d7[1], volume*gm_weights_d7[2], volume*gm_weights_d7[3], volume*gm_weights_d7[4]
        };
        const std::array<double, 4> volume_weights_d5 = {
            volume*gm_weights_d5[0], volume*gm_weights_d5[1], volume*gm_weights_d5[2], volume*gm_weights_d5[3]
        };

        const CodomainType val
                = volume_weights_d7[0]*gm_sum_1 + volume_weights_d7[1]*gm_sum_2
                + volume_weights_d7[2]*gm_sum_3 + volume_weights_d7[3]*gm_sum_4
                + volume_weights_d7[4]*gm_sum_5;
        
        const CodomainType test_val
                = volume_weights_d5[0]*gm_sum_1 + volume_weights_d5[1]*gm_sum_2
                + volume_weights_d5[2]*gm_sum_3 + volume_weights_d5[3]*gm_sum_4;
        
        CodomainType err = val - test_val;
        if constexpr (std::is_floating_point<CodomainType>::value)
            err = std::fabs(err);
        else
            std::ranges::transform(
                    err, err.begin(), static_cast<double(*)(double)>(std::fabs));

        const std::array<double, ndim> fourth_diff_normed
                = normed_fourth_difference(second_diff_2, second_diff_3);

        return {
            IntegralResult<CodomainType>{val, err},
            subdiv_axis(fourth_diff_normed)
        };
    }

    [[nodiscard]] static constexpr std::size_t points_count()
    {
        constexpr std::size_t ndim = std::tuple_size<DomainType>::value;
        return (1 << ndim) + 1 + 2*ndim*(1 + ndim);
    }

private:
    using DiffType
            = std::array<CodomainType, std::tuple_size<DomainType>::value>;
    using NormedDiffType
            = std::array<double, std::tuple_size<DomainType>::value>;

    [[nodiscard]] static constexpr std::size_t
    subdiv_axis(const NormedDiffType& fourth_diff_normed)
    {
        return std::size_t(std::distance(
                fourth_diff_normed.begin(),
                std::ranges::max_element(fourth_diff_normed)));
    }

    [[nodiscard]] static constexpr NormedDiffType 
    normed_fourth_difference(
        const DiffType& second_diff_2, const DiffType& second_diff_3)
    {
        constexpr std::size_t ndim = std::tuple_size<DomainType>::value;
        constexpr double ratio = (1.0/7.0);
        
        NormedDiffType fourth_diff_normed;
        for (size_t i = 0; i < ndim; ++i)
            fourth_diff_normed[i] = l1_norm(
                    second_diff_2[i] + ratio*second_diff_3[i]);

        return fourth_diff_normed;
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] static constexpr std::pair<CodomainType, DiffType> 
    symmetric_sum_1_var(
        FuncType f, const DomainType& center, const DomainType& half_lengths, const CodomainType& central_value, double gm_point)
    {
        constexpr std::size_t ndim = std::tuple_size<DomainType>::value;
        CodomainType val{};
        DomainType point = center;

        std::array<CodomainType, ndim> second_differences;
        second_differences.fill((-2.0)*central_value);
        for (std::size_t i = 0; i < ndim; ++i)
        {
            const double disp = gm_point*half_lengths[i];

            point[i] = center[i] + disp;
            const CodomainType fval_plus = f(point);
            val += fval_plus;
            second_differences[i] += fval_plus;
            
            point[i] = center[i] - disp;
            const CodomainType fval_minus = f(point);
            val += fval_minus;
            second_differences[i] += fval_minus;

            point[i] = center[i];
        }

        return {val, second_differences};
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] static constexpr CodomainType
    symmetric_sum_2_var(
        FuncType f, const DomainType& center, const DomainType& half_lengths)
    {
        constexpr double gm_point = 0.9486832980505137995996680633298155601160;
        constexpr std::size_t ndim = std::tuple_size<DomainType>::value;
        CodomainType val{};
        DomainType point = center;

        for (std::size_t i = 0; i < ndim; ++i)
        {
            const double disp1 = gm_point*half_lengths[i];

            point[i] = center[i] + disp1;
            for (std::size_t j = i + 1; j < ndim; ++j)
            {
                const double disp2 = gm_point*half_lengths[j];

                point[j] = center[j] + disp2;
                val += f(point);

                point[j] = center[j] - disp2;
                val += f(point);

                point[j] = center[j];
            }
            
            point[i] = center[i] - disp1;
            for (std::size_t j = i + 1; j < ndim; ++j)
            {
                const double disp2 = gm_point*half_lengths[j];

                point[j] = center[j] + disp2;
                val += f(point);

                point[j] = center[j] - disp2;
                val += f(point);

                point[j] = center[j];
            }

            point[i] = center[i];
        }

        return val;
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] static constexpr CodomainType
    symmetric_sum_n_var(
        FuncType f, const DomainType& center, const DomainType& half_lengths)
    {
        constexpr double gm_point = 0.6882472016116852977216287342936235251269;
        constexpr std::size_t ndim = std::tuple_size<DomainType>::value;
        DomainType point = center + gm_point*half_lengths;

        CodomainType val = f(point);
        unsigned int gray = 0;
        for (std::size_t i = 1; i < (1UL << ndim); ++i)
        {
            const unsigned int dim = ctz((unsigned int) i);
            const unsigned int flipped_bit = 1U << dim;
            gray ^= flipped_bit;
            point[dim] = (gray & flipped_bit) ?
                    center[dim] - gm_point*half_lengths[dim]
                    : center[dim] + gm_point*half_lengths[dim];
            
            val += f(point);
        }

        return val;
    }

    [[nodiscard]] static constexpr unsigned int
    ctz(unsigned int n)
    {
    #if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
        return (unsigned int)(__builtin_ctz(n));
    #else
        static constexpr unsigned int table[32] = {
            0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
        };

        return table[((n & ~n + 1)*0x077CB531U) >> 27];
    #endif
    }
};

[[nodiscard]] constexpr double uintpow(double x, unsigned int n)
{
    double res = 1;
    for (unsigned int i = 0; i < n; ++i)
        res *= x;
    return res;
}

}