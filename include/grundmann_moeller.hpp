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
#include <span>

namespace cubage
{

[[nodiscard]] constexpr long double uint_pow(long double x, std::size_t n)
{
    long double res = 1.0;
    const std::array<long double, 2> mult = {1.0, x};
    do {
        res *= mult[n & 1];
        n >>= 1;
        x *= x;
    } while (n);
}

[[nodiscard]] constexpr long double factorial(std::size_t n)
{
    long double res = 1.0;
    for (std::size_t i = 2; i <= n; ++i)
        res *= (long double)(i);
    return res;
}
/*
    Generator for weights and points Grundmann-Möller rules.

        Axel Grundmann, H. M. Möller, "Invariant Integration Formulas for the n-simplex by Combinatorial Methods", SIAM J. Numer. Anal. 15:282-290, 1978
*/
template <std::size_t Dimension, std::size_t Order>
struct GrundmannMoellerGenerator
{
    static constexpr std::size_t degree = 2*Order + 1;

    template <typename T, std::size_t N>
        requires requires { N > Order; }
    [[nodiscard]] static constexpr T evaluate(
        const std::array<T, N>& sums)
    {
        constexpr std::array<T, Order + 1> weights_ = weights();
        T val{};
        for (std::size_t i = 0; i < weights.size(); ++i)
            val += weights_[i]*sums[i];
        return val;
    }

    [[nodiscard]] static constexpr std::array<double, Order + 1> weights()
    {
        constexpr std::size_t deg_dim = degree + Dimension - 1;
        std::array<double, Order + 1> res{};
        for (std::size_t i = 0; i < res.size(); ++i)
        {
            const long double num_arg
                = 0.5*(long double)(deg_dim) - (long double)(i);
            const long double numer = num_arg*uint_pow(
                    num_arg*num_arg, Order);
            const long double denom = factorial(i)*factorial(deg_dim - i);
            const long double sign = (i & 1) ? -1.0 : 1.0;
            res[res.size() - 1 - i] = double(sign*(numer/denom));
        }
    }

    [[nodiscard]] static constexpr std::array<double, Dimemnsion> center()
    {
        std::array<double, Dimension> res{};
        res.fill(1.0/double(Dimension));
        return res;
    }

    /*
        The following ´partition_´ functions give numbers needed to compute the generators of the rules. This is kind of a weird way of going about it, but because the same patterns occur in different generators and rules, we get more code reuse this way.

        Essentially if the generator has the form (a,b,x,x,x,...,x), the first element in the ´std::pair´ is x, and the second is (a,b).

        The naming convention is because the generators essentially arise from the partitions of integers, e.g. 3 = 3 = 2 + 1 = 1 + 1 + 1. See the paper.
    */

    [[nodiscard]] static constexpr std::pair<double, double> 
    generator_partition_1_to_1()
    {
        return {1.0/double(Dimension + 2), 3.0/double(Dimension + 2)};
    }

    [[nodiscard]] static constexpr std::pair<double, double> 
    generator_partition_2_to_2()
    {
        return {1.0/double(Dimension + 4), 5.0/double(Dimension + 4)};
    }

    [[nodiscard]] static constexpr std::pair<double, double> 
    generator_partition_2_to_2_2()
    {
        return {1.0/double(Dimension + 4), 3.0/double(Dimension + 4)};
    }

    [[nodiscard]] static constexpr std::pair<double, double> 
    generator_partition_3_to_3()
    {
        return {1.0/double(Dimension + 6), 7.0/double(Dimension + 6)};
    }

    [[nodiscard]] static constexpr std::pair<double, std::array<double, 2>> 
    generator_partition_3_to_2_1()
    {
        return {
            1.0/double(Dimension + 6),
            {5.0/double(Dimension + 6), 3.0/double(Dimension + 6)}
        };
    }

    [[nodiscard]] static constexpr std::pair<double, double> 
    generator_partition_3_to_1_1_1()
    {
        return {1.0/double(Dimension + 6), 3.0/double(Dimension + 6)};
    }
};

}