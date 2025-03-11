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

[[nodiscard]] constexpr long double factorial(std::size_t n)
{
    long double res = 1.0;
    for (std::size_t i = 2; i <= n; ++i)
        res *= (long double)(i);
    return res;
}

/*
    Generator for points and weights of Stroud's degree 5 rule

        A. H. Stroud, "A Fifth Degree Integration Formula for the n-Simplex", SIAM J. Numer. Anal. 6:90-98, 1969
    
    but Genz and Cools give an explicit expression for the weights

        Alan Genz, Ronald Cools, "An Adaptive Numerical Cubature Algorithm for Simplices", 1997
*/
template <std::size_t Dimension, std::size_t Degree>
struct StroudGenerator
{
    [[nodiscard]] static constexpr std::array<double, Degree> weights()
    {
        if constexpr (Degree == 1)
            return {1.0/factorial(Dimension)};
        
        constexpr auto r = compute_ri();
        
        constexpr long double lambda1 = 1.0 - (long double)(Dimension)*r[0];
        constexpr long double lambda2 = 1.0 - (long double)(Dimension)*r[1];

        constexpr long double S2 = compute_S23(lambda1, lambda2);
        constexpr long double S3 = compute_S23(lambda2, lambda1);

        if constexpr (Degree == 3)
            return {
                1.0/factorial(Dimension - 1) - (long double)(Dimension)*(S2 + S3),
                compute_S23_prime(lambda1, lambda2),
                compute_S23_prime(lambda2, lambda1)
            };
        
        constexpr auto u = compute_ui();
        
        constexpr long double delta1 = 0.5*(1.0 - (long double)(Dimension)*u[0]);
        constexpr long double delta2 = 0.5*(1.0 - (long double)(Dimension)*u[1]);

        constexpr long double S4 = compute_S23(delta1, delta2);
        constexpr long double S5 = compute_S23(delta2, delta1);

        if constexpr (Degree == 5)
            return {
                1.0/factorial(Dimension - 1) - (long double)(Dimension)*(S2 + S3 + 0.5*(long double)(Dimension - 1)*(S4 + S5)),
                S2, S3, S4, S5
            };
    }

    [[nodiscard]] static constexpr std::array<double, Dimemnsion> center()
    {
        std::array<double, Dimension> res{};
        res.fill(1.0/double(Dimension));
        return res;
    }

    [[nodiscard]] static constexpr std::array<std::pair<double, double>, 2>
    generator_23()
    {
        constexpr auto r = compute_ri();
        return {
            {r[0], 1.0 - double(Dimension - 1)*r[0]},
            {r[1], 1.0 - double(Dimension - 1)*r[1]}
        };
    }

    [[nodiscard]] static constexpr std::array<std::pair<double, double>, 2>
    generator_45()
    {
        constexpr auto u = compute_ui();
        return {
            {u[0], 0.5*(1.0 - double(Dimension - 2)*u[0])},
            {u[1], 0.5*(1.0 - double(Dimension - 2)*u[1])}
        };
    }

private:
    [[nodiscard]] static constexpr std::array<long double, 2> compute_ri()
    {
        constexpr long double denom = (long double)(Dimension*Dimension + 6*(Dimension - 1));
        constexpr long double a = (long double)(Dimension + 3);
        constexpr long double b = std::sqrt(15.0L);
        return {(a - b)/denom, (a + b)/denom};
    }

    [[nodiscard]] static constexpr std::array<long double, 2> compute_ui()
    {
        constexpr long double denom = (long double)(Dimension*Dimension + 12*Dimension - 24);
        constexpr long double a = (long double)(Dimension + 6);
        constexpr long double b = std::sqrt(60.0L);
        return {(a - b)/denom, (a + b)/denom};
    }

    [[nodiscard]] static constexpr long double compute_S45(
        long double d1, long double d2)
    {
        const long double numer = 2.0 - d2*(long double)(Dimension + 4);
        const long double d1_2 = d1*d1;
        const long double d1_4 = d1_2*d1_2;
        const long double denom = d1_4*(d1 - d2)*factorial(Dimension + 4);
        return numer/denom;
    }

    [[nodiscard]] static constexpr long double compute_S23(
        long double l1, long double l2)
    {
        const long double numer = (long double)(2*(26 - Dimension)) - l2*(long double)((12 - Dimension)*(Dimension + 4));
        const long double l1_2 = l1*l1;
        const long double l1_4 = l1_2*l1_2;
        const long double denom = l1_4*(l1 - l2)*factorial(Dimension + 4);
        return numer/denom;
    }

    [[nodiscard]] static constexpr long double compute_S23_prime(
        long double l1, long double l2)
    {
        const long double numer = 2.0 - l2*(long double)(Dimension + 2);
        const long double l1_2 = l1*l1;
        const long double denom = l1_4*(l1 - l2)*factorial(Dimension + 2);
        return numer/denom;
    }
};

}