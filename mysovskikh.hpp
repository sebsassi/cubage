#pragma once

#include <array>
#include <cmath>

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

[[nodiscard]] constexpr std::pair<std::array<long double, 3>, std::size_t>
cubic_solution(const std::array<long double, 4>& coeffs)
{
    const long double defactor = 3.0*coeffs[3];
    const long double shift = coeffs[2]/defactor;
    const long double b2 = coeffs[2]*coeffs[2];
    const long double p = (defactor*coeffs[1] - b2)/(defactor*coeffs[3]);
    const long double df2 = defactor*defactor;
    const long double q = (coeffs[2]*(2.0*b2 - 9.0*coeffs[3]*coeffs[1]) + 3.0*df2*coeffs[0])/(df2*defactor);

    const double discriminant = -(4.0*p*p*p + 27.0*q*q);
    if (discriminant < 0)
    {
        if (p < 0)
        {
            const long double q_sign = (q > 0) ? 1.0 : -1.0;
            const long double abs_q = std::fabs(q);
            const long double g = -3.0/p;
            const long double arg = -0.5*g*abs_q*std::sqrt(g);
            const long double t = -2.0*q_sign*std::sqrt(-(1.0/3.0)*p)*std::cosh((1.0/3.0)*std::acosh(arg));
            return {{t - shift, 0.0, 0.0}, 1};
        }
        else
        {
            const long double g = 3.0/p;
            const long double arg = 0.5*q*g*std::sqrt(g);
            const long double t = -2.0*std::sqrt((1.0/3.0)*p)*std::sinh((1.0/3.0)*std::asinh(arg));
            return {{t - shift, 0.0, 0.0}, 1};
        }
    }
    else
    {
        const long double g = -3.0/p;
        const long double acos_arg = 0.5*q*g*std::sqrt(g);
        const long double cos_arg = (1.0/3.0)*acos(acos_arg);
        const long double prefactor = 2*std::sqrt((-1.0/3.0)*p);
        return {
            {
                prefactor*std::cos(cos_arg) - shift,
                prefactor*std::cos(cos_arg - (2.0*M_PI/3.0)) - shift,
                prefactor*std::cos(cos_arg - (4.0*M_PI/3.0)) - shift,
            }, 
            3
        };
    }
}

/*
    Generator for weights and points of a degree 7 rule by Mysovskikh.

        I. P. Mysovskikh, "On a cubature formula for the simplex" (in Russian), Vopros. Vycisl. i Prikl. Mat., Tashkent 51:74-90, 1975
    
    Old soviet era papers aren't generally avaiable on the internet, but Genz and Cools give a description of the rule

        Alan Genz, Ronald Cools, "An Adaptive Numerical Cubature Algorithm for Simplices", 1997
*/
template <std::size_t Dimension>
    requires requires { Dimension < 411; }
struct MysovskikhD7Generator
{
    static constexpr std::size_t degree = 7;

    [[nodiscard]] static constexpr std::array<double, 8> weights()
    {
        constexpr long double ndim = Dimension;
        constexpr long long int six3 = 6LL*6LL*6LL;
        constexpr long long int six4 = six3*6LL;
        constexpr long long int six5 = six4*6LL;

        constexpr long double U5 = (long double)(-six3*(52212LL - Dimension*(6353LL + Dimension*(1934LL - 27LL*Dimension))))/(23328.0*factorial(ndim + 6.0));
        constexpr long double U6 = (long double)(-six4*(7884L - Dimension*(1541L - 9L*Dimension)))/(23328.0*factorial(ndim + 6.0));
        constexpr long double U7 = (long double)(-six5*(8292L - Dimension*(1139L - 3L*Dimension)))/(23328.0*factorial(ndim + 7.0));

        constexpr std::array<long double, 3> a = compute_ai();
        std::array<double, 8> res{};
        res[1] = weight235(a[0], a[1], a[2]);
        res[2] = weight235(a[1], a[0], a[2]);
        res[3] = -2.0*unit_pow(0.5*ndim + 2.5, 7)/factorial(ndim + 6.0);
        res[4] = weight235(a[2], a[0], a[1]);
        res[5] = 30.0*unit_pow((1.0/3.0)*ndim + 7.0/3.0, 7)/factorial(ndim + 7.0);
        res[6] = -2.0*unit_pow(0.5*ndim + 3.5, 7)/factorial(ndim + 7.0);
        res[7] = (1.0/6.0)*unit_pow((2.0/3.0)*ndim + 14.0/3.0, 7)/factorial(ndim + 7.0);
        res[0] = 1.0/factorial(ndim) - (ndim + 1.0)*(res[1] + res[2] + res[4] + 0.5*ndim*(res[3] + res[5] + 2.0*res[7] + 0.5*(ndim - 1.0)*res[6]));
        return res;
    }

    [[nodiscard]] static constexpr std::array<double, Dimemnsion> center()
    {
        std::array<double, Dimension> res{};
        res.fill(1.0/double(Dimension));
        return res;
    }

    [[nodiscard]] static constexpr std::array<std::pair<double, double>, 3> 
    generator_235()
    {
        constexpr long double ndim = Dimension;
        constexpr std::array<double, 3> alpha = compute_alphai();

        return {
            {alpha[0], 1 - ndim*alpha[0]},
            {alpha[1], 1 - ndim*alpha[1]},
            {alpha[2], 1 - ndim*alpha[2]},
        };
    }

    [[nodiscard]] static constexpr std::pair<double, double>
    generator_4()
    {
        return {1.0/double(Dimension + 4), 3.0/double(Dimension + 4)};
    }

    [[nodiscard]] static constexpr std::pair<double, double>
    generator_6()
    {
        return {1.0/double(Dimension + 6), 4.0/double(Dimension + 6)};
    }

    [[nodiscard]] static constexpr std::pair<double, double>
    generator_7()
    {
        return {1.0/double(Dimension + 6), 3.0/double(Dimension + 6)};
    }

    [[nodiscard]] static constexpr std::pair<double, std::array<double, 2>> 
    generator_8()
    {
        return {
            1.0/double(Dimension + 6),
            {5.5/double(Dimension + 6), 2.5/double(Dimension + 6)}
        };
    }


private:
    [[nodiscard]] static constexpr long double weight235(
        long double a, long double b, long double c)
    {
        constexpr long double ndim = Dimension;
        constexpr long long int six3 = 6LL*6LL*6LL;
        constexpr long long int six4 = six3*6LL;
        constexpr long long int six5 = six4*6LL;

        constexpr long double U5 = (long double)(-six3*(52212LL - Dimension*(6353LL + Dimension*(1934LL - 27LL*Dimension))))/(23328.0*factorial(ndim + 6.0));
        constexpr long double U6 = (long double)(-six4*(7884L - Dimension*(1541L - 9L*Dimension)))/(23328.0*factorial(ndim + 6.0));
        constexpr long double U7 = (long double)(-six5*(8292L - Dimension*(1139L - 3L*Dimension)))/(23328.0*factorial(ndim + 7.0));

        const long double a2 = a*a;
        const long double a4 = a2*a2;
        const long double a5 = a4*a;
        return (U7 - (b + c)*U6 + b*c*U5)/(a5*(a2 - (b + c)*a + b*c));
    }

    [[nodiscard]] static constexpr std::array<long double, 3> compute_ai()
    {
        constexpr long double ndim = Dimension;
        constexpr std::array<long double, 4> coeffs = {
            -144.0*(142528.0 + ndim*(23073.0 - 115.0*ndim)),
            -12.0*(6690556.0 + ndim*(2641189.0 + ndim*(245378.0 - 1495.0*ndim))),
            -16.0*(6503401.0 + ndim*(4020794.0 + ndim*(787281.0 + ndim*(47323.0 - 385.0*ndim)))),
            -(6386660.0 + ndim*(4411997.0 + ndim*(951821.0 + ndim*(61659.0 - 665.0*ndim))))*(ndim + 7.0)
        };

        return cubic_solution(coeffs).first;
    }

    [[nodiscard]] static constexpr std::array<long double, 3> compute_alphai()
    {
        return (compute_ai() + 1.0)/(long double)(Dimension);
    }

};

}