#include <iostream>
#include <cmath>

#include "gauss_kronrod.hpp"

constexpr bool close(double a, double b, double tol)
{
    return std::fabs(a - b) < tol;
}

constexpr bool
gauss_kronrod_15_integrates_22nd_degree_polynomial_exactly()
{
    using Rule = cubage::GaussKronrod<double, double, 15>;
    constexpr cubage::Interval<double> limits = {0.0, 1.0};

    auto polynomial = [](double x)
    {
        const double x2 = x*x;
        const double x4 = x2*x2;
        const double x8 = x4*x4;
        const double x10 = x8*x2;
        const double x11 = x10*x;
        const double x22 = x11*x11;
        const double x21 = x11*x10;
        const double x7 = x4*x2*x;
        const double x14 = x10*x4;
        const double x3 = x2*x;
        return x22 + x21 + x14 + x11 + x7 + x4 + x3 + x;
    };

    const auto res = Rule::integrate(polynomial, limits);
    return close(res.val, 39891.0/30360.0, 1.0e-13);
}

constexpr bool
gauss_kronrod_21_integrates_31st_degree_polynomial_exactly()
{
    using Rule = cubage::GaussKronrod<double, double, 21>;
    constexpr cubage::Interval<double> limits = {0.0, 1.0};

    auto polynomial = [](double x)
    {
        const double x2 = x*x;
        const double x4 = x2*x2;
        const double x8 = x4*x4;
        const double x16 = x8*x8;
        const double x24 = x16*x8;
        const double x28 = x24*x4;
        const double x30 = x28*x2;
        const double x31 = x30*x;
        const double x10 = x8*x2;
        const double x5 = x4*x;
        return x31 + x30 + x28 + x24 + x10 + x5 + x4 + x2 + x;
    };

    const auto res = Rule::integrate(polynomial, limits);
    return close(res.val, 11304313.0/7911200.0, 1.0e-13);
}

constexpr bool
gauss_kronrod_31_integrates_46th_degree_polynomial_exactly()
{
    using Rule = cubage::GaussKronrod<double, double, 31>;
    constexpr cubage::Interval<double> limits = {0.0, 1.0};

    auto polynomial = [](double x)
    {
        const double x2 = x*x;
        const double x4 = x2*x2;
        const double x8 = x4*x4;
        const double x16 = x8*x8;
        const double x32 = x16*x16;
        const double x40 = x32*x8;
        const double x44 = x40*x4;
        const double x46 = x44*x2;
        const double x41 = x40*x;
        const double x35 = x32*x2*x;
        const double x21 = x16*x4*x;
        const double x10 = x8*x2;
        const double x6 = x4*x2;
        const double x7 = x6*x;
        return x46 + x41 + x35 + x21 + x10 + x8 + x7 + x6 + x2;
    };

    const auto res = Rule::integrate(polynomial, limits);
    return close(res.val, 34303.0/37224.0, 1.0e-13);
}



constexpr bool
gauss_kronrod_41_integrates_61st_degree_polynomial_exactly()
{
    using Rule = cubage::GaussKronrod<double, double, 41>;
    constexpr cubage::Interval<double> limits = {0.0, 1.0};

    auto polynomial = [](double x)
    {
        const double x2 = x*x;
        const double x4 = x2*x2;
        const double x8 = x4*x4;
        const double x16 = x8*x8;
        const double x32 = x16*x16;
        const double x48 = x32*x16;
        const double x56 = x48*x8;
        const double x60 = x56*x4;
        const double x61 = x60*x;
        const double x49 = x48*x;
        const double x40 = x32*x8;
        const double x24 = x16*x8;
        const double x26 = x24*x2;
        const double x27 = x26*x;
        const double x9 = x8*x;
        const double x3 = x2*x;
        return x61 + x56 + x49 + x40 + x32 + x27 + x24 + x9 + x3;
    };

    const auto res = Rule::integrate(polynomial, limits);
    return close(res.val, 297932454.0/557841900.0, 1.0e-13);
}



constexpr bool
gauss_kronrod_51_integrates_76th_degree_polynomial_exactly()
{
    using Rule = cubage::GaussKronrod<double, double, 51>;
    constexpr cubage::Interval<double> limits = {0.0, 1.0};

    auto polynomial = [](double x)
    {
        const double x2 = x*x;
        const double x4 = x2*x2;
        const double x8 = x4*x4;
        const double x16 = x8*x8;
        const double x32 = x16*x16;
        const double x64 = x32*x32;
        const double x72 = x64*x8;
        const double x76 = x72*x4;
        const double x48 = x32*x16;
        const double x40 = x32*x8;
        const double x24 = x16*x8;
        const double x25 = x24*x;
        const double x20 = x16*x4;
        const double x10 = x8*x2;
        const double x11 = x10*x;
        const double x3 = x2*x;
        return x76 + x72 + x48 + x40 + x25 + x20 + x11 + x4 + x3;
    };

    const auto res = Rule::integrate(polynomial, limits);
    return close(res.val, 869368702.0/1258317060.0, 1.0e-13);
}

constexpr bool
gauss_kronrod_61_integrates_91st_degree_polynomial_exactly()
{
    using Rule = cubage::GaussKronrod<double, double, 61>;
    constexpr cubage::Interval<double> limits = {0.0, 1.0};

    auto polynomial = [](double x)
    {
        const double x2 = x*x;
        const double x4 = x2*x2;
        const double x8 = x4*x4;
        const double x16 = x8*x8;
        const double x32 = x16*x16;
        const double x64 = x32*x32;
        const double x24 = x16*x8;
        const double x88 = x64*x24;
        const double x3 = x2*x;
        const double x91 = x88*x3;
        const double x35 = x32*x3;
        const double x19 = x16*x3;
        const double x11 = x8*x3;
        return x91 + x88 + x64 + x35 + x32 + x19 + x11 + x4 + x3;
    };

    const auto res = Rule::integrate(polynomial, limits);
    return close(res.val, 35771317.0/52689780.0, 1.0e-13);
}

static_assert(gauss_kronrod_15_integrates_22nd_degree_polynomial_exactly());
static_assert(gauss_kronrod_21_integrates_31st_degree_polynomial_exactly());
static_assert(gauss_kronrod_31_integrates_46th_degree_polynomial_exactly());
static_assert(gauss_kronrod_41_integrates_61st_degree_polynomial_exactly());
static_assert(gauss_kronrod_51_integrates_76th_degree_polynomial_exactly());
static_assert(gauss_kronrod_61_integrates_91st_degree_polynomial_exactly());

int main()
{
    
}