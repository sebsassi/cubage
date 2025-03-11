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
#include <iostream>
#include <cassert>

#include "array_arithmetic.hpp"
#include "hypercube_integrator.hpp"

constexpr bool close(double a, double b, double tol)
{
    return std::fabs(a - b) < tol;
}

bool gauss_kronrod_integrates_1d_gaussian()
{
    using Integrator = cubage::IntervalIntegrator<double, double>;
    constexpr double sigma = 0.01;
    auto function = [sigma](double x)
    {
        const double z = x/sigma;
        return std::exp(-0.5*z*z);
    };

    constexpr double abserr = 1.0e-13;
    constexpr double relerr = 0.0;
    Integrator::Limits limits = Integrator::Limits{-1.0, 1.0};
    const auto& [result, _] = Integrator().integrate(function, limits, abserr, relerr);
    std::cout << result.val << '\n';
    std::cout << result.err << '\n';
    return close(result.val, sigma*std::sqrt(2.0*M_PI), abserr);
}

bool genz_malik_integrates_2d_gaussian()
{
    using Integrator = cubage::HypercubeIntegrator<std::array<double, 2>, double>;
    constexpr double sigma = 0.01;
    auto function = [sigma](const std::array<double, 2>& x)
    {
        const auto z = (1.0/sigma)*x;
        const auto z2 = z*z;
        return std::exp(-0.5*(z2[0] + z2[1]));
    };

    constexpr double abserr = 1.0e-13;
    constexpr double relerr = 0.0;
    Integrator::Limits limits = Integrator::Limits{{-1.0, -1.0}, {1.0, 1.0}};
    const auto& [result, _] = Integrator().integrate(function, limits, abserr, relerr);
    std::cout << result.val << '\n';
    std::cout << result.err << '\n';
    return close(result.val, sigma*sigma*2.0*M_PI, abserr);
}

bool genz_malik_integrates_3d_gaussian()
{
    using Integrator = cubage::HypercubeIntegrator<std::array<double, 3>, double>;
    constexpr double sigma = 0.01;
    auto function = [sigma](const std::array<double, 3>& x)
    {
        const auto z = (1.0/sigma)*x;
        const auto z2 = z*z;
        return std::exp(-0.5*(z2[0] + z2[1] + z2[2]));
    };

    constexpr double abserr = 1.0e-13;
    constexpr double relerr = 0.0;
    const std::vector<Integrator::Limits> limits = {
        Integrator::Limits{{-1.0, -1.0, -1.0}, {1.0, 1.0, 1.0}}
    };
    const auto& [result, _] = Integrator().integrate(function, limits, abserr, relerr);
    std::cout << result.val << '\n';
    std::cout << result.err << '\n';
    return close(result.val, sigma*sigma*sigma*std::pow(2.0*M_PI, 1.5), abserr);
}


int main()
{
    assert(gauss_kronrod_integrates_1d_gaussian());
    assert(genz_malik_integrates_2d_gaussian());
    assert(genz_malik_integrates_3d_gaussian());
}