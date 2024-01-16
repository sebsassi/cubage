#include "gauss_kronrod.hpp"
#include "interval_region.hpp"
#include "genz_malik.hpp"
#include "box_region.hpp"
#include "cubage.hpp"

#include <iostream>
#include <cassert>


constexpr bool close(double a, double b, double tol)
{
    return std::fabs(a - b) < tol;
}

bool gauss_kronrod_integrates_1d_gaussian()
{
    using Rule = cubage::GaussKronrod<double, double, 15>;
    using Region = cubage::IntegrationInterval<Rule>;
    using Integrator = cubage::MultiIntegrator<Region>;
    constexpr double sigma = 0.01;
    auto function = [sigma](double x)
    {
        const double z = x/sigma;
        return std::exp(-0.5*z*z);
    };

    constexpr double abserr = 1.0e-13;
    constexpr double relerr = 0.0;
    const std::vector<Region::Limits> limits = {Region::Limits{-1.0, 1.0}};
    auto result = Integrator().integrate(function, limits, abserr, relerr);
    std::cout << result.val << '\n';
    std::cout << result.err << '\n';
    return close(result.val, sigma*std::sqrt(2.0*M_PI), abserr);
}

bool genz_malik_integrates_2d_gaussian()
{
    using Rule = cubage::GenzMalikD7<std::array<double, 2>, double>;
    using Region = cubage::IntegrationBox<Rule>;
    using Integrator = cubage::MultiIntegrator<Region>;
    constexpr double sigma = 0.01;
    auto function = [sigma](const std::array<double, 2>& x)
    {
        const auto z = (1.0/sigma)*x;
        const auto z2 = z*z;
        return std::exp(-0.5*(z2[0] + z2[1]));
    };

    constexpr double abserr = 1.0e-13;
    constexpr double relerr = 0.0;
    const std::vector<Region::Limits> limits = {
        Region::Limits{{-1.0, -1.0}, {1.0, 1.0}}
    };
    auto result = Integrator().integrate(function, limits, abserr, relerr);
    std::cout << result.val << '\n';
    std::cout << result.err << '\n';
    return close(result.val, sigma*sigma*2.0*M_PI, abserr);
}

int main()
{
    assert(gauss_kronrod_integrates_1d_gaussian());
    assert(genz_malik_integrates_2d_gaussian());
}