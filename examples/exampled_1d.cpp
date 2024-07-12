#include "cubage/hypercube_integrator.hpp"

int main()
{
    auto function = [](double x)
    {
        constexpr double sigma = 1.0;
        const auto z = (1.0/sigma)*x;
        const auto z2 = z*z;
        return std::exp(-0.5*z2);
    };

    using DomainType = double;
    using CodomainType = double;
    using Integrator = cubage::IntervalIntegrator<DomainType, CodomainType>;
    using Limits = typename Integrator::Limits;
    using Result = typename Integrator::Result;

    double a = -1.0;
    double b = +1.0;
    const Limits limits = {a, b};

    Integrator integrator{};

    constexpr double abserr = 1.0e-7;
    constexpr double relerr = 0.0;
    const Result res = integrator.integrate(
            function, limits, abserr, relerr);
    std::cout << "Value: " << res.val << '\n';
    std::cout << "Error: " << res.err << '\n';
}