#include "array_arithmetic.hpp"
#include "cubage/hypercube_integrator.hpp"

int main()
{
    constexpr std::size_t NDIM = 2;

    auto function = [](const std::array<double, NDIM>& x)
    {
        constexpr double sigma = 1.0;
        const auto z = (1.0/sigma)*x;
        const auto z2 = z*z;
        return std::exp(-0.5*(z2[0] + z2[1]));
    };

    using DomainType = std::array<double, NDIM>;
    using CodomainType = double;
    using Integrator = cubage::HypercubeIntegrator<DomainType, CodomainType>;
    using Limits = typename Integrator::Limits;

    std::array<double, NDIM> a = {-1.0, -1.0};
    std::array<double, NDIM> b = {+1.0, +1.0};
    const Limits limits = {a, b};

    Integrator integrator{};

    constexpr double abserr = 1.0e-7;
    constexpr double relerr = 0.0;
    constexpr std::size_t max_subdiv = 2000;
    const auto& [res, status] = integrator.integrate(
            function, limits, abserr, relerr, max_subdiv);
    
    if (status == cubage::Status::MAX_SUBDIV)
        std::cout << "Warning: reached maximum number of subdivisions\n";

    std::cout << "Value: " << res.val << '\n';
    std::cout << "Error: " << res.err << '\n';
}