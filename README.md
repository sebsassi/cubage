# Cubage

Cubage is a header-only template library for adaptive numerical integration in multiple dimensions. Currently, integration over hyperrectangular regions is supported using the Genz-Malik quadrature for dimensions > 1, and Gauss-Kronrod quadrature for the one-dimensional case. The composable template-based implementation, however, enables extension to other region shapes, mainly simplices, and other rules, such as Monte-Carlo.

Unlike most implementations of multidimensional integration, which have the dimensionality of the integral as a dynamic runtime parameter, in this library it is a compile-time constant parameter. The rationale for this is that the dimension of the integral rarely needs to change at runtime, and making it a compile-time parameter simplifies memory management and enables more optimization opportunities.

That is not the full extent of the template-based approach. More concretely, this library enables integration of any function `f(DomainType) -> CodomainType`, where `DomainType` and `CodomainType` are types which implement a finite-dimensional vector space (i.e., array-like types that implements addition, subtraction, and scalar multiplication).

## Examples

Example of integrating a 1D Gaussian over the interval `[-1, 1]`:
```cpp
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

    double a = -1.0;
    double b = +1.0;
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
```
Here `function` may be a lambda, function pointer, or any object which has a `CodomainType operator()(DomainType x)` method. In the `limits` parameter, multiple intervals are also accepted, e.g., as a `std::vector<Limits>`.

Example of integrating a 2D Gaussian over the box `[-1, 1]^2`:
```cpp
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
```
Since the `Integrator` expects `DomainType` to support basic vector algebra operations, the `array_arithmetic.hpp` header is provided as a convenience with implementations of the relevant operations for `std::array`.

