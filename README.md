# Cubage

Cubage is a header-only template library for adaptive numerical integration in multiple dimensions. Currently, integration over hyperrectangular regions is supported using the Genz-Malik quadrature for dimensions > 1, and Gauss-Kronrod quadrature for the one-dimensional case.

This library differs from most multidimensional integration libraries in that the dimension of the integration domain is a compile-time template parameter instead of a runtime parameter. The rationale for this decision is that in any typical intergation problem, the dimension is almost always known at compile-time. Therefore providing the dimension as a compile-time parameter simplifies the memory management and enables more optimization opportunities.

The template-based approach enables the integration of any function `func` with signature `CodomainType func(DomainType)`. Here `DomainType` and `CodomainType` are any types, which model finite-dimensional vectors. That is, they are floating point scalars, or array-like types with floating point elements, which implement addition, subtraction, and scalar multiplication.

## Installation

As a header-only library, the only necessary step for using this library is to copy the contents of the `include` directory into your project:
```bash
cd cubage
cp -R include/. <path to project>/cubage
```
However, a CMake-based installation is also available. This installs cubage as a package, which can be used by other CMake projects:
```bash
cd cubage
cmake -S . -B build
cmake --install build --prefix <install directory>
```

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
#include "cubage/array_arithmetic.hpp"
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

