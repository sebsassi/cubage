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

#include "array_arithmetic.hpp"
#include "genz_malik.hpp"


constexpr bool close(double a, double b, double tol)
{
    return std::fabs(a - b) < tol;
}

constexpr double upow(double x, std::size_t n)
{
    double res = 1.0;
    const std::array<double, 2> mult = {1.0, x};
    while (1)
    {
        res *= mult[n & 1];     // if n is odd, multiply result by x
        n >>= 1;                // Divide n by 2
        if (!n) break;          // Stop iteration when n reaches 0
        x *= x;                 // Square x
    }
    return res;
}

constexpr bool
constant_unity_function_in_3d_null_box_integrates_to_zero()
{
    using Rule = cubage::GenzMalikD7<std::array<double, 3>, double>;
    constexpr cubage::Box<std::array<double, 3>> limits = {
        std::array<double, 3>{0.0, 0.0, 0.0},
        std::array<double, 3>{0.0, 0.0, 0.0}
    };
    const auto& [res, axis] = Rule::integrate([]([[maybe_unused]] std::array<double, 3> x){ return 1.0; }, limits);
    return close(res.val, 0.0, 1.0e-13);
}

constexpr bool
constant_zero_function_in_3d_unit_box_integrates_to_zero()
{
    using Rule = cubage::GenzMalikD7<std::array<double, 3>, double>;
    constexpr cubage::Box<std::array<double, 3>> limits = {
        std::array<double, 3>{0.0, 0.0, 0.0},
        std::array<double, 3>{1.0, 1.0, 1.0}
    };

    auto zero_function = []([[maybe_unused]] const std::array<double, 3>& x)
    {
        return 0.0;
    };
    const auto& [res, axis] = Rule::integrate(zero_function, limits);
    return close(res.val, 0.0, 1.0e-13);
}

constexpr bool
constant_unity_function_in_3d_unit_box_integrates_to_unity()
{
    using Rule = cubage::GenzMalikD7<std::array<double, 3>, double>;
    constexpr cubage::Box<std::array<double, 3>> limits = {
        std::array<double, 3>{0.0, 0.0, 0.0},
        std::array<double, 3>{1.0, 1.0, 1.0}
    };

    auto unity_function = []([[maybe_unused]] const std::array<double, 3>& x){
        return 1.0;
    };
    const auto& [res, axis] = Rule::integrate(unity_function, limits);
    return close(res.val, 1.0, 1.0e-13);
}

constexpr bool 
linear_function_with_unit_coeffs_in_3d_unit_box_integrates_to_3_halves()
{
    using Rule = cubage::GenzMalikD7<std::array<double, 3>, double>;
    constexpr cubage::Box<std::array<double, 3>> limits = {
        std::array<double, 3>{0.0, 0.0, 0.0},
        std::array<double, 3>{1.0, 1.0, 1.0}
    };

    auto linear_function = [](std::array<double, 3> x)
    {
        return x[0] + x[1] + x[2];
    };
    const auto& [res, axis] = Rule::integrate(linear_function, limits);
    return close(res.val, 3.0/2.0, 1.0e-13);
}

constexpr bool 
radius_squared_in_3d_unit_box_integrates_to_unity()
{
    using Rule = cubage::GenzMalikD7<std::array<double, 3>, double>;
    constexpr cubage::Box<std::array<double, 3>> limits = {
        std::array<double, 3>{0.0, 0.0, 0.0},
        std::array<double, 3>{1.0, 1.0, 1.0}
    };

    auto radius_squared = [](std::array<double, 3> x)
    {
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2]; 
    };
    const auto& [res, axis] = Rule::integrate(radius_squared, limits);
    return close(res.val, 1.0, 1.0e-13);
}

constexpr bool
seventh_degree_polynomial_integrates_exactly()
{
    using Rule = cubage::GenzMalikD7<std::array<double, 3>, double>;
    constexpr cubage::Box<std::array<double, 3>> limits = {
        std::array<double, 3>{0.0, 0.0, 0.0},
        std::array<double, 3>{1.0, 1.0, 1.0}
    };
    auto polynomial = [](std::array<double, 3> x)
    {
        return x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]
            + x[0]*x[0]*x[0]*x[1]*x[1]*x[1]*x[2]
            + x[0]*x[0]*x[2]*x[2]
            + x[0]*x[1]*x[2];
    };
    const auto& [res, axis] = Rule::integrate(polynomial, limits);
    return close(res.val, 113.0/288.0, 1.0e-13);
}

constexpr bool
error_of_fift_degree_polynomial_integral_is_zero()
{
    using Rule = cubage::GenzMalikD7<std::array<double, 3>, double>;
    constexpr cubage::Box<std::array<double, 3>> limits = {
        std::array<double, 3>{0.0, 0.0, 0.0},
        std::array<double, 3>{1.0, 1.0, 1.0}
    };
    auto polynomial = [](std::array<double, 3> x)
    {
        return x[0]*x[0]*x[0]*x[0]*x[0]
            + x[0]*x[0]*x[1]*x[1]*x[2]
            + x[0]*x[0]*x[2]*x[2]
            + x[0]*x[1]*x[2];
    };
    const auto& [res, axis] = Rule::integrate(polynomial, limits);
    return close(res.err, 0.0, 1.0e-13);
}

constexpr bool
subdiv_axis_is_in_nonconst_direction()
{
    using Rule = cubage::GenzMalikD7<std::array<double, 3>, double>;
    constexpr cubage::Box<std::array<double, 3>> limits = {
        std::array<double, 3>{0.0, 0.0, 0.0},
        std::array<double, 3>{1.0, 1.0, 1.0}
    };
    auto polynomial = [](std::array<double, 3> x)
    {
        return std::cos(x[2]);
    };
    const auto& [res, axis] = Rule::integrate(polynomial, limits);
    return axis == 2;
}

static_assert(constant_unity_function_in_3d_null_box_integrates_to_zero());
static_assert(constant_zero_function_in_3d_unit_box_integrates_to_zero());
static_assert(constant_unity_function_in_3d_unit_box_integrates_to_unity());
static_assert(linear_function_with_unit_coeffs_in_3d_unit_box_integrates_to_3_halves());
static_assert(radius_squared_in_3d_unit_box_integrates_to_unity());
static_assert(seventh_degree_polynomial_integrates_exactly());
static_assert(error_of_fift_degree_polynomial_integral_is_zero());
static_assert(subdiv_axis_is_in_nonconst_direction());

int main()
{
    
}