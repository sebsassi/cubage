#include <cmath>

#include "array_arithmetic.hpp"
#include "box_region.hpp"

constexpr bool close(double a, double b, double tol)
{
    return std::fabs(a - b) < tol;
}

constexpr bool volume_of_3d_unit_box_is_1()
{
    const std::array<double, 3> xmin = {0.0, 0.0, 0.0};
    const std::array<double, 3> xmax = {1.0, 1.0, 1.0};
    const double volume = cubage::Box<std::array<double, 3>>{xmin, xmax}.volume();
    return close(volume, 1.0, 1.0e-13);
}

static_assert(volume_of_3d_unit_box_is_1());

int main()
{

}