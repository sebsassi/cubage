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