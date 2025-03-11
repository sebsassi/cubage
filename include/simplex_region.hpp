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
#pragma once

#include <array>

#include "concepts.hpp"

template <typename T>
struct Simplex
{
    std::array<T, std::tuple_size<T>::value + 1> vertices;

    [[nodiscard]] std::pair<Simplex, Simplex> subdivide_longest_edge()
    {
        std::pair<std::size_t, std::size_t> inds{};
        double max_length_sqr = 0.0;
        for (std::size_t i = 0; i < std::tuple_size<T>::value; ++i)
        {
            for (std::size_t j = i + 1; j < std::tuple_size<T>::value; ++j)
            {
                const T disp = vertices[j] - vertices[i];
                const T disp2 = disp*disp;
                const double length_sqr = std::accumulate(
                        disp2.begin(), disp2.end(), 0.0);
                if (length_sqr > max_length_sqr)
                {
                    max_length_sqr = length_sqr;
                    inds = {i, j};
                }
            }
        }

        const T center = 0.5*(vertices[inds.first] + vertices[inds.second]);

        std::pair<Simplex, Simplex> res = {*this, *this};
        res.first[inds.first] = center;
        res.second[inds.second] = center;

        return res;
    }
};

template <typename Rule>
class IntegrationSimplex
{
public:
    using DomainType = typename Rule::DomainType;
    using CodomainType = typename Rule::CodomainType;
    using Limits = Simplex<DomainType>;
    using Result = IntegralResult<CodomainType>;

    IntegrationSimplex() = default;

    IntegrationSimplex(const Limits& p_limits)
    : limits(p_limits), result_(), maxerr_() {}

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] constexpr std::pair<IntegrationSimplex, IntegrationSimplex>
    subdivide(FuncType f) const
    {
        const std::pair<Simplex, Simplex> verts
            = limits.subdivide_longest_edge();

        std::pair<IntegrationSimplex, IntegrationSimplex> simplices = {
            IntegrationSimplex(verts.first),
            IntegrationSimplex(verts.second)
        };
        simplices.first.integrate(f);
        simplices.second.integrate(f);

        return simplices;
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    const IntegralResult<CodomainType>& integrate(FuncType f)
    {
        result_ = Rule::integrate(f, limits);
        if constexpr (std::is_floating_point<CodomainType>::value)
            maxerr_ = res.err;
        else
            maxerr_ = *std::ranges::max_element(res.err);
        return result_;
    }

private:
    Limits limits;
    IntegralResult<CodomainType> result_;
    double maxerr_;
};