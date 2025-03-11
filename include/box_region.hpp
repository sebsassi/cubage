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
#include <ranges>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <compare>

#include "concepts.hpp"

namespace cubage
{

template <typename FieldType>
    requires ArrayLike<FieldType> && FloatingPointVectorOperable<FieldType>
struct Box
{
    using value_type = typename FieldType::value_type;
    FieldType xmin;
    FieldType xmax;

    [[nodiscard]] constexpr FieldType side_lengths() const noexcept
    {
        return xmax - xmin;
    }

    [[nodiscard]] constexpr value_type volume() const noexcept
    {
        FieldType lengths = side_lengths();
        return std::accumulate(
                lengths.begin(), lengths.end(), 1.0,
                std::multiplies<value_type>());
    }

    [[nodiscard]] constexpr FieldType center() const noexcept
    {
        return 0.5*(xmax + xmin);
    }

    [[nodiscard]] constexpr std::pair<FieldType, FieldType>
    subdivide(std::size_t subdiv_axis) const noexcept
    {
        FieldType xmax_first = xmax;
        FieldType xmin_second = xmin;

        const auto mid = 0.5*(xmin[subdiv_axis] + xmax[subdiv_axis]);
        xmax_first[subdiv_axis] = mid;
        xmin_second[subdiv_axis] = mid;

        return {xmax_first, xmin_second};
    }
};

template <typename FieldType>
concept BoxIntegratorSignature
= requires (typename FieldType::CodomainType (*f)(typename FieldType::DomainType), typename FieldType::Limits limits)
{
    { FieldType::integrate(f, limits) } -> std::same_as<std::pair<IntegralResult<typename FieldType::CodomainType>, std::size_t>>;
};
    

template <typename Domain>
    requires ArrayLike<Domain>
class SubdivisibleBox
{
public:
    using DomainType = Domain;
    using Limits = Box<DomainType>;

    constexpr SubdivisibleBox() = default;

    constexpr SubdivisibleBox(
        const DomainType& p_xmin, const DomainType& p_xmax):
        SubdivisibleBox(Limits{p_xmin, p_xmax}) {}

    explicit constexpr SubdivisibleBox(const Limits& p_limits):
        m_limits(p_limits)
    {
        const auto sides = m_limits.side_lengths();
        for (const auto side : sides)
            if (side <= 0)
                throw std::invalid_argument(
                        "invalid integration limits: max <= min");
    }

    [[nodiscard]] constexpr const Limits&
    limits() const noexcept { return m_limits; }

    [[nodiscard]] constexpr std::pair<SubdivisibleBox, SubdivisibleBox>
    subdivide() const noexcept
    {
        const auto& [xmax_first, xmin_second] = m_limits.subdivide(m_subdiv_axis);

        std::pair<SubdivisibleBox, SubdivisibleBox> boxes = {
            SubdivisibleBox(m_limits.xmin, xmax_first),
            SubdivisibleBox(xmin_second, m_limits.xmax)
        };

        return boxes;
    }

    template <typename Rule, typename FuncType>
        requires MapsAs<FuncType, DomainType, typename Rule::CodomainType>
        && BoxIntegratorSignature<Rule>
    constexpr const IntegralResult<typename Rule::CodomainType> integrate(FuncType f) noexcept
    {
        const auto& [res, axis] = Rule::integrate(f, m_limits);
        m_subdiv_axis = axis;
        return res;
    }

private:
    Limits m_limits;
    std::size_t m_subdiv_axis;
};

}