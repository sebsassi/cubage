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

template <typename T>
    requires ArrayLike<T> && FloatingPointVectorOperable<T>
struct Box
{
    using value_type = typename T::value_type;
    T xmin;
    T xmax;

    [[nodiscard]] constexpr T side_lengths() const noexcept
    {
        return xmax - xmin;
    }

    [[nodiscard]] constexpr value_type volume() const noexcept
    {
        T lengths = side_lengths();
        return std::accumulate(
                lengths.begin(), lengths.end(), 1.0,
                std::multiplies<value_type>());
    }

    [[nodiscard]] constexpr T center() const noexcept
    {
        return 0.5*(xmax + xmin);
    }

    [[nodiscard]] constexpr std::pair<T, T>
    subdivide(std::size_t subdiv_axis) const noexcept
    {
        T xmax_first = xmax;
        T xmin_second = xmin;

        const auto mid = 0.5*(xmin[subdiv_axis] + xmax[subdiv_axis]);
        xmax_first[subdiv_axis] = mid;
        xmin_second[subdiv_axis] = mid;

        return {xmax_first, xmin_second};
    }
};

template <typename T>
concept BoxIntegratorSignature
= requires (typename T::CodomainType (*f)(typename T::DomainType), typename T::Limits limits)
{
    { T::integrate(f, limits) } -> std::same_as<std::pair<IntegralResult<typename T::CodomainType>, std::size_t>>;
};
    

template <typename Rule>
    requires ArrayLike<typename Rule::DomainType>
    && BoxIntegratorSignature<Rule>
class IntegrationBox
{
public:
    using RuleType = Rule;
    using DomainType = typename Rule::DomainType;
    using CodomainType = typename Rule::CodomainType;
    using Limits = Box<DomainType>;
    using Result = IntegralResult<CodomainType>;

    constexpr IntegrationBox() = default;

    constexpr IntegrationBox(
        const DomainType& p_xmin, const DomainType& p_xmax):
        IntegrationBox(Limits{p_xmin, p_xmax}) {}

    explicit constexpr IntegrationBox(const Limits& p_limits):
        m_limits(p_limits)
    {
        const auto sides = m_limits.side_lengths();
        for (const auto side : sides)
            if (side <= 0)
                throw std::invalid_argument(
                        "invalid integration limits: max <= min");
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] constexpr std::pair<IntegrationBox, IntegrationBox>
    subdivide(FuncType f) const noexcept
    {
        const auto& [xmax_first, xmin_second] = m_limits.subdivide(m_subdiv_axis);

        std::pair<IntegrationBox, IntegrationBox> boxes = {
            IntegrationBox(m_limits.xmin, xmax_first),
            IntegrationBox(xmin_second, m_limits.xmax)
        };
        boxes.first.integrate(f);
        boxes.second.integrate(f);

        return boxes;
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    constexpr const IntegralResult<CodomainType>& integrate(FuncType f) noexcept
    {
        const auto& [res, axis] = Rule::integrate(f, m_limits);
        m_result = res;
        m_subdiv_axis = axis;
        if constexpr (std::is_floating_point<CodomainType>::value)
            m_maxerr = res.err;
        else
            m_maxerr = *std::ranges::max_element(res.err);
        return m_result;
    }

    [[nodiscard]] constexpr const IntegralResult<CodomainType>&
    result() const noexcept
    {
        return m_result;
    }

    [[nodiscard]] constexpr double maxerr() const noexcept
    {
        return m_maxerr;
    }

    constexpr auto operator<=>(const IntegrationBox& b) const noexcept
    {
        return maxerr() <=> b.maxerr();
    }

    constexpr bool operator==(const IntegrationBox& b) const noexcept
    {
        return maxerr() == b.maxerr();
    }

private:
    Limits m_limits;
    std::size_t m_subdiv_axis;
    IntegralResult<CodomainType> m_result;
    double m_maxerr;
};

}