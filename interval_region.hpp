#pragma once

#include <array>
#include <ranges>
#include <algorithm>
#include <iterator>
#include <numeric>

#include "concepts.hpp"

namespace cubage
{

template <std::floating_point T>
struct Interval
{
    T xmin;
    T xmax;

    [[nodiscard]] constexpr inline T length() const noexcept
    {
        return xmax - xmin;
    }

    [[nodiscard]] constexpr inline T center() const noexcept
    {
        return 0.5*(xmax + xmin);
    }
};


template <typename T>
concept IntervalIntegratorSignature
= requires (typename T::CodomainType (*f)(typename T::DomainType), typename T::Limits limits)
{
    { T::integrate(f, limits) } -> std::same_as<IntegralResult<typename T::CodomainType>>;
};

template <typename Rule>
    requires std::floating_point<typename Rule::DomainType>
    && IntervalIntegratorSignature<Rule>
class IntegrationInterval
{
public:
    using RuleType = Rule;
    using DomainType = typename Rule::DomainType;
    using CodomainType = typename Rule::CodomainType;
    using Limits = Interval<DomainType>;
    using Result = IntegralResult<CodomainType>;

    constexpr IntegrationInterval() = default;

    constexpr IntegrationInterval(const DomainType& p_xmin, const DomainType& p_xmax)
    : limits{p_xmin, p_xmax}, m_result(), m_maxerr() {}

    explicit constexpr IntegrationInterval(const Limits& p_limits)
    : limits(p_limits), m_result(), m_maxerr() {}

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] constexpr std::pair<IntegrationInterval, IntegrationInterval>
    subdivide(FuncType f) const noexcept
    {
        const DomainType mid = limits.center();

        std::pair<IntegrationInterval, IntegrationInterval> intervals = {
            IntegrationInterval(limits.xmin, mid),
            IntegrationInterval(mid, limits.xmax)
        };
        intervals.first.integrate(f);
        intervals.second.integrate(f);

        return intervals;
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    constexpr const IntegralResult<CodomainType>& integrate(FuncType f) noexcept
    {
        m_result = Rule::integrate(f, limits);
        if constexpr (std::is_floating_point<CodomainType>::value)
            m_maxerr = m_result.err;
        else
            m_maxerr = *std::ranges::max_element(m_result.err);
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

    constexpr auto operator<=>(const IntegrationInterval& b) const noexcept
    {
        return maxerr() <=> b.maxerr();
    }
    
    constexpr bool operator==(const IntegrationInterval& b) const noexcept
    {
        return maxerr() == b.maxerr();
    }

private:
    Limits limits;
    IntegralResult<CodomainType> m_result;
    double m_maxerr;
};

}