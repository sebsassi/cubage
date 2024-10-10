#pragma once

#include <array>
#include <ranges>
#include <algorithm>
#include <iterator>
#include <numeric>

#include "concepts.hpp"

namespace cubage
{

template <std::floating_point FieldType>
struct Interval
{
    FieldType xmin;
    FieldType xmax;

    [[nodiscard]] constexpr inline FieldType length() const noexcept
    {
        return xmax - xmin;
    }

    [[nodiscard]] constexpr inline FieldType center() const noexcept
    {
        return 0.5*(xmax + xmin);
    }
};


template <typename FieldType>
concept IntervalIntegratorSignature
= requires (typename FieldType::CodomainType (*f)(typename FieldType::DomainType), typename FieldType::Limits limits)
{
    { FieldType::integrate(f, limits) } -> std::same_as<IntegralResult<typename FieldType::CodomainType>>;
};

template <typename Domain>
    requires std::floating_point<Domain>
class SubdivisibleInterval
{
public:
    using DomainType = Domain;
    using Limits = Interval<DomainType>;

    constexpr SubdivisibleInterval() = default;

    constexpr SubdivisibleInterval(
        const DomainType& p_xmin, const DomainType& p_xmax):
        SubdivisibleInterval(Limits{p_xmin, p_xmax}) {}

    explicit constexpr SubdivisibleInterval(const Limits& p_limits):
        m_limits(p_limits)
    {
        if (m_limits.length() <= 0)
            throw std::invalid_argument(
                    "invalid integration limits: max <= min");
    }

    [[nodiscard]] constexpr const Limits&
    limits() const noexcept { return m_limits; }

    [[nodiscard]] constexpr std::pair<SubdivisibleInterval, SubdivisibleInterval>
    subdivide() const noexcept
    {
        const DomainType mid = m_limits.center();

        std::pair<SubdivisibleInterval, SubdivisibleInterval> intervals = {
            SubdivisibleInterval(m_limits.xmin, mid),
            SubdivisibleInterval(mid, m_limits.xmax)
        };

        return intervals;
    }

    template <typename Rule, typename FuncType>
        requires MapsAs<FuncType, DomainType, typename Rule::CodomainType> && IntervalIntegratorSignature<Rule>
    constexpr const IntegralResult<typename Rule::CodomainType> integrate(FuncType f) noexcept
    {
        return Rule::integrate(f, m_limits);
    }

private:
    Limits m_limits;
};

}