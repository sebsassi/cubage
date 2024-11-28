#pragma once

#include <vector>
#include <algorithm>
#include <ranges>
#include <utility>
#include <cmath>
#include <type_traits>
#include <concepts>
#include <span>

#include <iostream>

#include "integral_result.hpp"
#include "concepts.hpp"

namespace cubage
{

struct NormIndividual {};

template <typename FieldType>
concept BiSubdivisible = requires (FieldType x, typename FieldType::CodomainType (*f)(typename FieldType::DomainType))
{
    { x.subdivide(f) } -> std::same_as<std::pair<FieldType, FieldType>>;
};

template <typename FieldType>
concept Limited = requires { typename FieldType::Limits; };

template <typename FieldType>
concept ResultStoring = requires (FieldType x)
{
    { x.result() } -> std::same_as<const IntegralResult<typename FieldType::CodomainType>&>;
};

template <typename FieldType>
concept WeaklyOrdered = requires (FieldType x, FieldType y)
{
    x < y;
};

template <typename FieldType>
concept Integrating =
requires (FieldType x, typename FieldType::CodomainType (*f)(typename FieldType::DomainType))
{
    { x.integrate(f) } -> std::same_as<const IntegralResult<typename FieldType::CodomainType>&>;
};

template <typename Rule>
class IntegrationRegion
{
public:
    using RuleType = Rule;
    using RegionType = Rule::RegionType;
    using DomainType = typename Rule::DomainType;
    using CodomainType = typename Rule::CodomainType;
    using Limits = RuleType::Limits;
    using Result = IntegralResult<CodomainType>;

    explicit constexpr IntegrationRegion(const Limits& p_limits):
        m_region(p_limits) {}

    explicit constexpr IntegrationRegion(const RegionType& p_region):
        m_region(p_region) {}

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] constexpr std::pair<IntegrationRegion, IntegrationRegion>
    subdivide(FuncType f) const noexcept
    {
        const auto& [left, right] = m_region.subdivide();

        std::pair<IntegrationRegion, IntegrationRegion> regions = {
            IntegrationRegion(left), IntegrationRegion(right)
        };
        regions.first.integrate(f);
        regions.second.integrate(f);

        return regions;
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    constexpr const IntegralResult<CodomainType>& integrate(FuncType f) noexcept
    {
        m_result = m_region.template integrate<RuleType>(f);
        if constexpr (std::is_floating_point<CodomainType>::value)
            m_maxerr = m_result.err;
        else
            m_maxerr = *std::ranges::max_element(m_result.err);
        return m_result;
    }

    [[nodiscard]] constexpr const IntegralResult<CodomainType>&
    result() const noexcept { return m_result; }

    [[nodiscard]] constexpr double
    maxerr() const noexcept { return m_maxerr; }

    constexpr auto operator<=>(const IntegrationRegion& b) const noexcept
    {
        return maxerr() <=> b.maxerr();
    }
    
    constexpr bool operator==(const IntegrationRegion& b) const noexcept
    {
        return maxerr() == b.maxerr();
    }

    [[nodiscard]] constexpr const Limits&
    limits() const noexcept { return m_region.limits(); }

private:
    RegionType m_region;
    IntegralResult<CodomainType> m_result;
    double m_maxerr;
};

template <typename FieldType>
concept SubdivisionIntegrable
    = WeaklyOrdered<FieldType> && Limited<FieldType> && Integrating<FieldType> && BiSubdivisible<FieldType>
    && ResultStoring<FieldType>;

template <typename T, typename ValueType>
concept SizedRangeOf = std::ranges::sized_range<T>
    && std::same_as<std::ranges::range_value_t<T>, ValueType>;

template <typename T, typename ValueType>
concept ValueOrSizedRangeOf
    = std::same_as<T, ValueType> || SizedRangeOf<T, ValueType>;

template <typename RuleType, typename NormType = NormIndividual>
class MultiIntegrator
{
public:
    using RegionType = IntegrationRegion<RuleType>;
    using Limits = typename RegionType::Limits;
    using CodomainType = typename RegionType::CodomainType;
    using DomainType = typename RegionType::DomainType;
    using ResultType = IntegralResult<CodomainType>;

    MultiIntegrator() = default;

    template <typename FuncType, typename LimitsType>
        requires MapsAs<FuncType, DomainType, CodomainType>
            && ValueOrSizedRangeOf<LimitsType, Limits>
    [[nodiscard]] Result<ResultType, Status> integrate(
            FuncType f, LimitsType&& integration_domain,
            double abserr, double relerr,
            std::size_t max_subdiv = std::numeric_limits<std::size_t>::max())
    {
        if constexpr (SizedRangeOf<LimitsType, Limits>)
            m_region_eval_count = std::ranges::size(integration_domain);
        else
            m_region_eval_count = 1;
        generate_region_heap(integration_domain);
        ResultType res = integrate_initial_regions(f);

        while (!has_converged(res, abserr, relerr) && m_region_heap.size() < max_subdiv)
            subdivide_top_region(f, res);
        
        // resum to minimize spooky floating point error accumulation
        res = ResultType{};
        for (const auto& region : m_region_heap)
            res += region.result();
        
        Status status = (m_region_heap.size() >= max_subdiv) ? 
            Status::MAX_SUBDIV : Status::SUCCESS;
        return {res, status};
    }

    [[nodiscard]] std::size_t func_eval_count() const noexcept
    {
        return m_region_eval_count*RuleType::points_count();
    }

    [[nodiscard]] std::size_t region_eval_count() const noexcept
    {
        return m_region_eval_count;
    }
    [[nodiscard]] std::size_t region_count() const noexcept
    {
        return m_region_heap.size();
    }

    [[nodiscard]] std::span<const RegionType> regions() const noexcept
    {
        return std::span(m_region_heap);
    }

    [[nodiscard]] std::size_t capacity() const noexcept
    {
        return m_region_heap.capacity();
    }

private:
    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] inline ResultType integrate_initial_regions(FuncType f)
    {
        ResultType res{};
        for (auto& region : m_region_heap)
            res += region.integrate(f);
        std::ranges::make_heap(m_region_heap);

        return res;
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    inline void subdivide_top_region(FuncType f, ResultType& res)
    {
        m_region_eval_count += 2;
        const RegionType top_region = pop_top_region();

        const std::pair<RegionType, RegionType> new_regions
            = top_region.subdivide(f);
        
        res += new_regions.first.result() + new_regions.second.result()
            - top_region.result();

        push_to_heap(new_regions.first);
        push_to_heap(new_regions.second);
    }

    [[nodiscard]] inline bool has_converged(
        const ResultType& res, double abserr, double relerr) const noexcept
    {
        if constexpr (std::floating_point<CodomainType>)
            return res.err <= abserr || res.err <= res.val*relerr;
        else
        {
            if constexpr (std::is_same_v<NormType, NormIndividual>)
            {
                for (std::size_t i = 0; i < res.ndim(); ++i)
                {
                    if (res.err[i] > abserr
                            && res.err[i] > std::fabs(res.val[i])*relerr)
                        return false;
                }
                return true;
            }
            else
            {
                const double norm_val = NormType::norm(res.val);
                const double norm_err = NormType::norm(res.err);

                return norm_err <= abserr || norm_err <= norm_val*relerr;
            }
        }
    }

    inline void push_to_heap(const RegionType& region)
    {
        m_region_heap.push_back(region);
        std::ranges::push_heap(m_region_heap);
    }

    [[nodiscard]] inline RegionType pop_top_region()
    {
        std::ranges::pop_heap(m_region_heap);
        RegionType top_region = m_region_heap.back();
        m_region_heap.pop_back();
        return top_region;
    }

    template <typename LimitsRange>
        requires SizedRangeOf<LimitsRange, Limits>
    void generate_region_heap(LimitsRange&& limits)
    {
        m_region_heap.clear();
        m_region_heap.reserve(std::ranges::size(limits));
        for (const auto& limit : limits)
            m_region_heap.emplace_back(limit);
    }

    void generate_region_heap(const Limits& limits)
    {
        m_region_heap.clear();
        m_region_heap.reserve(1);
        m_region_heap.emplace_back(limits);
    }

private:
    std::vector<RegionType> m_region_heap;
    std::size_t m_region_eval_count;
};

}