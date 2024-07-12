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

template <typename T>
concept BiSubdivisible = requires (T x, typename T::CodomainType (*f)(typename T::DomainType))
{
    { x.subdivide(f) } -> std::same_as<std::pair<T, T>>;
};

template <typename T>
concept Limited = requires { typename T::Limits; };

template <typename T>
concept ResultStoring = requires (T x)
{
    { x.result() } -> std::same_as<const IntegralResult<typename T::CodomainType>&>;
};

template <typename T>
concept WeaklyOrdered = requires (T x, T y)
{
    x < y;
};

template <typename T>
concept Integrating =
requires (T x, typename T::CodomainType (*f)(typename T::DomainType))
{
    { x.integrate(f) } -> std::same_as<const IntegralResult<typename T::CodomainType>&>;
};

template <typename T>
concept SubdivisionIntegrable
    = WeaklyOrdered<T> && Limited<T> && Integrating<T> && BiSubdivisible<T>
    && ResultStoring<T>;

template <SubdivisionIntegrable RegionType, typename NormType = NormIndividual>
class MultiIntegrator
{
public:
    using Region = RegionType;
    using Limits = typename Region::Limits;
    using CodomainType = typename Region::CodomainType;
    using DomainType = typename Region::DomainType;
    using Result = IntegralResult<CodomainType>;

    MultiIntegrator() = default;

    template <typename FuncType, typename LimitsType>
        requires MapsAs<FuncType, DomainType, CodomainType>
        &&  (std::same_as<std::remove_cvref_t<LimitsType>, Limits> || std::convertible_to<LimitsType, std::span<const Limits>>)
    [[nodiscard]] Result integrate(
            FuncType f, LimitsType&& integration_domain,
            double abserr, double relerr)
    {
        if constexpr (std::same_as<std::remove_cvref_t<LimitsType>, Limits>)
            m_region_eval_count = 1;
        else
            m_region_eval_count = integration_domain.size();
        generate(integration_domain);
        Result res = initialize(f);

        while (!has_converged(res, abserr, relerr))
            subdivide_top_region(f, res);
        
        // resum to minimize spooky floating point error accumulation
        res = Result{};
        for (const auto& region : m_region_heap)
            res += region.result();
        return res;
    }

    [[nodiscard]] std::size_t func_eval_count() const noexcept
    {
        return m_region_eval_count*Region::RuleType::points_count();
    }
    
    [[nodiscard]] std::size_t region_eval_count() const noexcept
    {
        return m_region_eval_count;
    }

    [[nodiscard]] std::size_t region_count() const noexcept
    {
        return m_region_heap.size();
    }

    [[nodiscard]] std::size_t capacity() const noexcept
    {
        return m_region_heap.capacity();
    }

private:
    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] inline Result initialize(FuncType f)
    {
        Result res{};
        for (auto& region : m_region_heap)
            res += region.integrate(f);
        std::ranges::make_heap(m_region_heap);

        return res;
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    inline void subdivide_top_region(FuncType f, Result& res)
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
        const Result& res, double abserr, double relerr) const noexcept
    {
        if constexpr (std::floating_point<CodomainType>)
            return res.err <= abserr || res.err <= res.val*relerr;
        else
        {
            if constexpr (std::is_same_v<NormType, NormIndividual>)
            {
#if (__GNUC__ > 12)
                for (const auto& [val, err] : std::ranges::views::zip(res.val, res.err))
                    if (err > abserr && err > std::fabs(val)*relerr) return false;
                return true;
#else
                for (std::size_t i = 0; i < res.ndim(); ++i)
                {
                    if (res.err[i] > abserr
                            && res.err[i] > std::fabs(res.val[i])*relerr)
                        return false;
                }
                return true;
#endif
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

    void generate(std::span<const Limits> limits)
    {
        m_region_heap.clear();
        m_region_heap.reserve(limits.size());
        for (const auto& limit : limits)
            m_region_heap.emplace_back(limit);
    }

    void generate(const Limits& limits)
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