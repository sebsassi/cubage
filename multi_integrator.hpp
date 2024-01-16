#pragma once

#include <vector>
#include <algorithm>
#include <ranges>
#include <utility>
#include <cmath>
#include <type_traits>
#include <concepts>

#include <iostream>

#include "array_arithmetic.hpp"
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
    using CodomainType = typename RegionType::CodomainType;
    using DomainType = typename RegionType::DomainType;
    using Result = IntegralResult<CodomainType>;

    MultiIntegrator() = default;

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] Result integrate(
            FuncType f,
            const std::vector<typename RegionType::Limits>& integration_domain,
            double abserr, double relerr)
    {
        generate(region_heap, integration_domain);
        Result res = initialize(f);

        while (!has_converged(res, abserr, relerr))
        {
            subdivide_top_region(f, res);
        }
        
        // resum to minimize spooky floating point error accumulation
        res = Result{};
        for (const auto& region : region_heap)
            res += region.result();
        return res;
    }

private:
    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] inline Result initialize(FuncType f)
    {
        Result res{};
        for (auto& region : region_heap)
            res += region.integrate(f);
        std::ranges::make_heap(region_heap);

        return res;
    }

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    inline void subdivide_top_region(FuncType f, Result& res)
    {
        const RegionType top_region = pop_top_region();

        const std::pair<RegionType, RegionType> new_regions
            = top_region.subdivide(f);
        
        res += new_regions.first.result() + new_regions.second.result()
            - top_region.result();

        push_to_heap(new_regions.first);
        push_to_heap(new_regions.second);
    }

    [[nodiscard]] inline bool has_converged(
        const Result& res, double abserr, double relerr) const
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
        region_heap.push_back(region);
        std::ranges::push_heap(region_heap);
    }

    [[nodiscard]] inline RegionType pop_top_region()
    {
        std::ranges::pop_heap(region_heap);
        RegionType top_region = region_heap.back();
        region_heap.pop_back();
        return top_region;
    }

private:
    std::vector<RegionType> region_heap;
};

template <SubdivisionIntegrable RegionType, typename LimitType>
    requires std::is_same_v<typename RegionType::Limits, LimitType>
static void generate(
    std::vector<RegionType>& regions, const std::vector<LimitType>& limits)
{
    regions.resize(limits.size());
#if (__GNUC__ > 12)
    for (auto& [region, limit] : std::ranges::views::zip(regions, limits))
        region = RegionType(limit);
#else
    for (std::size_t i = 0; i < limits.size(); ++i)
        regions[i] = RegionType(limits[i]);
#endif
}

}