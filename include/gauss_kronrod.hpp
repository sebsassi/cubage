#pragma once

#include <cmath>
#include <ranges>
#include <algorithm>

#include "integral_result.hpp"
#include "gauss_kronrod_data.hpp"
#include "interval_region.hpp"

namespace cubage
{

template <std::size_t Degree>
concept GKIsImplemented = requires
{ 
    GK<Degree>::degree == Degree;
    { GK<Degree>::gauss_points() }
        -> std::same_as<std::array<double, (Degree - 1)/4>>;
    { GK<Degree>::kronrod_points() }
        -> std::same_as<std::array<double, (Degree + 1)/4>>;
    { GK<Degree>::gauss_weights() }
        -> std::same_as<std::array<double, (Degree + 1)/4>>;
    
    { GK<Degree>::kronrod_weights() }
        -> std::same_as<
                std::tuple<
                    double,
                    std::array<double, (Degree - 1)/4>,
                    std::array<double, (Degree + 1)/4>
                >
            >;
};

template <std::floating_point DomainTypeParam, typename CodomainTypeParam, std::size_t Degree>
    requires GKIsImplemented<Degree>
    && (std::floating_point<CodomainTypeParam>
        || (FloatingPointVectorOperable<CodomainTypeParam>
            && ArrayLike<CodomainTypeParam>))
struct GaussKronrod
{
    using DomainType = DomainTypeParam;
    using CodomainType = CodomainTypeParam;
    using ReturnType = IntegralResult<CodomainType>;
    using RuleData = GK<Degree>;
    using Limits = Interval<DomainType>;
    using RegionType = SubdivisibleInterval<DomainType>;

    template <typename FuncType>
        requires MapsAs<FuncType, DomainType, CodomainType>
    [[nodiscard]] static constexpr
    ReturnType integrate(FuncType f, const Limits& limits) noexcept
    {
        constexpr auto gauss_points = RuleData::gauss_points();
        constexpr auto kronrod_points = RuleData::kronrod_points();

        constexpr auto gauss_weights = RuleData::gauss_weights();
        constexpr auto kronrod_weights = RuleData::kronrod_weights();
        constexpr auto center_weight = std::get<0>(kronrod_weights);
        constexpr auto kronrod_weights_g = std::get<1>(kronrod_weights);
        constexpr auto kronrod_weights_k = std::get<2>(kronrod_weights);

        const DomainType center = limits.center();
        const DomainType half_length = 0.5*limits.length();

        const CodomainType central_value = f(center);

        std::array<std::pair<CodomainType, CodomainType>, gauss_points.size()> gauss_values{};
        for (size_t i = 0; i < gauss_points.size(); ++i)
        {
            const DomainType disp = half_length*gauss_points[i];
            gauss_values[i].first = f(center + disp);
            gauss_values[i].second = f(center - disp);
        }

        std::array<std::pair<CodomainType, CodomainType>, kronrod_points.size()> kronrod_values{};
        for (size_t i = 0; i < kronrod_points.size(); ++i)
        {
            const DomainType disp = half_length*kronrod_points[i];
            kronrod_values[i].first = f(center + disp);
            kronrod_values[i].second = f(center - disp);
        }

        CodomainType val = central_value*center_weight;
        CodomainType val_gauss{};

        for (size_t i = 0; i < gauss_points.size(); ++i)
        {
            const auto& [first, second] = gauss_values[i];
            const CodomainType sum = first + second;
            val += sum*kronrod_weights_g[i];
            val_gauss += sum*gauss_weights[i];
        }

        if constexpr (((Degree - 1)/2) & 1)
            val_gauss += central_value*gauss_weights.back();

        for (size_t i = 0; i < kronrod_points.size(); ++i)
        {
            const auto& [first, second] = kronrod_values[i];
            val += (first + second)*kronrod_weights_k[i];
        }

        const CodomainType err_null_1 = vfabs(val - val_gauss);
        
        const CodomainType favg = 0.5*val;
        CodomainType err_null_0 = vfabs(central_value - favg)*center_weight;
        for (size_t i = 0; i < gauss_points.size(); ++i)
        {
            const auto& [first, second] = gauss_values[i];
            err_null_0 += (vfabs(first - favg) + vfabs(second - favg))*kronrod_weights_g[i];
        }
        for (size_t i = 0; i < kronrod_points.size(); ++i)
        {
            const auto& [first, second] = kronrod_values[i];
            err_null_0 += (vfabs(first - favg) + vfabs(second - favg))*kronrod_weights_k[i];
        }

        CodomainType err = err_null_0;
        if constexpr (std::is_floating_point<CodomainType>::value)
            err *= half_length*berntsen_espelid_estimate(
                    err_null_1, err_null_0);
        else
        {
            for (std::size_t i = 0; i < std::tuple_size<CodomainType>::value; ++i)
                err[i] *= half_length*berntsen_espelid_estimate(
                        err_null_1[i], err_null_0[i]);
        }

        return IntegralResult<CodomainType>{half_length*val, err};
    }

    [[nodiscard]] static constexpr std::size_t num_points() noexcept
    {
        return Degree;
    }

private:
    [[nodiscard]] static constexpr CodomainType
    vfabs(const CodomainType& x)
    {
        if constexpr (std::is_floating_point<CodomainType>::value)
            return std::fabs(x);
        else
#if (__GNUC__ > 13)
            return x
                | std::views::transform(std::fabs)
                | std::ranges::to<CodomainType>();
#else
        {
            CodomainType res{};
            std::views::transform(x, res, std::fabs);
            return res;
        }
#endif
    }

    /*
        Error estimation method of Berntsen and Espelid

            Jarle Berntsen, Terje O. Espelid, "Error Estimation in Automatic Quadrature Routines", ACM Trans. Math. Softw. 17:233-252, 1991
    */
    [[nodiscard]] static constexpr double berntsen_espelid_estimate(
        double err_null_1, double err_null_0)
    {
        const double ratio = err_null_1/err_null_0;
        if (200.0*ratio <= 1.0)
        {
            const double factor = std::sqrt(200.0*ratio);
            return factor*factor*factor;
        }
        else
            return 1.0;
    }
};

}