#pragma once

#include "multi_integrator.hpp"
#include "box_region.hpp"
#include "genz_malik.hpp"
#include "gauss_kronrod.hpp"

namespace cubage
{
template <std::floating_point DomainType, typename CodomainType, std::size_t Degree = 15>
using IntervalIntegrator = MultiIntegrator<GaussKronrod<DomainType, CodomainType, Degree>>;

template <GenzMalikIntegrable DomainType, typename CodomainType>
using HypercubeIntegrator = MultiIntegrator<GenzMalikD7<DomainType, CodomainType>>;
}