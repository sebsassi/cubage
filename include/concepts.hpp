#pragma once

#include <concepts>
#include "integral_result.hpp"

namespace cubage
{

template <typename F, typename DomainType, typename CodomainType>
concept MapsAs = requires (F f, DomainType x)
{
    { f(x) } -> std::same_as<CodomainType>; 
};
}