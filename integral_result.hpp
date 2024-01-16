#pragma once

#include <concepts>
#include <utility>

namespace cubage
{

template <typename T>
concept VectorValued = requires (T a, T b)
{
    a += b;
    a -= b;
    a + b;
    a - b;
}
&& requires (T a, typename T::value_type c)
{
    a *= c;
    a*c;
    c*a;
};

template <typename T>
concept FloatingPointVectorOperable
    = VectorValued<T> && std::floating_point<typename T::value_type>;

template <typename T>
concept ArrayLike = requires (T x, std::size_t i)
{
    std::tuple_size<T>::value;
    x[i];
};

template <typename T>
    requires std::floating_point<T>
        || (ArrayLike<T> && FloatingPointVectorOperable<T>)
struct IntegralResult
{
    T val;
    T err;

    [[nodiscard]] constexpr std::size_t ndim() const noexcept
    {
        if constexpr (std::is_floating_point<T>::value)
            return 1;
        else
            return std::tuple_size<T>::value;
    }

    constexpr IntegralResult& operator+=(const IntegralResult& x)
    {
        val += x.val;
        err += x.err;
        return *this;
    }

    constexpr IntegralResult& operator-=(const IntegralResult& x)
    {
        val -= x.val;
        err -= x.err;
        return *this;
    }

    [[nodiscard]] constexpr IntegralResult
    operator+(const IntegralResult& x) const
    {
        IntegralResult res = *this;
        res += x;
        return res;
    }

    [[nodiscard]] constexpr IntegralResult
    operator-(const IntegralResult& x) const
    {
        IntegralResult res = *this;
        res -= x;
        return res;
    }
};

}