/*
Copyright (c) 2024 Sebastian Sassi

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do 
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
SOFTWARE.
*/
#pragma once

#include <concepts>
#include <utility>

namespace cubage
{

template <typename T>
concept VectorValued = requires (T& a, T b)
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
concept RealVector = std::floating_point<T>
        || (FloatingPointVectorOperable<T> && ArrayLike<T>);

template <typename ValueType, typename StatusType>
struct Result
{
    using value_type = ValueType;
    using status_type = StatusType;

    value_type value;
    status_type status;
};

enum class Status
{
    SUCCESS,
    MAX_SUBDIV
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

    constexpr IntegralResult& operator+=(const IntegralResult& x) noexcept
    {
        val += x.val;
        err += x.err;
        return *this;
    }

    constexpr IntegralResult& operator-=(const IntegralResult& x) noexcept
    {
        val -= x.val;
        err -= x.err;
        return *this;
    }

    [[nodiscard]] constexpr IntegralResult
    operator+(const IntegralResult& x) const noexcept
    {
        return IntegralResult{this->val + x.val, this->err + x.err};
    }

    [[nodiscard]] constexpr IntegralResult
    operator-(const IntegralResult& x) const noexcept
    {
        return IntegralResult{this->val - x.val, this->err - x.err};
    }
};

}