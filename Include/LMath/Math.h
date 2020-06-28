/*
Copyright (c) 2019 Lior Lahav

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

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

#ifndef __LMATH_MATH_H__
#define __LMATH_MATH_H__ 1

#include <type_traits>
#include <cmath>
#if defined(_WIN32)
	#include <intrin.h> // for bsr instruction.
#endif
namespace LMath
{
	namespace Constans
	{
		using ConstantType = double;
		constexpr ConstantType E = 2.71828182845904523536;// e
		constexpr ConstantType Log2e = 1.44269504088896340736;// log2(e)
		constexpr ConstantType Log10e = 0.434294481903251827651;// log10(e)
		constexpr ConstantType Ln2 = 0.693147180559945309417;// ln(2)
		constexpr ConstantType Ln10 = 2.30258509299404568402;// ln(10)
		constexpr ConstantType Pi = 3.14159265358979323846;// pi
		constexpr ConstantType TwoPi = 2.0 * Pi;
		constexpr ConstantType PiOver2 = 1.57079632679489661923;// pi/2
		constexpr ConstantType PiOver4 = 0.785398163397448309616;// pi/4
		constexpr ConstantType OneOverPi = 0.318309886183790671538;// 1/pi
		constexpr ConstantType TwoOverPi = 0.636619772367581343076;// 2/pi
		constexpr ConstantType TwoOverSqrtPi = 1.12837916709551257390;// 2/sqrt(pi)
		constexpr ConstantType Sqrt2 = 1.41421356237309504880;// sqrt(2)
		constexpr ConstantType OneOverSqrtw = 0.707106781186547524401;// 1/sqrt(2)
		constexpr ConstantType DegToRad = Pi / 180.0; 
		constexpr ConstantType RadToDeg = 180.0 / Pi;
		
	}

	template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	constexpr T Sign(T val)
	{
		return (static_cast<T>(0) < val) - (val < static_cast<T>(0));
	}

	// Mathematical Modulo for unsigned integral numbers - just operator %
	template <typename T, typename std::enable_if_t<std::is_integral_v<T> && std::is_unsigned_v<T>, int> = 0>
	constexpr T Modulo(T val, T mod)
	{
		return val % mod;
	}

	// Mathematical Modulo for signed integral numbers - [operator %] X 2
	template <typename T, typename std::enable_if_t<std::is_integral_v<T> && std::is_signed_v<T>, int> = 0>
	constexpr T Modulo(T val, T mod)
	{
		return (mod + (val % mod)) % mod;
	}

	template <typename T, typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
	constexpr T Modulo(T val, T mod)
	{
		return fmod(mod + fmod(val, mod), mod);
	}

	// Mathematical Modulo for unsigned integral numbers - just operator %
	template <typename T, typename std::enable_if_t<std::is_integral_v<T>, int> = 0>
	constexpr T ModuloOperator(T val, T mod)
	{
		return val % mod;
	}

	template <typename T, typename std::enable_if_t<std::is_floating_point_v<T>,int> = 0>
	constexpr T ModuloOperator(T val, T mod)
	{
		return fmod(val, mod);
	}


	template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	constexpr T ToDegrees(T val)
	{
		return val * Constans::RadToDeg;
	}

	template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
	constexpr  T ToRadians(T val)
	{
		return val * Constans::DegToRad;
	}

	template <class T>
	bool IsPowerOfTwo(T x)
	{
		static_assert(std::is_unsigned_v<T>, "Argument to IsPowerOfTwo must be unsigned.");
		static_assert(std::is_integral_v<T>, "Argument to IsPowerOfTwo must be an integral numeric type.");
		return (x & (~x + 1)) == x;
	}

	template <typename T>
	T NextPowerOfTwo(T num)
	{
		static_assert(std::is_unsigned_v<T>, "Argument to IsPowerOfTwo must be unsigned.");
		static_assert(std::is_integral_v<T>, "Argument to IsPowerOfTwo must be an integral numeric type.");

#if defined(_WIN32)
		unsigned long bitIndex;
		unsigned char res = _BitScanReverse(reinterpret_cast<unsigned long*>(&bitIndex), num - 1);
		return static_cast<T>(static_cast<size_t>(1) << (res * (bitIndex + 1)));
#elif defined(__GNUC__) ||  defined(__clang__)

		if (num < 2)
			return num;

		return  1 << (((sizeof(T) * 8 - 1) ^ __builtin_clz(num - 1)) + 1);
#else
		num--;
		num |= num >> 1;
		num |= num >> 2;
		num |= num >> 4;
		num |= num >> 8;
		num |= num >> 16;
		num++;
		return num;
#endif

	}
}
#endif