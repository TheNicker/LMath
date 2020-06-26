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


#ifndef __LMATH_BOUNDS_H__
#define __LMATH_BOUNDS_H__ 1

#include "Vector.h"
#include "Math.h"

namespace LMath
{

	template <class T, size_t DIM>
	class BoundsBase
	{
	public:
		enum class Relation { Disjoint, Adjacent, Overlap };
		using VectorType = VectorBase<T, DIM>;
		using ElementType = typename VectorBase<T, DIM>::ElementType;
		inline static const size_t Dimensions = DIM;

		BoundsBase()
		{
			
		}

		BoundsBase(const VectorType& min, const VectorType& max)
		{
			minValue = min;
			maxValue = max;
		}
		

		VectorType GetMin() const
		{
			return minValue;
		}

		VectorType GetMax() const
		{
			return maxValue;
		}

		BoundsBase Intersection(const BoundsBase& rhs) const
		{
			return { minValue.Max(rhs.minValue), maxValue.Min(rhs.maxValue) };
		}

		BoundsBase Inflate(ElementType delta) const
		{
			return { minValue - static_cast<ElementType>(delta) , maxValue + static_cast <ElementType>(delta) };
		}



		bool IntersectsIncludeBoundry(const BoundsBase& rhs) const
		{
			bool intersects = true;
			size_t i = 0;
			while (intersects == true && i < DIM)
			{
				intersects &= maxValue.at(i) >= rhs.minValue.at(i) && minValue.at(i) <= rhs.maxValue.at(i);
				i++;
			}

			return intersects;
		}

		void computeIntersection(ElementType sign, Relation& currentState) const
		{
			switch (static_cast<int8_t>(sign))
			{
			case -1:
				//Disjoint
				currentState = Relation::Disjoint;
				break;
			case 0:
				// Adjacent
				currentState = (std::min)(currentState, Relation::Adjacent);
				break;
				//overlap
			case 1:
				currentState = (std::min)(currentState, Relation::Overlap);
				break;
			}

		}

		Relation RelationTo(const BoundsBase& rhs) const
		{
			Relation result = Relation::Overlap;
			size_t i = 0;
			while (result != Relation::Disjoint && i < DIM)
			{
				const ElementType sign = Sign((maxValue.at(i) - rhs.minValue.at(i)) * (rhs.maxValue.at(i) - minValue.at(i)));
				computeIntersection(sign, result);
				i++;
			}

			return result;
		}


		Relation RelationTo(const VectorType& rhs) const
		{
			Relation result = Relation::Overlap;
			size_t i = 0;
			while (result != Relation::Disjoint && i < DIM)
			{
				const ElementType sign = Sign((maxValue.at(i) - rhs.at(i)) * (rhs.at(i) - minValue.at(i)));
				computeIntersection(sign, result);
				i++;
			}

			return result;
		}

		template<typename... Args>
		bool IsInside(ElementType first, const Args... args) const
		{
			auto v = VectorType(first, args...);
			return IsInside(v);
		}



		bool IsInside(const VectorType& vec) const
		{
			bool inside = true;
			size_t i = 0;
			while (inside == true && i < DIM)
			{
				inside &= minValue.at(i) < vec.at(i) && maxValue.at(i) > vec.at(i);
				i++;
			}

			return inside;
		}

		bool IsInsideOrBoundry(const VectorType& vec) const
		{
			bool inside = true;
			size_t i = 0;
			while (inside == true && i < DIM)
			{
				inside &= minValue.at(i) <= vec.at(i) && maxValue.at(i) >= vec.at(i);
				i++;
			}

			return inside;
		}
		

		ElementType GetSpan(size_t dim) const
		{
			return maxValue.at(dim) - minValue.at(dim);
		}

		BoundsBase Round()
		{
			return { minValue.Round(), maxValue.Round() };
		}

		//check that min value is equals or larger than max value.
		bool IsValid() const
		{
			for (size_t i = 0; i < DIM; i++)
				if (minValue.at(i) > maxValue.at(i))
					return false;

			return true;
		}


		bool operator ==(const BoundsBase& rhs) const
		{
			return minValue == rhs.minValue && maxValue == rhs.maxValue;
		}

		BoundsBase Join(const BoundsBase& rhs) const
		{
			return { minValue.Min(rhs.minValue) , { maxValue.Max(rhs.maxValue)} };
		}

	private:
		VectorType minValue;
		VectorType maxValue;
	};
}
#endif