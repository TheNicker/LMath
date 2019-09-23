#pragma once

#include "Vector.h"

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>, int> = 0>
static constexpr T Sign(T val)
{
	return (static_cast<T>(0) < val) - (val < static_cast<T>(0));
}

template <class T,size_t dim>
class AreaBase
{
public:
	enum class IntersectState { Disjoint, Touch, Overlap };
	using VectorType = VectorBase<T, dim>;
	using ElementType = typename VectorBase<T, dim>::ElementType;
	static const AreaBase Uninitialized;
	static const size_t Dimensions = dim;


	AreaBase() 
	{
		Reset();
	}

	AreaBase(const std::initializer_list<VectorType>& initializer ) 
	{
		if (initializer.size() != 2)
		{
			throw "Sf";
		}
		
		minValue = *(initializer.begin());
		maxValue = *(initializer.begin() + 1);
	}


	AreaBase (const VectorType& min, const VectorType& max)
	{
		minValue = min;
		maxValue = max;
	}

	void Reset()
	{
		minValue = VectorType();
		maxValue = VectorType();
	}

	VectorType GetMin() const
	{
		return minValue;
	}

	VectorType GetMax() const
	{
		return maxValue;
	}

	AreaBase Intersection(const AreaBase& rhs) const
	{
		return { minValue.Max(rhs.minValue), maxValue.Min(rhs.maxValue) };
	}

	AreaBase Inflate(ElementType delta) const
	{
		return { minValue - static_cast<ElementType>(delta) , maxValue + static_cast <ElementType>(delta) };
	}

	

	bool IntersectsIncludeBoundry(const AreaBase& rhs) const
	{
		bool intersects = true;
		size_t i = 0;
		while (intersects == true && i < dim)
		{
			intersects &= maxValue.at(i) >= rhs.minValue.at(i) && minValue.at(i) <= rhs.maxValue.at(i);
			i++;
		}

		return intersects;
	}

	void computeIntersection(int sign, IntersectState& currentState) const
	{
		switch (sign)
		{
		case -1:
			//Disjoint
			currentState = IntersectState::Disjoint;
			break;
		case 0:
			// Touch
			currentState = std::min(currentState, IntersectState::Touch);
			break;
			//overlap
		case 1:
			currentState = std::min(currentState, IntersectState::Overlap);
			break;
		}

	}

	IntersectState Contains(const AreaBase& rhs) const
	{
		IntersectState result = IntersectState::Overlap;
		size_t i = 0;
		while (result != IntersectState::Disjoint && i < dim)
		{
			const ElementType sign = Sign( (maxValue.at(i) - rhs.minValue.at(i)) * (rhs.maxValue.at(i) - minValue.at(i)));
			computeIntersection(sign, result);
			i++;
		}

		return result;
	}


	IntersectState Contains(const VectorType& rhs) const
	{
		IntersectState result = IntersectState::Overlap;
		size_t i = 0;
		while (result != IntersectState::Disjoint && i < dim)
		{
			const ElementType sign = Sign((maxValue.at(i) - rhs.at(i)) * (rhs.at(i) - minValue.at(i)));
			computeIntersection(sign, result);
			i++;
		}

		return result;
	}

	[[deprecated("deprectead use Contains instead")]]
	bool Intersects(const AreaBase& rhs) const
	{
		bool intersects = true;
		size_t i = 0;
		while (intersects == true && i < dim)
		{
			intersects &= maxValue.at(i) > rhs.minValue.at(i) && minValue.at(i) < rhs.maxValue.at(i);
			i++;
		}
		
		return intersects;
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
		while (inside == true && i < dim)
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
		while (inside == true && i < dim)
		{
			inside &= minValue.at(i) <= vec.at(i) && maxValue.at(i) >= vec.at(i);
			i++;
		}

		return inside;
	}


	bool IsInitialized() const
	{
		return minValue.IsInitialzied() && maxValue.IsInitialzied();
	}
	
	ElementType GetSpan(size_t dim) const
	{
		return maxValue.at(dim) - minValue.at(dim);
	}

	AreaBase Round()
	{
		return { minValue.Round(), maxValue.Round() };
	}

	//check that min value is equals or larger than max value.
	bool IsValid() const
	{
		for (size_t i = 0; i < dim; i++)
			if (minValue.at(i) > maxValue.at(i))
				return false;

		return true;
	}

	bool operator ==(const AreaBase& rhs) const
	{
		return minValue == rhs.minValue && maxValue == rhs.maxValue;
	}

	AreaBase Join(const AreaBase& rhs) const
	{
		return { minValue.Min(rhs.minValue) , { maxValue.Max(rhs.maxValue)} };
	}

private:
	VectorType minValue;
	VectorType maxValue;
};


template <class T, size_t dim>
const AreaBase<T, dim> AreaBase<T, dim>::Uninitialized;