#pragma once

#include <type_traits>
#include <algorithm>
#include <stdexcept>
#include <limits>

//#define LL_VECTOR_OGRE_EXTENSIONS
#define LMATH_VECTOR_WINDOWS_EXTENSIONS

#ifdef LMATH_VECTOR_WINDOWS_EXTENSIONS
	#include <Windows.h>
	#pragma push_macro("min")
	#pragma push_macro("max")
	#undef min
	#undef max
#endif


#ifdef LMATH_VECTOR_OGRE_EXTENSIONS
	#include "OgreVector2.h"
	#include "OgreVector3.h"
#endif



template <typename T, size_t dim>
class VectorBase
{
public:
	using ElementType = typename T;
	inline static ElementType L_Epsilon = std::numeric_limits<ElementType>::epsilon();
	inline static constexpr ElementType L_Half = static_cast<ElementType>(0.5);
	inline static constexpr ElementType L_Zero = static_cast<ElementType>(0);
	inline static constexpr ElementType L_One = static_cast<ElementType>(1);
	inline static constexpr ElementType L_Two = static_cast<ElementType>(2);
	inline static constexpr ElementType L_Four = static_cast<ElementType>(4);
	
	
	inline static ElementType UninitializedValue = std::numeric_limits< ElementType>::min();
	static const VectorBase Uninitialized;
	static const VectorBase Zero;
	static const VectorBase Unit;

	
#pragma region Contructors

	template <typename... Args, typename = typename std::enable_if_t<sizeof...(Args) == dim - 1>>
	VectorBase(ElementType first, Args... args)
	{
		_Assign(0, static_cast<ElementType>(first), static_cast<ElementType>(args)...);
	}

	VectorBase(const ElementType* data)
	{
		memcpy(&mElements, data, sizeof(ElementType) * dim);
	}

	template <typename RHS_U, size_t RHS_DIM>
	//Convert from any vector to any vector 
	VectorBase(const VectorBase<RHS_U, RHS_DIM> & rhs)
	{
		const size_t elementsToConvert = std::min(dim, RHS_DIM);

		for (size_t i = 0; i < elementsToConvert; i++)
			at(i) = rhs.at(i);


		for (size_t i = elementsToConvert; i < dim; i++)
			at(i) = 0;
	}

	VectorBase() : VectorBase(UninitializedValue)
	{

	}

	VectorBase(ElementType initializationValue)
	{
		for (size_t i = 0; i < dim; i++)
			at(i) = initializationValue;
	}


	VectorBase(const std::initializer_list<ElementType>& initList)
	{
		if (initList.size() != dim)
			throw std::runtime_error("Error, wrong numbe of arguments");

		for (size_t i = 0; i < initList.size(); i++)
			at(i) = *(initList.begin() + i);
		//auto it = initList.begin();
		//size_t i = 0;
		
		/*while (it != initList.end())
		{
			at(i) = *it;
			it++;
			i++;
		}*/
	}
#pragma endregion
	
	ElementType& at(size_t idx)
	{
		return mElements[idx];
	}

	const ElementType& at(size_t idx) const
	{
		return mElements[idx];
	}

	bool IsInitialzied() const
	{
		for (size_t i = 0; i < dim; i++)
			if (at(i) == UninitializedValue)
				return false;

		return true;
	}

	ElementType DotProduct(const VectorBase& rhs) const
	{
		ElementType ret = 0;
		for (size_t i = 0; i < dim; i++)
			ret += at(i) * rhs.at(i);
		return ret;
	}

	template< typename = typename std::enable_if_t<std::is_floating_point_v<ElementType> >>
	VectorBase Round() const
	{
		VectorBase ret;

		for (size_t i = 0; i < dim; i++)
			ret.at(i) = std::round(at(i));

		return ret;
	}


	template< typename = typename std::enable_if_t< (dim == 3)>>
	VectorBase Cross(const VectorBase& rhs)
	{
		return
		{
			at(1) * rhs.at(2) - at(2) * rhs.at(1),
			at(2) * rhs.at(0) - at(0) * rhs.at(2),
			at(0) * rhs.at(1) - at(1) * rhs.at(0)
		};
	}

	template< typename = typename std::enable_if_t< (dim == 3)>>
	VectorBase Orthogonal() const
	{
		return at(2) < at(0) ? VectorBase(at(1), -at(0), 0) : VectorBase(0, -at(2), at(1));
	}



	
	
	
	ElementType NormSquared() const { return DotProduct(*this); }
	ElementType SquaredLength() const { return NormSquared(); }
	
	double Norm() const { return std::sqrt(DotProduct(*this)); }
	double Length() const { return Norm(); }

	
	bool isZeroLength() const
	{
		return SquaredLength() < std::numeric_limits<ElementType>::epsilon();
	}

	bool operator ==(const VectorBase& rhs) const
	{
		for (size_t i = 0; i < dim; i++)
			if (at(i) != rhs.at(i))
				return false;

		return true;
	}

	bool operator !=(const VectorBase& rhs) const
	{
		return !(*this == rhs);
	}

#pragma region Arithmetic

	VectorBase operator -() const
	{
		VectorBase vec;

		for (size_t i = 0; i < dim; i++)
			vec.at(i) = -at(i);
		return vec;
	}


	VectorBase operator -(ElementType value) const
	{
		VectorBase vec;

		for (size_t i = 0; i < dim; i++)
			vec.at(i) = at(i)  - value;
		return vec;
	}

	friend VectorBase operator-(ElementType value, const VectorBase& vec)
	{
		VectorBase ret;

		for (size_t i = 0; i < dim; i++)
			ret.at(i) = value - vec.at(i);
		return ret;
	}

	VectorBase operator -(const VectorBase& value)  const
	{
		VectorBase vec;
		for (size_t i = 0; i < dim; i++)
			vec.at(i) = at(i) - value.at(i);
		return vec;
	}

	VectorBase operator +(ElementType value) const
	{
		VectorBase vec;
		for (size_t i = 0; i < dim; i++)
			vec.at(i) = at(i) + value;

		return vec;
	}

	friend VectorBase operator +(ElementType value, const VectorBase& vec)
	{
		return  vec + value;
	}

	VectorBase operator +(const VectorBase& value) const
	{
		VectorBase vec;
		for (size_t i = 0; i < dim; i++)
			vec.at(i) = at(i) + value.at(i);

		return vec;
	}


	friend VectorBase operator*(ElementType value, const VectorBase& vec)
	{
		return vec * value;
	}


	VectorBase operator *(ElementType value) const
	{
		VectorBase vec;
		for (size_t i = 0; i < dim; i++)
			vec.at(i) =  at(i) * value;

		return vec;
	}

	VectorBase operator *(const VectorBase& value) const
	{
		VectorBase vec;
		for (size_t i = 0; i < dim; i++)
			vec.at(i) = at(i) * value.at(i);

		return vec;
	}



	friend VectorBase operator/(ElementType value, const VectorBase& vec)
	{
		VectorBase ret;
		for (size_t i = 0; i < dim; i++)
			ret.at(i) = value / vec.at(i);

		return ret;
	}


	VectorBase operator /(ElementType value) const
	{
		VectorBase vec;
		for (size_t i = 0; i < dim; i++)
			vec.at(i) = at(i) / value;

		return vec;
	}

	VectorBase operator /(const VectorBase& value) const
	{
		VectorBase vec;
		for (size_t i = 0; i < dim; i++)
			vec.at(i) = at(i) / value.at(i);

		return vec;
	}
#pragma endregion

	VectorBase Max(const VectorBase& rhs) const 
	{
		VectorBase ret;
		for (size_t i = 0; i < dim; i++)
			ret.at(i) = std::max(at(i), rhs.at(i));
		
		return ret;
	}

	VectorBase Min(const VectorBase& rhs) const
	{
		VectorBase ret;
		for (size_t i = 0; i < dim; i++)
			ret.at(i) = std::min(at(i), rhs.at(i));

		return ret;
	}


	VectorBase Normalize() const
	{
		double norm = Norm();
		return norm > 0 ? *this / norm : Zero;
	}

	VectorBase Reflect(const VectorBase& normal) const
	{
		return *this - static_cast<VectorBase::ElementType>(2) * DotProduct(normal) * normal;
	}


#ifdef LMATH_VECTOR_WINDOWS_EXTENSIONS
	static constexpr bool Is2DimIntegral = dim == 2 && std::is_integral_v<ElementType>;

	template< typename = typename std::enable_if_t<Is2DimIntegral>>
	operator POINT() const
	{
		return { at(0), at(1)};
	}

	template< typename = typename std::enable_if_t<Is2DimIntegral>>
	VectorBase(const POINT& p)
	{
		at(0) = p.x;
		at(1) = p.y;
	}

	template< typename = typename std::enable_if_t<Is2DimIntegral>>
	operator SIZE() const
	{
		return { at(0), at(1) };
	}

	template< typename = typename std::enable_if_t<Is2DimIntegral>>
	VectorBase(const SIZE& p)
	{
		at(0) = p.cx;
		at(1) = p.cy;
	}
	
#endif


	
#ifdef LMATH_VECTOR_OGRE_EXTENSIONS
	template< typename = typename std::enable_if_t< (dim == 2)>>
	Ogre::Vector2 ToVector2() const
	{

		return { at(0) ,at(1) };
	}

	template< typename = typename std::enable_if_t< (dim == 3)>>
	Ogre::Vector2 ToVector3() const
	{
		return { at(0) ,at(1), at(2) };
	}

	template< typename = typename std::enable_if_t< (dim == 4)>>
	Ogre::Vector2 ToVector4() const
	{
		return { at(0) ,at(1), at(2),at(3) };
	}
#endif

#pragma region EuclidSpaceNotation
	template< typename = typename std::enable_if_t< (dim >= 1)>>
	ElementType& X() 
	{
		return at(0);
	}

	template< typename = typename std::enable_if_t< (dim >= 1)>>
	const ElementType& X() const
	{
		return at(0);
	}

	template< typename = typename std::enable_if_t< (dim >= 2)>>
	ElementType & Y()
	{
		return at(1);
	}

	template< typename = typename std::enable_if_t< (dim >= 2)>>
	const ElementType & Y() const
	{
		return at(1);
	}

	template< typename = typename std::enable_if_t< (dim >= 3)>>
	ElementType & Z()
	{
		return at(2);
	}

	template< typename = typename std::enable_if_t< (dim >= 3)>>
	const ElementType & Z() const
	{
		return at(2);
	}

	template< typename = typename std::enable_if_t< (dim >= 4)>>
	ElementType& W()
	{
		return at(3);
	}

	template< typename = typename std::enable_if_t< (dim >= 4)>>
	const ElementType& W() const
	{
		return at(3);
	}
#pragma endregion

#pragma region Private helper methods
	private:
		//Vector consturctor helper, faster than using initializer list.
		void _Assign(size_t index, ElementType element)
		{
			at(index) = element;
		}

		template <typename ElementType, typename ...ARGS>
		void _Assign(size_t index, ElementType element, ARGS... args)
		{
			at(index) = element;
			_Assign(index + 1, args...);
		}

#pragma endregion
		

private:
	ElementType mElements[dim];
};

template <int dim>
using VectorI32 = VectorBase<int32_t, dim>;

template <int dim>
using VectorI64 = VectorBase<int64_t, dim>;

template <int dim>
using VectorF64 = VectorBase<double, dim>;

template <int dim>
using VectorF32 = VectorBase<float, dim>;

template <class T, size_t dim>
const VectorBase<T, dim> VectorBase<T, dim>::Uninitialized;


template <class T, size_t dim>
const VectorBase<T, dim> VectorBase<T, dim>::Zero(VectorBase<T, dim>::L_Zero);

template <class T, size_t dim>
const VectorBase<T, dim> VectorBase<T, dim>::Unit(VectorBase<T, dim>::L_One);


#ifdef LMATH_VECTOR_WINDOWS_EXTENSIONS
#pragma pop_macro("min")
#pragma pop_macro("max")
#endif