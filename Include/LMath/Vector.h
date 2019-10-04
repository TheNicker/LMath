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


#ifndef __LMATH_VECTOR_H__
#define __LMATH_VECTOR_H__ 1

#include "LMathConfig.h"
#include "Element.h"
#include "Math.h"

#if LMATH_VECTOR_WINDOWS_EXTENSIONS == 1
	#include <Windows.h>
#endif


#if LMATH_VECTOR_OGRE_EXTENSIONS == 1
#include <OgreVector2.h>
#include <OgreVector3.h>
#endif

namespace LMath
{
	template <typename T, size_t DIM >
	class VectorBaseCartesian : public ElementBase<T, DIM>
	{
		using BaseClass = ElementBase<T, DIM>;
		using BaseClass::BaseClass;
	};

#if LMATH_VECTOR_CARTESIAN_COMPONENT & LMATH_VECTOR_CARTESIAN_COMPONENT_PLAIN_VARIABLE_BIG_CASE
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_0_UPPER_CASE ElementType& X = this->at(0);
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_1_UPPER_CASE ElementType& Y = this->at(1);
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_2_UPPER_CASE ElementType& Z = this->at(2);
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_3_UPPER_CASE ElementType& W = this->at(3);
#else
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_0_UPPER_CASE 
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_1_UPPER_CASE 
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_2_UPPER_CASE 
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_3_UPPER_CASE 
#endif

#if LMATH_VECTOR_CARTESIAN_COMPONENT & LMATH_VECTOR_CARTESIAN_COMPONENT_PLAIN_VARIABLE_LOW_CASE
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_0_LOWER_CASE ElementType& x = this->at(0);
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_1_LOWER_CASE ElementType& y = this->at(1);
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_2_LOWER_CASE ElementType& z = this->at(2);
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_3_LOWER_CASE ElementType& w = this->at(3);
#else
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_0_LOWER_CASE 
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_1_LOWER_CASE 
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_2_LOWER_CASE 
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_3_LOWER_CASE 
#endif

#if LMATH_VECTOR_CARTESIAN_COMPONENT & LMATH_VECTOR_CARTESIAN_COMPONENT_GETTER
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_0_GETTER ElementType& X() {return this->at(0);} const ElementType& X() const {return this->at(0);}
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_1_GETTER ElementType& Y() {return this->at(1);} const ElementType& Y() const {return this->at(1);}
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_2_GETTER ElementType& Z() {return this->at(2);} const ElementType& Z() const {return this->at(2);}
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_3_GETTER ElementType& W() {return this->at(3);} const ElementType& W() const {return this->at(3);}

#else
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_0_GETTER
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_1_GETTER
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_2_GETTER
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_3_GETTER
#endif


#define LMATH_DEFINE_CARTESIAN_COMPONENT_0 \
LMATH_DEFINE_CARTESIAN_COMPONENT_0_LOWER_CASE \
LMATH_DEFINE_CARTESIAN_COMPONENT_0_UPPER_CASE \
LMATH_DEFINE_CARTESIAN_COMPONENT_0_GETTER

#define LMATH_DEFINE_CARTESIAN_COMPONENT_1 \
LMATH_DEFINE_CARTESIAN_COMPONENT_1_LOWER_CASE \
LMATH_DEFINE_CARTESIAN_COMPONENT_1_UPPER_CASE \
LMATH_DEFINE_CARTESIAN_COMPONENT_1_GETTER

#define LMATH_DEFINE_CARTESIAN_COMPONENT_2 \
LMATH_DEFINE_CARTESIAN_COMPONENT_2_LOWER_CASE \
LMATH_DEFINE_CARTESIAN_COMPONENT_2_UPPER_CASE \
LMATH_DEFINE_CARTESIAN_COMPONENT_2_GETTER

#define LMATH_DEFINE_CARTESIAN_COMPONENT_3 \
LMATH_DEFINE_CARTESIAN_COMPONENT_3_LOWER_CASE \
LMATH_DEFINE_CARTESIAN_COMPONENT_3_UPPER_CASE \
LMATH_DEFINE_CARTESIAN_COMPONENT_3_GETTER


	template <typename T>
	class VectorBaseCartesian<T, 1> : public ElementBase<T, 1>
	{
	public:
		using BaseClass = ElementBase<T, 1>;
		using BaseClass::BaseClass;
		using ElementType = typename BaseClass::ElementType;
		VectorBaseCartesian(const BaseClass& base) : BaseClass(base) {}
		LMATH_DEFINE_CARTESIAN_COMPONENT_0
	};


	template <typename T>
	class VectorBaseCartesian<T, 2> : public ElementBase<T, 2>
	{
	public:
		using BaseClass = ElementBase<T, 2>;
		using BaseClass::BaseClass;
		using ElementType = typename BaseClass::ElementType;
		VectorBaseCartesian(const BaseClass& base) : BaseClass(base) {}
		LMATH_DEFINE_CARTESIAN_COMPONENT_0
		LMATH_DEFINE_CARTESIAN_COMPONENT_1

		static const VectorBaseCartesian Up;
		static const VectorBaseCartesian Down;
		static const VectorBaseCartesian Left;
		static const VectorBaseCartesian Right;
	};

	template <typename T>
	class VectorBaseCartesian<T, 3> : public ElementBase<T, 3>
	{
	public:
		using BaseClass = ElementBase<T, 3>;
		using BaseClass::BaseClass;
		using BaseClass::at;
		using ElementType = typename BaseClass::ElementType;
		VectorBaseCartesian(const BaseClass& base) : BaseClass(base){}

		LMATH_DEFINE_CARTESIAN_COMPONENT_0
		LMATH_DEFINE_CARTESIAN_COMPONENT_1
		LMATH_DEFINE_CARTESIAN_COMPONENT_2

		static const VectorBaseCartesian Down;
		static const VectorBaseCartesian Up;
		static const VectorBaseCartesian Left;
		static const VectorBaseCartesian Right;
		static const VectorBaseCartesian Forward;
		static const VectorBaseCartesian Backward;
	};

	template <typename T>
	class VectorBaseCartesian<T, 4> : public ElementBase<T, 4>
	{
	public:
		using BaseClass = ElementBase<T, 4>;
		using BaseClass::BaseClass;
		using ElementType = typename BaseClass::ElementType;
		VectorBaseCartesian(const BaseClass& base) : BaseClass(base) {}

		LMATH_DEFINE_CARTESIAN_COMPONENT_0
		LMATH_DEFINE_CARTESIAN_COMPONENT_1
		LMATH_DEFINE_CARTESIAN_COMPONENT_2
		LMATH_DEFINE_CARTESIAN_COMPONENT_3

	};


	template <typename T, size_t DIM>
	class VectorBase : public VectorBaseCartesian<T, DIM>
	{
	public:
		using BaseClass = VectorBaseCartesian<T, DIM>;
		using VectorBaseCartesian<T, DIM>::VectorBaseCartesian;
		using BaseClass::at;
		using ElementType = typename BaseClass::ElementType;
		using BaseClass::L_One;
		using BaseClass::L_Two;
		
		static const VectorBase Zero;
		static const VectorBase Unit;

        VectorBase(const BaseClass& base) : BaseClass(base)
        {

        }


		bool operator ==(const VectorBase& rhs) const
		{
			for (size_t i = 0; i < DIM; i++)
				if (at(i) != rhs.at(i))
					return false;

			return true;
		}

		bool operator !=(const VectorBase& rhs) const
		{
			return !(*this == rhs);
		}

		ElementType Dot(const VectorBase& rhs) const
		{
			ElementType ret = 0;
			for (size_t i = 0; i < DIM; i++)
				ret += at(i) * rhs.at(i);
			return ret;
		}

		VectorBase Round() const
		{
			VectorBase ret;

			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = std::round(at(i));

			return ret;
		}


		

		ElementType NormSquared() const { return Dot(*this); }
		ElementType SquaredLength() const { return NormSquared(); }


		double Norm() const { return std::sqrt(NormSquared()); }
		double Length() const { return Norm(); }


		bool isZeroLength() const
		{
			return SquaredLength() < std::numeric_limits<ElementType>::epsilon();
		}

		VectorBase Normalized() const
		{
			double norm = Norm();
			return norm > 0 ? *this / norm : Zero;
		}
	

#pragma region Arithmetic

		
		VectorBase operator -() const
		{
			VectorBase vec;

			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = -at(i);
			return vec;
		}

		//operator -
		VectorBase operator -(ElementType value) const
		{
			VectorBase vec;

			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) - value;
			return vec;
		}

		friend VectorBase operator-(ElementType value, const VectorBase& vec)
		{
			VectorBase ret;

			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = value - vec.at(i);
			return ret;
		}

		VectorBase operator -(const VectorBase& value)  const
		{
			VectorBase vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) - value.at(i);
			return vec;
		}

		VectorBase& operator -=(ElementType value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) - value;

			return *this;
		}

		VectorBase& operator -=(const VectorBase& value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) - value.at(i);

			return *this;
		}


		// opeator  +
		VectorBase operator +(ElementType value) const
		{
			VectorBase vec;
			for (size_t i = 0; i < DIM; i++)
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
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) + value.at(i);

			return vec;
		}


		VectorBase& operator +=(ElementType value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) + value;

			return *this;
		}

		VectorBase& operator +=(const VectorBase& value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) + value.at(i);

			return *this;
		}


		//operator *
		friend VectorBase operator*(ElementType value, const VectorBase& vec)
		{
			return vec * value;
		}


		VectorBase operator *(ElementType value) const
		{
			VectorBase vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) * value;

			return vec;
		}

		VectorBase& operator *=(ElementType value) 
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) * value;

			return *this;
		}

		VectorBase& operator *=(const VectorBase& value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) * value.at(i);

			return *this;
		}


		VectorBase operator *(const VectorBase& value) const
		{
			VectorBase vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) * value.at(i);

			return vec;
		}


		// Operator /
		friend VectorBase operator/(ElementType value, const VectorBase& vec)
		{
			VectorBase ret;
			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = value / vec.at(i);

			return ret;
		}

		VectorBase& operator /=(const VectorBase& value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) / value.at(i);

			return *this;
		}

		VectorBase& operator /=(const ElementType value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) / value;

			return *this;
		}


		VectorBase operator /(ElementType value) const
		{
			VectorBase vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) / value;

			return vec;
		}

		VectorBase operator /(const VectorBase& value) const
		{
			VectorBase vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) / value.at(i);

			return vec;
		}
		
		// Operator %
		VectorBase operator %(const VectorBase& value) const
		{
			VectorBase vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) & value.at(i);
			return vec;
		}

		VectorBase operator %(ElementType value) const
		{
			VectorBase vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) & value;
			return vec;
		}

		friend VectorBase operator %(ElementType value, const VectorBase& vec)
		{
			VectorBase ret;
			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = value % at(i);
			return ret;
		}

		VectorBase& operator %=(ElementType value) 
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) &= value;
			return this;
		}

		VectorBase& operator %=(const VectorBase& value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) &= value.at(i);
			return this;
		}

	


		// mathematic modulo operation
		VectorBase Modulo(const VectorBase& value)
		{
			VectorBase vec;
			for (size_t i = 0; i < DIM; i++)
				at(i) = LMath::Modulo(at(i), value.at(i));
			
			return vec;
		}

		VectorBase Modulo(ElementType value)
		{
			VectorBase vec;
			for (size_t i = 0; i < DIM; i++)
				at(i) = LMath::Modulo(at(i), value);
			return vec;
		}

#pragma endregion

		VectorBase Max(const VectorBase& rhs) const
		{
			VectorBase ret;
			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = (std::max)(at(i), rhs.at(i));

			return ret;
		}

		VectorBase Min(const VectorBase& rhs) const
		{
			VectorBase ret;
			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = (std::min)(at(i), rhs.at(i));

			return ret;
		}


		ElementType Distance(const VectorBase& rhs) const
		{
			return (*this - rhs).Length();
		}


		VectorBase Lerp(const VectorBase& rhs , double t) const
		{
			return (rhs - *this) * t + *this;
		}

		VectorBase Slerp(const VectorBase& rhs, double t) const
		{
			VectorBase a = *this;
			VectorBase b = rhs;
			double magA = a.Length();
			double magB = b.Length();
			a /= magA;
			b /= magB;
			double dot = std::clamp(a.Dot(b), -L_One, L_One);
			double theta = acos(dot) * t;
			VectorBase relativeVec = (b - a * dot).Normalized();
			VectorBase newVec = a * std::cos(theta) + relativeVec * std::sin(theta);
			return newVec * (magA + (magB - magA) * t);
		}


		double Angle(const VectorBase& rhs) const
		{
			return std::acos(std::clamp(Dot(rhs) / (Length() * rhs.Length()), -L_One, L_One));
		}

#pragma region Vector2/3 specific

		


		VectorBase Cross(const VectorBase& rhs) const
		{
			static_assert(DIM == 3, "not a vector with 3 DIMensions");
			
			return
			{
				at(1) * rhs.at(2) - at(2) * rhs.at(1),
				at(2) * rhs.at(0) - at(0) * rhs.at(2),
				at(0) * rhs.at(1) - at(1) * rhs.at(0)
			};
		}

		VectorBase Orthogonal() const
		{
			static_assert(DIM == 3, "not a vector with 3 DIMensions");
			return at(2) < at(0) ? VectorBase(at(1), -at(0), 0) : VectorBase(0, -at(2), at(1));
		}

		

		VectorBase Reflect(const VectorBase& normal) const
		{
			static_assert(DIM == 3 || DIM == 2, "not a vector with two or three DIMensions");
			return *this -  L_Two * Dot(normal)* normal;
		}

		

#pragma endregion


#if LMATH_VECTOR_WINDOWS_EXTENSIONS == 1
		operator POINT() const
		{
            static_assert(DIM == 2, "not a vector with 2 DIMensions");
			return { at(0), at(1) };
		}

		VectorBase(const POINT& p)
		{
            static_assert(DIM == 2, "not a vector with 2 DIMensions");
			at(0) = p.x;
			at(1) = p.y;
		}
		
		operator SIZE() const
		{
            static_assert(DIM == 2, "not a vector with 2 DIMensions");
			return { at(0), at(1) };
		}
		
		VectorBase(const SIZE& p)
		{
            static_assert(DIM == 2, "not a vector with 2 DIMensions");
			at(0) = p.cx;
			at(1) = p.cy;
		}

#endif



#if LMATH_VECTOR_OGRE_EXTENSIONS == 1
        VectorBase(const Ogre::Vector2 & v) 
        {
            static_assert(DIM == 2, "not a vector with 2 DIMensions");
            return (v.x, v.y);
        }

        operator Ogre::Vector2() const
        {
            static_assert(DIM == 2, "not a vector with 2 DIMensions");
            return { at(0) ,at(1) };
        }

        VectorBase(const Ogre::Vector3& v)
        {
            static_assert(DIM == 3, "not a vector with 3 DIMensions");
            return (v.x, v.y, v.z);
        }

        operator Ogre::Vector3() const
        {
            static_assert(DIM == 3, "not a vector with 3 DIMensions");
            return { at(0) ,at(1) , at(2) };
        }

        VectorBase(const Ogre::Vector4& v)
        {
            static_assert(DIM == 4, "not a vector with 4 DIMensions");
            return (v.x, v.y,v.z,v.w);
        }

        operator Ogre::Vector4() const
        {
            static_assert(DIM == 4, "not a vector with 4 DIMensions");
            return { at(0) ,at(1), at(2), at(3) };
        }
#endif



	};
	template <class T>	const VectorBaseCartesian<T, 2> VectorBaseCartesian<T, 2>::Up (0, 1);
	template <class T>	const VectorBaseCartesian<T, 2> VectorBaseCartesian<T, 2>::Down(0,-1);
	template <class T>	const VectorBaseCartesian<T, 2> VectorBaseCartesian<T, 2>::Left(-1, 0);
	template <class T>	const VectorBaseCartesian<T, 2> VectorBaseCartesian<T, 2>::Right(1, 0);


	template <class T>	const VectorBaseCartesian<T, 3> VectorBaseCartesian<T, 3>::Up(0, 1, 0);
	template <class T>	const VectorBaseCartesian<T, 3> VectorBaseCartesian<T, 3>::Down(0, -1, 0);
	template <class T>	const VectorBaseCartesian<T, 3> VectorBaseCartesian<T, 3>::Forward(0, 0, 1);
	template <class T>	const VectorBaseCartesian<T, 3> VectorBaseCartesian<T, 3>::Backward(0, 0, -1);
	template <class T>	const VectorBaseCartesian<T, 3> VectorBaseCartesian<T, 3>::Left(-1, 0, 0);
	template <class T>	const VectorBaseCartesian<T, 3> VectorBaseCartesian<T, 3>::Right(1, 0 ,0);

	template <class T, size_t DIM>
	const VectorBase<T, DIM> VectorBase<T, DIM>::Zero(VectorBase<T, DIM>::L_Zero);

	template <class T, size_t DIM>
	const VectorBase<T, DIM> VectorBase<T, DIM>::Unit(VectorBase<T, DIM>::L_One);
}

#endif