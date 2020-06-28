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
	#define LMATH_DEFINE_CARTESIAN_COMPONENT_0_GETTER ElementType& X() {static_assert(DIM > 0, "improper use"); return this->at(0);} const ElementType& X() const {static_assert(DIM > 0, "improper use");return this->at(0);}
#define LMATH_DEFINE_CARTESIAN_COMPONENT_1_GETTER ElementType& Y() {static_assert(DIM > 1, "improper use"); return this->at(1);} const ElementType& Y() const {static_assert(DIM > 1, "improper use");return this->at(1);}
#define LMATH_DEFINE_CARTESIAN_COMPONENT_2_GETTER ElementType& Z() {static_assert(DIM > 2, "improper use"); return this->at(2);} const ElementType& Z() const {static_assert(DIM > 2, "improper use");return this->at(2);}
#define LMATH_DEFINE_CARTESIAN_COMPONENT_3_GETTER ElementType& W() {static_assert(DIM > 3, "improper use"); return this->at(3);} const ElementType& W() const {static_assert(DIM > 3, "improper use");return this->at(3);}

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


	template <typename T, size_t DIM>
	class VectorBaseTemplate : public ElementBase<T, DIM>
	{
	public:
		using BaseClass = ElementBase<T, DIM>;
		using ElementBase<T, DIM>::ElementBase;
		using BaseClass::at;
		using ElementType = typename BaseClass::ElementType;
		using BaseClass::L_One;
		using BaseClass::L_Two;
		
		static const VectorBaseTemplate Zero;
		static const VectorBaseTemplate Unit;

		
		LMATH_DEFINE_CARTESIAN_COMPONENT_0
		LMATH_DEFINE_CARTESIAN_COMPONENT_1
		LMATH_DEFINE_CARTESIAN_COMPONENT_2
		LMATH_DEFINE_CARTESIAN_COMPONENT_3


        VectorBaseTemplate(const BaseClass& base) : BaseClass(base)
        {

        }

	

		bool operator ==(const VectorBaseTemplate& rhs) const
		{
			for (size_t i = 0; i < DIM; i++)
				if (at(i) != rhs.at(i))
					return false;

			return true;
		}

		bool operator !=(const VectorBaseTemplate& rhs) const
		{
			return !(*this == rhs);
		}

		ElementType Dot(const VectorBaseTemplate& rhs) const
		{
			ElementType ret = 0;
			for (size_t i = 0; i < DIM; i++)
				ret += at(i) * rhs.at(i);
			return ret;
		}

		VectorBaseTemplate Abs() const
		{
			VectorBaseTemplate ret;

			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = std::abs(at(i));

			return ret;
		}

		VectorBaseTemplate Round() const
		{
			VectorBaseTemplate ret;

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

		VectorBaseTemplate Normalized() const
		{
			double norm = Norm();
			return norm > 0 ? *this / norm : Zero;
		}
	

#pragma region Arithmetic

		
		VectorBaseTemplate operator -() const
		{
			VectorBaseTemplate vec;

			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = -at(i);
			return vec;
		}

		//operator -
		VectorBaseTemplate operator -(ElementType value) const
		{
			VectorBaseTemplate vec;

			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) - value;
			return vec;
		}

		friend VectorBaseTemplate operator-(ElementType value, const VectorBaseTemplate& vec)
		{
			VectorBaseTemplate ret;

			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = value - vec.at(i);
			return ret;
		}

		VectorBaseTemplate operator -(const VectorBaseTemplate& value)  const
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) - value.at(i);
			return vec;
		}

		VectorBaseTemplate& operator -=(ElementType value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) - value;

			return *this;
		}

		VectorBaseTemplate& operator -=(const VectorBaseTemplate& value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) - value.at(i);

			return *this;
		}


		// opeator  +
		VectorBaseTemplate operator +(ElementType value) const
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) + value;

			return vec;
		}

		friend VectorBaseTemplate operator +(ElementType value, const VectorBaseTemplate& vec)
		{
			return  vec + value;
		}

		VectorBaseTemplate operator +(const VectorBaseTemplate& value) const
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) + value.at(i);

			return vec;
		}


		VectorBaseTemplate& operator +=(ElementType value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) + value;

			return *this;
		}

		VectorBaseTemplate& operator +=(const VectorBaseTemplate& value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) + value.at(i);

			return *this;
		}


		//operator *
		friend VectorBaseTemplate operator*(ElementType value, const VectorBaseTemplate& vec)
		{
			return vec * value;
		}


		VectorBaseTemplate operator *(ElementType value) const
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) * value;

			return vec;
		}

		VectorBaseTemplate& operator *=(ElementType value) 
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) * value;

			return *this;
		}

		VectorBaseTemplate& operator *=(const VectorBaseTemplate& value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) * value.at(i);

			return *this;
		}


		VectorBaseTemplate operator *(const VectorBaseTemplate& value) const
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) * value.at(i);

			return vec;
		}


		// Operator /
		friend VectorBaseTemplate operator/(ElementType value, const VectorBaseTemplate& vec)
		{
			VectorBaseTemplate ret;
			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = value / vec.at(i);

			return ret;
		}

		VectorBaseTemplate& operator /=(const VectorBaseTemplate& value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) / value.at(i);

			return *this;
		}

		VectorBaseTemplate& operator /=(const ElementType value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = at(i) / value;

			return *this;
		}


		VectorBaseTemplate operator /(ElementType value) const
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) / value;

			return vec;
		}

		VectorBaseTemplate operator /(const VectorBaseTemplate& value) const
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = at(i) / value.at(i);

			return vec;
		}
		
		// Operator %
		VectorBaseTemplate operator %(const VectorBaseTemplate& value) const
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = ModuloOperator(at(i), value.at(i));
			return vec;
		}

		VectorBaseTemplate operator %(ElementType value) const
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				vec.at(i) = ModuloOperator(at(i), value);
			return vec;
		}

		friend VectorBaseTemplate operator %(ElementType value, const VectorBaseTemplate& vec)
		{
			VectorBaseTemplate ret;
			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = ModuloOperator(value, vec.at(i));
			return ret;
		}

		VectorBaseTemplate& operator %=(ElementType value) 
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = ModuloOperator(at(i), value);
			return this;
		}

		VectorBaseTemplate& operator %=(const VectorBaseTemplate& value)
		{
			for (size_t i = 0; i < DIM; i++)
				at(i) = ModuloOperator(at(i), value.at(i));
			return this;
		}

	


		// mathematic modulo operation
		VectorBaseTemplate Modulo(const VectorBaseTemplate& value)
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				at(i) = Modulo(at(i), value.at(i));
			
			return vec;
		}

		VectorBaseTemplate Modulo(ElementType value)
		{
			VectorBaseTemplate vec;
			for (size_t i = 0; i < DIM; i++)
				at(i) = Modulo(at(i), value);
			return vec;
		}

#pragma endregion

		VectorBaseTemplate Max(const VectorBaseTemplate& rhs) const
		{
			VectorBaseTemplate ret;
			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = (std::max)(at(i), rhs.at(i));

			return ret;
		}

		VectorBaseTemplate Min(const VectorBaseTemplate& rhs) const
		{
			VectorBaseTemplate ret;
			for (size_t i = 0; i < DIM; i++)
				ret.at(i) = (std::min)(at(i), rhs.at(i));

			return ret;
		}


		ElementType Distance(const VectorBaseTemplate& rhs) const
		{
			return (*this - rhs).Length();
		}


		VectorBaseTemplate Lerp(const VectorBaseTemplate& rhs , double t) const
		{
			return (rhs - *this) * t + *this;
		}

		VectorBaseTemplate Slerp(const VectorBaseTemplate& rhs, double t) const
		{
			VectorBaseTemplate a = *this;
			VectorBaseTemplate b = rhs;
			double magA = a.Length();
			double magB = b.Length();
			a /= magA;
			b /= magB;
			double dot = std::clamp(a.Dot(b), -L_One, L_One);
			double theta = acos(dot) * t;
			VectorBaseTemplate relativeVec = (b - a * dot).Normalized();
			VectorBaseTemplate newVec = a * std::cos(theta) + relativeVec * std::sin(theta);
			return newVec * (magA + (magB - magA) * t);
		}


		double Angle(const VectorBaseTemplate& rhs) const
		{
			return std::acos(std::clamp(Dot(rhs) / (Length() * rhs.Length()), -L_One, L_One));
		}

#pragma region Vector2/3 specific

		


		VectorBaseTemplate Cross(const VectorBaseTemplate& rhs) const
		{
			static_assert(DIM == 3, "not a vector with 3 DIMensions");
			
			return
			{
				at(1) * rhs.at(2) - at(2) * rhs.at(1),
				at(2) * rhs.at(0) - at(0) * rhs.at(2),
				at(0) * rhs.at(1) - at(1) * rhs.at(0)
			};
		}

		VectorBaseTemplate Orthogonal() const
		{
			static_assert(DIM == 3, "not a vector with 3 DIMensions");
			return at(2) < at(0) ? VectorBaseTemplate(at(1), -at(0), 0) : VectorBaseTemplate(0, -at(2), at(1));
		}

		

		VectorBaseTemplate Reflect(const VectorBaseTemplate& normal) const
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

		VectorBaseTemplate(const POINT& p)
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
		
		VectorBaseTemplate(const SIZE& p)
		{
            static_assert(DIM == 2, "not a vector with 2 DIMensions");
			at(0) = p.cx;
			at(1) = p.cy;
		}

#endif



#if LMATH_VECTOR_OGRE_EXTENSIONS == 1
        VectorBaseTemplate(const Ogre::Vector2 & v) 
        {
            static_assert(DIM == 2, "not a vector with 2 DIMensions");
            return (v.x, v.y);
        }

        operator Ogre::Vector2() const
        {
            static_assert(DIM == 2, "not a vector with 2 DIMensions");
            return { at(0) ,at(1) };
        }

        VectorBaseTemplate(const Ogre::Vector3& v)
        {
            static_assert(DIM == 3, "not a vector with 3 DIMensions");
            return (v.x, v.y, v.z);
        }

        operator Ogre::Vector3() const
        {
            static_assert(DIM == 3, "not a vector with 3 DIMensions");
            return { at(0) ,at(1) , at(2) };
        }

        VectorBaseTemplate(const Ogre::Vector4& v)
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
	
	template <typename T, size_t DIM >
	class VectorBaseTemplateX : public VectorBaseTemplate<T, DIM>
	{
		using BaseClass = VectorBaseTemplate<T, DIM>;
		using BaseClass::BaseClass;
	};


	template <typename T>
	class VectorBaseTemplateX<T,2> : public VectorBaseTemplate<T, 2>
	{
	public:
		using BaseClass = VectorBaseTemplate<T, 2>;
		using BaseClass::BaseClass;
		using ElementType = typename BaseClass::ElementType;
		VectorBaseTemplateX(const BaseClass& base) : BaseClass(base) {}
		static const VectorBaseTemplateX Down;
		static const VectorBaseTemplateX Up;
		static const VectorBaseTemplateX Left;
		static const VectorBaseTemplateX Right;
	};



	template <typename T>
	class VectorBaseTemplateX<T,3> : public VectorBaseTemplate<T, 3>
	{
	public:
		using BaseClass = VectorBaseTemplate<T, 3>;
		using BaseClass::BaseClass;
		using ElementType = typename BaseClass::ElementType;

		VectorBaseTemplateX(const BaseClass& base) : BaseClass(base) {}
		static const VectorBaseTemplateX Down;
		static const VectorBaseTemplateX Up;
		static const VectorBaseTemplateX Left;
		static const VectorBaseTemplateX Right;
		static const VectorBaseTemplateX Forward;
		static const VectorBaseTemplateX Backward;
	};



	template <typename T>
	class VectorBaseTemplateX<T,4> : public VectorBaseTemplate<T, 4>
	{
	public:
		using BaseClass = VectorBaseTemplate<T, 4>;
		using BaseClass::BaseClass;
		using ElementType = typename BaseClass::ElementType;
		VectorBaseTemplateX(const BaseClass& base) : BaseClass(base) {}
	};



	template <typename T, size_t DIM>
	using VectorBase = VectorBaseTemplateX<T, DIM>;


	
	////Predefined vectors in R3
	template <class T> const VectorBaseTemplateX<T,3> VectorBaseTemplateX<T,3>::Backward(0, 0, -1);
	template <class T> const VectorBaseTemplateX<T,3> VectorBaseTemplateX<T,3>::Forward(0, 0, 1);
	template <class T> const VectorBaseTemplateX<T,3> VectorBaseTemplateX<T,3>::Right(1, 0, 0);
	template <class T> const VectorBaseTemplateX<T,3> VectorBaseTemplateX<T,3>::Left (-1,0,0);
	template <class T> const VectorBaseTemplateX<T,3> VectorBaseTemplateX<T,3>::Up(0, 1, 0);
	template <class T> const VectorBaseTemplateX<T,3> VectorBaseTemplateX<T,3>::Down(0,-1,0);

	//Predefined vectors in R2
	template <class T> const VectorBaseTemplateX<T,2> VectorBaseTemplateX<T,2>::Down(0, -1);
	template <class T> const VectorBaseTemplateX<T,2> VectorBaseTemplateX<T,2>::Up(0, 1);
	template <class T> const VectorBaseTemplateX<T,2> VectorBaseTemplateX<T,2>::Left(-1, 0);
	template <class T> const VectorBaseTemplateX<T,2> VectorBaseTemplateX<T,2>::Right (1,0);
	

	

	template <class T, size_t DIM>
	const VectorBaseTemplate<T, DIM> VectorBaseTemplate<T, DIM>::Zero(VectorBaseTemplate<T, DIM>::L_Zero);

	template <class T, size_t DIM>
	const VectorBaseTemplate<T, DIM> VectorBaseTemplate<T, DIM>::Unit(VectorBaseTemplate<T, DIM>::L_One);


}


#endif