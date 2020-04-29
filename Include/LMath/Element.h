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

#ifndef __LMATH_ELEMENT_H__
#define __LMATH_ELEMENT_H__ 1

#include <type_traits>
#include <algorithm>
#include <limits>
#include <cmath>
#include <sstream>



namespace LMath
{
	template <typename T, size_t dim>
	class ElementBase
	{

	public:
		using ElementType = T;
		inline static constexpr size_t DIM = dim;
		inline static ElementType L_Epsilon = std::numeric_limits<ElementType>::epsilon();
		inline static constexpr ElementType L_Zero = static_cast<ElementType>(0);
		inline static constexpr ElementType L_Half = static_cast<ElementType>(0.5);
		inline static constexpr ElementType L_One = static_cast<ElementType>(1);
		inline static constexpr ElementType L_Two = static_cast<ElementType>(2);
		inline static constexpr ElementType L_Four = static_cast<ElementType>(4);

		

		T& at(size_t idx)
		{
			return mElements[idx];
		}

		const T& at(size_t idx) const
		{
			return mElements[idx];
		}

		template <size_t INDEX = 0, typename... Args>
		ElementBase(ElementType first, Args... args)
		{
			_Assign<INDEX>(static_cast<ElementType>(first), args...);
		}


		template <size_t INDEX = 0, size_t RHS_DIM, typename... Args>
		ElementBase(const ElementBase<ElementType, RHS_DIM>& rhs, Args... args)
		{
			_Assign<INDEX>(rhs, args...);
		}


		ElementBase(const ElementType* data)
		{
			memcpy(&mElements, data, sizeof(ElementType) * dim);
		}

		ElementBase(const std::array<ElementType,dim>& arr)
		{
			memcpy(&mElements, arr.data(), sizeof(ElementType) * dim);
		}


		// Convert from the same type different dimension.
		template <size_t RHS_DIM >
		explicit ElementBase(const ElementBase<ElementType, RHS_DIM>& rhs)
		{
			_AssignFromElementBase(rhs);
		}

		// Convert from different type and dimension.
		template <typename RHS_U, size_t RHS_DIM >
		explicit ElementBase(const ElementBase<RHS_U, RHS_DIM>& rhs)
		{
			_AssignFromElementBase(rhs);
		}
		
		

		ElementBase() 
#if LMATH_ALWAYS_INITIALIZE_VARIABLES == 1
			: ElementBase(L_Zero)
#endif
		{

		}

		ElementBase(ElementType initializationValue)
		{
			for (size_t i = 0; i < dim; i++)
				at(i) = initializationValue;
		}


		std::string ToString()
		{
			std::stringstream stream;
			stream << '(' << at(0);

			for (size_t i = 1; i < dim; i++)
				stream << ", " << at(i);

			stream << ')';

			return stream.str();
		}



	private:
#pragma region Private helper methods
	private:
		template <size_t INDEX>
		void _Assign(ElementType element)
		{

			static_assert(INDEX == dim - 1, "Error, wrong number of arguments passed to VectorBase constructor");
			at(INDEX) = element;
		}

		template <size_t INDEX, typename ...ARGS>
		void _Assign(ElementType element, ARGS... args)
		{
			at(INDEX) = element;
			_Assign<INDEX + 1>(args...);
		}

		template <size_t INDEX, size_t RHS_DIM>
		void _Assign(const ElementBase<ElementType, RHS_DIM>& rhs)
		{
			static_assert(INDEX + RHS_DIM == dim, "Error, wrong number of arguments passed to VectorBase constructor");
			for (size_t i = 0; i < RHS_DIM; i++)
				at(i + INDEX) = rhs.at(i);
		}

		template <size_t INDEX, size_t RHS_DIM, typename ...ARGS>
		void _Assign(const ElementBase<ElementType, RHS_DIM>& rhs, ARGS... args)
		{
			for (size_t i = 0; i < RHS_DIM; i++)
				at(i + INDEX) = rhs.at(i);

			_Assign<INDEX + RHS_DIM>(args...);
		}

		template <typename RHS_TYPE, size_t RHS_DIM>
		void _AssignFromElementBase(const ElementBase<RHS_TYPE, RHS_DIM>& rhs)
		{
			const size_t elementsToConvert = (std::min)(dim, RHS_DIM);

			for (size_t i = 0; i < elementsToConvert; i++)
				at(i) = rhs.at(i);


			for (size_t i = elementsToConvert; i < dim; i++)
				at(i) = 0;

		}

#pragma endregion



	private:
		T mElements[dim];
	};
}
#endif