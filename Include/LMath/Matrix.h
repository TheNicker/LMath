#ifndef __LMATH_MATRIX__
#define __LMATH_MATRIX__

#include "Vector.h"
#include "Quaternion.h"

namespace LMath
{
	//Row major base matrix 
	template <typename T,size_t ROWS, size_t COLS>
	class MatrixBase
	{
	public:
		static bool constexpr  IsSquareMatrix() {return ROWS == COLS;}
		using VectorType =  VectorBase<T, COLS>;
		using ElementType = typename VectorType::ElementType;

		static const MatrixBase Zero;
		static const MatrixBase Identity;
		#pragma region Constructors

		
		template <typename... Args>
		MatrixBase(ElementType first, Args... args)
		{
			_Assign<0, ROWS * COLS>(first, args...);
		}

		MatrixBase() {}
		MatrixBase(std::array<VectorType,ROWS>  vectors) :mVectors(vectors) {}
		MatrixBase(VectorType  vectors[ROWS]) :mVectors(vectors) {}
		MatrixBase( std::array<typename VectorType::ElementType, ROWS * COLS> elements)
		{
			std::memcpy(mVectors, elements.data(), sizeof(ElementType) * ROWS * COLS);
		}

#pragma endregion
		//Operator *
		template <size_t RHS_ROWS, size_t RHS_COLS>
		MatrixBase <ElementType, ROWS, RHS_COLS> operator*(const MatrixBase<ElementType, RHS_ROWS, RHS_COLS>& rhs)
		{
			using MatrixOutputType = MatrixBase <ElementType, ROWS, RHS_COLS>;
			MatrixOutputType result = MatrixOutputType::Zero;
			for (size_t col = 0; col < RHS_COLS ; col++)
				for (size_t row = 0; row < ROWS ; row++)
					for (size_t idx = 0; idx < ROWS ; idx++)
						result.at(row).at(col) += this->at(row).at(idx) * rhs.at(idx).at(col);

			return result;
		}


		
		MatrixBase operator*(ElementType value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++) 
					result.at(row).at(col) = at(row).at(col) * value;

			return result;
		}

		VectorType operator*(const VectorType& value) const
		{
			VectorType  result = VectorType::Zero;
			for (size_t row = 0; row < ROWS; row++)
				for (size_t col = 0; col < COLS; col++)
					result.at(row) += at(row).at(col) * value.at(col);
					

			return result;
		}



		//operator +

		MatrixBase operator+(ElementType value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row).at(col) = at(row).at(col) + value;

			return result;
		}

		MatrixBase operator+(MatrixBase value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row).at(col) = at(row).at(col) + value.at(row).at(col);

			return result;
		}

		//operator - 
		MatrixBase operator-(ElementType value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row).at(col) = at(row).at(col) - value;

			return result;
		}

		MatrixBase operator-(MatrixBase value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row).at(col) = at(row).at(col) - value.at(row).at(col);

			return result;
		}

		//operator /

		MatrixBase operator/(ElementType value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row).at(col) = at(row).at(col) / value;

			return result;
		}


		bool operator ==(const MatrixBase& value) const
		{
			return std::memcmp(this, &value, sizeof(value)) == 0;
		}

		bool operator !=(const MatrixBase& value) const
		{
			return !(*this == value);
		}


#pragma region accessors
		VectorType& at(size_t idx)
		{
			return mVectors[idx];
		}

		const VectorType& at(size_t idx) const
		{
			return mVectors[idx];
		}

		ElementType& at(size_t row,size_t col)
		{
			return mVectors[row].at(col);
		}

		const ElementType& at(size_t row, size_t col) const
		{
			return mVectors[row].at(col);
		}
#pragma endregion

		ElementType Determinant() const
		{
			static_assert(IsSquareMatrix() , "Error, Determinant is only applicable to square matrices.");
			return _ComputeDeterminant(*this);
		}

		MatrixBase<T, COLS, ROWS>  Transpose()
		{
			MatrixBase<T, COLS, ROWS> result;

			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(col).at(row) = this->at(row).at(col);

			return result;
		}


		MatrixBase Scale(const MatrixBase& value)
		{
			static_assert( IsSquareMatrix(), "Scale applies only to square matrices.");

			MatrixBase result = MatrixBase::Zero;
			for (size_t pos = 0; pos < ROWS; pos++)
				result.at(pos).at(pos) = this->at(pos).at(pos) * value.at(pos).at(pos);

			return result;
		}



		MatrixBase<T, ROWS - 1, COLS - 1> SubMartix(size_t rowToRemove, size_t coltoRemove) const
		{
			MatrixBase<T, ROWS - 1, COLS - 1> result;
			size_t currentCol = 0;
			size_t currentRow = 0;
			for (size_t col = 0; col < COLS; col++)
			{
				if (col == coltoRemove)
					continue;
				for (size_t row = 0; row < ROWS; row++)
				{
					if (row == rowToRemove)
						continue;
					result.at(currentRow).at(currentCol) = at(row).at(col);
					currentRow++;
				}
				currentCol++;
				currentRow = 0;
			}

			return result;
			
		}

		
		bool IsInvertible() const
		{
			if constexpr (IsSquareMatrix() == false)
			{
				return false;
			}
			else
			{
				return std::abs(Determinant()) > VectorType::L_Zero;
			}
		}

		MatrixBase Inverse() const
		{
			static_assert(IsSquareMatrix(), "Scale applies only to square matrices.");

			MatrixBase minors;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					minors.at(row).at(col) = _ComputeDeterminant(SubMartix(row, col));

			return minors / minors.Determinant();
		}

		
		
		QuaternionBase<T> ToQuaternion() const
		{
			static_assert(ROWS == 3 && COLS == 3, "Only Matrix 3X3 can be converted to quanternoin.");
			

			
			auto  Add = [&](size_t idx1, size_t idx2) -> ElementType
			{
				return at(idx1, idx2) + at(idx2, idx1);
			};

			auto Sub = [&](size_t idx1, size_t idx2)-> ElementType
			{
				return at(idx1, idx2) - at(idx2, idx1);
			};

			
			QuaternionBase<T> q;
			double trace =  at(0).at(0) + at(1).at(1) + at(2).at(2);
			if (trace > 0)
			{
				double s = 0.5 / std::sqrt(trace + 1);
				q.W() = 0.25 / s;
				q.X() = Sub(2,1) * s;
				q.Y() = Sub(0,2) * s;
				q.Z() = Sub(1,0) * s;
			}
			else
			{
				if (at(0,0)  > at(1,1) && at(0,0) > at(2,2))
				{
					double s = 2 * sqrt(1 + at(0).at(0) - at(1).at(1) - at(2).at(2));
					q.W() = Sub(2,1) / s;
					q.X() = 0.25 * s;
					q.Y() = Add(0,1) / s;
					q.Z() = Add(0,2) / s;
				}
				else if (at(1).at(1) > at(2).at(2))
				{
					double s = 2 * sqrt(1 + at(1).at(1) - at(0).at(0) - at(2).at(2));
					q.W() = Sub(0,2) / s;
					q.X() = Add(0,1) / s;
					q.Y() = 0.25 * s;
					q.Z() = Add(1,2) / s;
				}
				else
				{
					double s = 2 * sqrt(1 + at(2).at(2) - at(0).at(0) - at(1).at(1));
					q.W() = Sub(1,0) / s;
					q.X() = Add(0,2) / s;
					q.Y() = Add(1,2) / s;
					q.Z() = 0.25 * s;
				}
			}
			return q;
		}

		//public class methods
		static MatrixBase  CreateScaleMatrix(const VectorType& scaleVec)
		{
			for (size_t pos = 0; i < ROWS; pos++)
				at(pos).at(pos) = scaleVec.at(pos);
			
		}

		template <typename... Args>
		static MatrixBase  CreateScaleMatrix(ElementType first, Args... args)
		{
			MatrixBase scaleMatrix = MatrixBase::Zero;
			scaleMatrix._AssignScale<0, ROWS> (first, args...);
			return scaleMatrix;
		}

		static MatrixBase<T,3,3> FromQuaternion(const QuaternionBase<T>& rotation)
		{
			MatrixBase<T, 3,3> m = MatrixBase<T, 3,3>::Zero;
			
			double sqw = rotation.W() * rotation.W();
			double sqx = rotation.X() * rotation.X();
			double sqy = rotation.Y() * rotation.Y();
			double sqz = rotation.Z() * rotation.Z();

			double invSqr = 1 / (sqx + sqy + sqz + sqw);
			m.at(0).at(0) = (sqx - sqy - sqz + sqw) * invSqr;
			m.at(1).at(1) = (-sqx + sqy - sqz + sqw) * invSqr;
			m.at(2).at(2) = (-sqx - sqy + sqz + sqw) * invSqr;

			double tmp1 = rotation.X() * rotation.Y();
			double tmp2 = rotation.Z() * rotation.W();
			m.at(1).at(0) = 2.0 * (tmp1 + tmp2) * invSqr;
			m.at(0).at(1) = 2.0 * (tmp1 - tmp2) * invSqr;

			tmp1 = rotation.X() * rotation.Z();
			tmp2 = rotation.Y() * rotation.W();

			m.at(2,0)  = 2.0 * (tmp1 - tmp2) * invSqr;
			m.at(0,2)  = 2.0 * (tmp1 + tmp2) * invSqr;
			tmp1 = rotation.Y() * rotation.Z();
			tmp2 = rotation.X() * rotation.W();
			m.at(2,1) = 2.0 * (tmp1 + tmp2) * invSqr;
			m.at(1,2) = 2.0 * (tmp1 - tmp2) * invSqr;
			return m;
		}
		

		


	private: //Private helper methods
		
		

		static MatrixBase  CreateIdentityMatrix()
		{
			MatrixBase identity = MatrixBase::Zero;
			
			static_assert(IsSquareMatrix(), "Error, identity matrix applie only to square matrices");
			
			for (size_t row = 0; row < ROWS; row++)
				identity.at(row).at(row) = VectorType::L_One;
			return identity;

		}

		static MatrixBase  CreateZeroMatrix()
		{
			MatrixBase zero;
			for (size_t i = 0; i < ROWS; i++)
				zero.mVectors[i] = VectorType(VectorType::L_Zero);
		
			return zero;
		}

		template <size_t INDEX,size_t NUM_ELEMENTS>
		void _Assign(ElementType element)
		{

			static_assert(INDEX == NUM_ELEMENTS - 1, "Error, wrong number of arguments passed to VectorBase constructor");
			at(INDEX / COLS).at(INDEX % COLS) = element;
		}

		template <size_t INDEX, size_t NUM_ELEMENTS ,typename ...ARGS>
		void _Assign(ElementType element, ARGS... args)
		{
			at(INDEX / COLS).at(INDEX % COLS) = element;
			_Assign<INDEX + 1, NUM_ELEMENTS>(args...);
		}


		template <size_t POS, size_t NUM_ELEMENTS>
		void _AssignScale(ElementType element)
		{

			static_assert(POS == NUM_ELEMENTS -1, "Error, wrong number of arguments passed to VectorBase constructor");
			at(POS).at(POS) = element;
		}

		template <size_t POS, size_t NUM_ELEMENTS , typename ...ARGS>
		void _AssignScale(ElementType element, ARGS... args)
		{
			at(POS).at(POS) = element;
			_Assign<POS + 1, NUM_ELEMENTS>(args...);
		}


		template <size_t MATRIX_SIZE>
		static ElementType _ComputeDeterminant(const MatrixBase <T, MATRIX_SIZE, MATRIX_SIZE>& matrix)
		{
			size_t ARBITRARY_ROW = 0;
			if constexpr (MATRIX_SIZE == 2)
				return matrix.at(0).at(0) * matrix.at(1).at(1) - matrix.at(0).at(1) * matrix.at(1).at(0);
			else
			{
				ElementType determinant = 0;
				for (size_t col = 0; col < MATRIX_SIZE; col++)
				{

					ElementType sign = col % 2 == 0 ? 1 : -1;
					determinant += sign * matrix.at(ARBITRARY_ROW).at(col) * _ComputeDeterminant(matrix.SubMartix(ARBITRARY_ROW, col));
				}

				return determinant;
			}
		}



	private:
		VectorType mVectors[ROWS];
	};

	template <typename T, size_t rows,size_t cols>
	const MatrixBase<T, rows, cols> MatrixBase<T, rows, cols>::Zero = MatrixBase<T, rows, cols>::CreateZeroMatrix();
	
	template <typename T, size_t rows, size_t cols>
	const MatrixBase<T, rows, cols> MatrixBase<T, rows, cols>::Identity = MatrixBase<T, rows, cols>::CreateIdentityMatrix();
}

#endif