#ifndef __LMATH_MATRIX__
#define __LMATH_MATRIX__

#include "Vector.h"
#include "Quaternion.h"

namespace LMath
{
	enum class MatrixVectors
	{
		  Column
		, Row
	};

	//Row major base matrix 
	template <typename T,size_t ROWS, size_t COLS, MatrixVectors VECTORS = MatrixVectors::Row>
	class MatrixBase
	{
	public:
		static bool constexpr  IsSquareMatrix() { return ROWS == COLS; }
		using VectorType = VectorBase<T, COLS>;
		using ElementType = typename VectorType::ElementType;
		using Base = MatrixBase<T, ROWS, COLS, VECTORS>;
		inline static constexpr size_t Rows = ROWS;
		inline static constexpr size_t Cols = COLS;


		static const MatrixBase Zero;
		static const MatrixBase Identity;

#pragma region Constructors

		template <typename... Args>
		constexpr MatrixBase(ElementType first, Args... args)
		{
			_Assign<0, ROWS* COLS>(first, args...);
		}

		constexpr MatrixBase()
		{
#if LMATH_ALWAYS_INITIALIZE_VARIABLES == 1
			* this = Matrix::Zero;
#endif
			/// default constructor - do nothing
			/// 
		}

#pragma endregion

		// Convert matrix to any matrix
		template <typename U, size_t RHS_ROWS, size_t RHS_COLS>
		explicit MatrixBase(const MatrixBase<U, RHS_ROWS, RHS_COLS, VECTORS>& rhs)
		{

			const size_t rowsToConvert = (std::min)(ROWS, RHS_ROWS);
			const size_t colsToConvert = (std::min)(COLS, RHS_COLS);
			*this = Zero;

			for (size_t row = 0; row < rowsToConvert; row++)
				for (size_t col = 0; col < colsToConvert; col++)
					at(row, col) = static_cast<ElementType>(rhs.at(row, col));
		}

		ElementType* data() { return &mElements[0]; }
		const ElementType* data() const { return &mElements[0]; }

		//Operator *
		template <size_t RHS_ROWS, size_t RHS_COLS>
		MatrixBase <ElementType, ROWS, RHS_COLS, VECTORS> operator*(const MatrixBase<ElementType, RHS_ROWS, RHS_COLS, VECTORS>& rhs) const
		{
			using MatrixOutputType = MatrixBase <ElementType, ROWS, RHS_COLS, VECTORS>;
			MatrixOutputType result = MatrixOutputType::Zero;

			for (size_t col = 0; col < RHS_COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					for (size_t idx = 0; idx < ROWS; idx++)
					{
						if constexpr (VECTORS == MatrixVectors::Row)
						{
							result.at(row, col) += this->at(idx, col) * rhs.at(row, idx);
						}
						else
						{
							result.at(row, col) += this->at(row, idx) * rhs.at(idx, col);
						}
					}
			return result;
		}



		MatrixBase operator*(ElementType value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row, col) = at(row, col) * value;

			return result;
		}

		VectorType operator*(const VectorType& value) const
		{
			VectorType  result = VectorType::Zero;

			if constexpr (VECTORS == MatrixVectors::Row)
			{

				for (size_t col = 0; col < COLS; col++)
					for (size_t row = 0; row < ROWS; row++)
						result.at(col) += at(row, col) * value.at(row);
			}
			else if constexpr (VECTORS == MatrixVectors::Column)
			{
				for (size_t row = 0; row < ROWS; row++)
					for (size_t col = 0; col < COLS; col++)
						result.at(row) += at(row, col) * value.at(col);
			}

			return result;
		}

#if LMATH_ENABLE_MATRIX4_MUL_IN_VECTOR3 == 1
		////special optional common case when multiplying 4X4 matrices with vector3, (add 1.0 to the last componenet).
		template <typename Vector3 = VectorBase<ElementType, 3>>
		Vector3 operator*(const Vector3& valueVec3) const
		{
			static_assert(ROWS == 4 && COLS == 4, "Special case only for multiplying Matrix4 in vector 3 ");
			VectorType result = *this * VectorType(valueVec3, VectorType::L_One);
			return static_cast<Vector3>(result / result.at(3));
		}
#endif



		//operator +

		MatrixBase operator+(ElementType value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row, col) = at(row, col) + value;

			return result;
		}

		MatrixBase operator+(MatrixBase value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row, col) = at(row, col) + value.at(row, col);

			return result;
		}

		//operator - 
		MatrixBase operator-(ElementType value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row, col) = at(row, col) - value;

			return result;
		}

		MatrixBase operator-()  const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row, col) = -this->at(row, col);

			return result;
		}


		MatrixBase operator-(MatrixBase value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row, col) = at(row, col) - value.at(row, col);

			return result;
		}

		//operator /

		MatrixBase operator/(ElementType value) const
		{
			MatrixBase result;
			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(row, col) = at(row, col) / value;

			return result;
		}


		bool operator ==(const MatrixBase& value) const
		{
			for (size_t i = 0; i < ROWS * COLS; i++)
				if (mElements[i] != value.mElements[i])
					return false;

			return true;
		}

		bool operator !=(const MatrixBase& value) const
		{
			return !(*this == value);
		}


#pragma region accessors
		ElementType& at(size_t idx)
		{
			return  mElements[idx];
		}

		const ElementType& at(size_t idx) const
		{
			return mElements[idx];
		}

		ElementType& at(size_t row, size_t col)
		{
			return mMatrix[row][col];
		}

		const ElementType& at(size_t row, size_t col) const
		{
			return mMatrix[row][col];
		}

		void assignVector(const VectorType& vector, size_t row)
		{
			if constexpr (VECTORS == MatrixVectors::Row)
			{
				for (size_t col = 0; col < COLS; col++)
					at(row, col) = vector.at(col);
			}
			else if constexpr (VECTORS == MatrixVectors::Column)
			{
				for (size_t col = 0; col < COLS; col++)
					at(col, row) = vector.at(col);
			}
			else
			{
				static_assert(false, "Invalid state");
			}
		}

#pragma endregion

		ElementType Determinant() const
		{
			static_assert(IsSquareMatrix(), "Error, Determinant is only applicable to square matrices.");
			return _ComputeDeterminant(*this);
		}

		MatrixBase<T, COLS, ROWS, VECTORS>  Transpose() const
		{
			MatrixBase<T, COLS, ROWS, VECTORS> result;

			for (size_t col = 0; col < COLS; col++)
				for (size_t row = 0; row < ROWS; row++)
					result.at(col, row) = this->at(row, col);

			return result;
		}


		MatrixBase MultiplyPerComponent(const MatrixBase& value)
		{
			MatrixBase result = *this;

			for (size_t row = 0; row < ROWS; row++)
				for (size_t col = 0; col < COLS; col++)
					result.at(row, col) *= value.at(row, col);

			return result;
		}



		MatrixBase<T, ROWS - 1, COLS - 1, VECTORS> SubMartix(size_t rowToRemove, size_t coltoRemove) const
		{
			MatrixBase<T, ROWS - 1, COLS - 1, VECTORS> result;
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
					result.at(currentRow, currentCol) = at(row, col);
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


		MatrixBase Adjoint() const
		{
			static_assert(IsSquareMatrix(), "Adjoint applies only to square matrices.");

			MatrixBase minors;

			for (size_t row = 0; row < ROWS; row++)
				for (size_t col = 0; col < COLS; col++)
					minors.at(row, col) = GetMinorSign(row, col) * _ComputeDeterminant(SubMartix(row, col));

			return minors;
		}

		MatrixBase Inverse() const
		{
			static_assert(IsSquareMatrix(), "Inverse applies only to square matrices.");
			return Adjoint().Transpose() / Determinant();
		}

		QuaternionBase<T> ToQuaternion() const
		{
			static_assert(ROWS == 3 && COLS == 3, "Only Matrix 3X3 can be converted to quaternion.");



			auto  Add = [&](size_t idx1, size_t idx2) -> ElementType
			{
				return at(idx1, idx2) + at(idx2, idx1);
			};

			auto Sub = [&](size_t idx1, size_t idx2)-> ElementType
			{
				return at(idx1, idx2) - at(idx2, idx1);
			};


			QuaternionBase<T> q;
			double trace = at(0, 0) + at(1, 1) + at(2, 2);
			if (trace > 0)
			{
				double s = 0.5 / std::sqrt(trace + 1);
				q.W() = 0.25 / s;
				q.X() = Sub(2, 1) * s;
				q.Y() = Sub(0, 2) * s;
				q.Z() = Sub(1, 0) * s;
			}
			else
			{
				if (at(0, 0) > at(1, 1) && at(0, 0) > at(2, 2))
				{
					double s = 2 * sqrt(1 + at(0, 0) - at(1, 1) - at(2, 2));
					q.W() = Sub(2, 1) / s;
					q.X() = 0.25 * s;
					q.Y() = Add(0, 1) / s;
					q.Z() = Add(0, 2) / s;
				}
				else if (at(1, 1) > at(2, 2))
				{
					double s = 2 * sqrt(1 + at(1, 1) - at(0, 0) - at(2, 2));
					q.W() = Sub(0, 2) / s;
					q.X() = Add(0, 1) / s;
					q.Y() = 0.25 * s;
					q.Z() = Add(1, 2) / s;
				}
				else
				{
					double s = 2 * sqrt(1 + at(2, 2) - at(0, 0) - at(1, 1));
					q.W() = Sub(1, 0) / s;
					q.X() = Add(0, 2) / s;
					q.Y() = Add(1, 2) / s;
					q.Z() = 0.25 * s;
				}
			}
			return q;
		}



		void SetScale(const VectorType& scaleVec)
		{
			for (size_t pos = 0; pos < ROWS; ++pos)
				at(pos, pos) = scaleVec.at(pos);
		}

		template <typename... Args>
		void  SetScale(ElementType first, Args... args)
		{
			_AssignScale<0, ROWS>(first, args...);
		}

		static MatrixBase  CreateScaleMatrix(const VectorType& scaleVec)
		{
			MatrixBase scaleMatrix = MatrixBase::Identity;
			scaleMatrix.SetScale(scaleVec);
			return scaleMatrix;
		}

		template <typename... Args>
		static MatrixBase  CreateScaleMatrix(ElementType first, Args... args)
		{
			MatrixBase scaleMatrix = MatrixBase::Identity;
			scaleMatrix.SetScale(first, args...);
			return scaleMatrix;
		}

		static MatrixBase<ElementType, 4, 4, VECTORS>  CreateTranslationMatrix(const VectorBase<ElementType, 3>& trans)
		{
			MatrixBase<ElementType, 4, 4, VECTORS> transMatrix = MatrixBase<ElementType, 4, 4>::Identity;
			transMatrix.at(3, 0) = trans.at(0);
			transMatrix.at(3, 1) = trans.at(1);
			transMatrix.at(3, 2) = trans.at(2);
			return transMatrix;
		}


		static MatrixBase FromQuaternion(const QuaternionBase<T>& rotation)
		{
			static_assert(IsSquareMatrix(), " Matrix dimensions must be square");
			static_assert(ROWS >= 3, "Matrix dimensions must be equal or greater to 3");

			MatrixBase rotationMatrix = static_cast<MatrixBase>(Create3X3RotationMatrix(rotation));

			for (size_t pos = 3; pos < ROWS; pos++)
				rotationMatrix.at(pos, pos) = VectorType::L_One;

			return rotationMatrix;
		}


		static MatrixBase<ElementType, 3, 3, VECTORS > Create3X3RotationMatrix(const QuaternionBase<T>& rotation)
		{
#if LMATH_ENABLE_AUTO_VECTOR_NORMALIZATION == 0
			const ElementType fTx = rotation.X() + rotation.X();
			const ElementType fTy = rotation.Y() + rotation.Y();
			const ElementType fTz = rotation.Z() + rotation.Z();
			const ElementType fTwx = fTx * rotation.W();
			const ElementType fTwy = fTy * rotation.W();
			const ElementType fTwz = fTz * rotation.W();
			const ElementType fTxx = fTx * rotation.X();
			const ElementType fTxy = fTy * rotation.X();
			const ElementType fTxz = fTz * rotation.X();
			const ElementType fTyy = fTy * rotation.Y();
			const ElementType fTyz = fTz * rotation.Y();
			const ElementType fTzz = fTz * rotation.Z();
			if constexpr (VECTORS == MatrixVectors::Row)
			{
				return
				{
					// Left
					  VectorType::L_One - (fTyy + fTzz)
					, fTxy + fTwz
					, fTxz - fTwy

					// Up
					, fTxy - fTwz
					, VectorType::L_One - (fTxx + fTzz)
					, fTyz + fTwx

					//Forward
					, fTxz + fTwy
					, fTyz - fTwx
					, VectorType::L_One - (fTxx + fTyy)

				};
			}
			else if constexpr (VECTORS == MatrixVectors::Column)
			{
				return
				{
					  VectorType::L_One - (fTyy + fTzz)
					, fTxy - fTwz
					, fTxz + fTwy
					, fTxy + fTwz
					, VectorType::L_One - (fTxx + fTzz)
					, fTyz - fTwx
					, fTxz - fTwy
					, fTyz + fTwx
					, VectorType::L_One - (fTxx + fTyy)
				};
			}
#else
			MatrixBase<T, 3, 3, VECTORS> m = MatrixBase<T, 3, 3, VECTORS>::Zero;
			double sqw = rotation.W() * rotation.W();
			double sqx = rotation.X() * rotation.X();
			double sqy = rotation.Y() * rotation.Y();
			double sqz = rotation.Z() * rotation.Z();
			double invSqr = 1 / (sqx + sqy + sqz + sqw);
			m.at(0, 0) = (sqx - sqy - sqz + sqw) * invSqr;
			m.at(1, 1) = (-sqx + sqy - sqz + sqw) * invSqr;
			m.at(2, 2) = (-sqx - sqy + sqz + sqw) * invSqr;



			double tmp1 = rotation.X() * rotation.Y();
			double tmp2 = rotation.Z() * rotation.W();


			if constexpr (VECTORS == MatrixVectors::Row)
			{
				m.at(0, 1) = 2.0 * (tmp1 + tmp2) * invSqr;
				m.at(1, 0) = 2.0 * (tmp1 - tmp2) * invSqr;
			}
			else
			{
				m.at(1, 0) = 2.0 * (tmp1 + tmp2) * invSqr;
				m.at(0, 1) = 2.0 * (tmp1 - tmp2) * invSqr;
			}


			tmp1 = rotation.X() * rotation.Z();
			tmp2 = rotation.Y() * rotation.W();

			if constexpr (VECTORS == MatrixVectors::Row)
			{
				m.at(0, 2) = 2.0 * (tmp1 - tmp2) * invSqr;
				m.at(2, 0) = 2.0 * (tmp1 + tmp2) * invSqr;
			}
			else
			{
				m.at(2, 0) = 2.0 * (tmp1 - tmp2) * invSqr;
				m.at(0, 2) = 2.0 * (tmp1 + tmp2) * invSqr;
			}
			
			tmp1 = rotation.Y() * rotation.Z();
			tmp2 = rotation.X() * rotation.W();

			if constexpr (VECTORS == MatrixVectors::Row)
			{
				m.at(1, 2) = 2.0 * (tmp1 + tmp2) * invSqr;
				m.at(2, 1) = 2.0 * (tmp1 - tmp2) * invSqr;
			}
			else
			{
				m.at(2, 1) = 2.0 * (tmp1 + tmp2) * invSqr;
				m.at(1, 2) = 2.0 * (tmp1 - tmp2) * invSqr;
			}

			return m;
#endif
		}
		
		static MatrixBase<T, 4, 4, VECTORS>  CreateViewMatrix(const VectorBase<T, 3>& position, const QuaternionBase<T>& orientation)
		{
			using Matrix4 = MatrixBase<T, 4, 4, VECTORS>;
			using Matrix3 = MatrixBase<T, 3, 3, VECTORS>;

			// Column major order view matrix:
			//  [ Lx  Ly  Lz  0   ]
			//  [ Ux  Uy  Uz  0   ]
			//  [ Dx  Dy  Dz  0   ]
			//  [ Tx  Ty  Tz  1   ]

			Matrix3 rotation = Matrix3::FromQuaternion(orientation);
			Matrix3 transpose = rotation.Transpose();
			Matrix4 viewMatrix = static_cast<Matrix4>(transpose);
			
			viewMatrix.assignVector({ -transpose * position, 1.0 }, 3);
			return viewMatrix;

		}

		static MatrixBase<T, 4, 4, VECTORS>  CreateProjectionMatrix(
			ElementType left
			, ElementType bottom
			, ElementType top
			, ElementType right
			, ElementType n
			, ElementType f
			, ElementType z_clip
			, bool leftHanded
		)
		{

			MatrixBase<T, 4, 4, VECTORS> projectionMatrix = MatrixBase<T, 4, 4>::Zero;

			auto inv_width = static_cast<ElementType>(1) / (right - left);
			auto inv_height = static_cast<ElementType>(1) / (top - bottom);
			auto inv_depth = static_cast<ElementType>(1) / (f - n);
			auto near2 = static_cast<ElementType>(2) * n;
			auto s = leftHanded == true ? static_cast<ElementType>(1) : static_cast<ElementType>(-1);

			if (z_clip == -1) 
			{
				projectionMatrix.at(2, 2) = s * (f + n) * inv_depth;
				projectionMatrix.at(3, 2) = static_cast<ElementType>(-2) * f * n * inv_depth;
			}
			else 
			{ // z_clip == z_clip_zero
				projectionMatrix.at(2, 2) = s * f * inv_depth;
				projectionMatrix.at(3, 2) = -s * n * projectionMatrix.at(2, 2);
			}

			projectionMatrix.at(0, 0) = near2 * inv_width;
			projectionMatrix.at(1, 1) = near2 * inv_height;
			projectionMatrix.at(2, 0) = -s * (right + left) * inv_width;
			projectionMatrix.at(2, 1) = -s * (top + bottom) * inv_height;
			projectionMatrix.at(2, 3) = s;
			 
			return  projectionMatrix;
		}

		


	private: //Private helper methods
		
		static ElementType GetMinorSign(size_t row, size_t col)
		{
			return ((row % 2) == 0 ? 1 : -1) * ((col % 2) == 0 ? 1 : -1);
		}
		

		constexpr static MatrixBase  CreateIdentityMatrix()
		{
			MatrixBase identity = MatrixBase::Zero;
			
			static_assert(IsSquareMatrix(), "Error, identity matrix applie only to square matrices");
			
			for (size_t row = 0; row < ROWS; row++)
				identity.at(row, row) = VectorType::L_One;
			return identity;

		}

		constexpr static MatrixBase  CreateZeroMatrix()
		{
			MatrixBase zero;
			for (size_t i = 0; i < ROWS * COLS; i++)
				zero.mElements[i] = VectorType::L_Zero;
		
			return zero;
		}

		template <size_t INDEX,size_t NUM_ELEMENTS>
		constexpr void _Assign(ElementType element)
		{

			static_assert(INDEX == NUM_ELEMENTS - 1, "Error, wrong number of arguments passed to VectorBase constructor");
			
			at(INDEX) = element;
		}

		template <size_t INDEX, size_t NUM_ELEMENTS ,typename ...ARGS>
		constexpr void _Assign(ElementType element, ARGS... args)
		{
			at(INDEX) = element;
			_Assign<INDEX + 1, NUM_ELEMENTS>(args...);
		}


		template <size_t POS, size_t NUM_ELEMENTS>
		constexpr void _AssignScale(ElementType element)
		{

			static_assert(POS == NUM_ELEMENTS -1, "Error, wrong number of arguments passed to VectorBase constructor");
			at(POS, POS) = element;
		}

		template <size_t POS, size_t NUM_ELEMENTS , typename ...ARGS>
		constexpr void _AssignScale(ElementType element, ARGS... args)
		{
			at(POS, POS) = element;
			_AssignScale<POS + 1, NUM_ELEMENTS>(args...);
		}


		template <size_t MATRIX_SIZE, MatrixVectors VECTORS>
		static ElementType _ComputeDeterminant(const MatrixBase <T, MATRIX_SIZE, MATRIX_SIZE, VECTORS>& matrix)
		{
			size_t ARBITRARY_ROW = 0;
			if constexpr (MATRIX_SIZE == 2)
				return matrix.at(0, 0) * matrix.at(1, 1) - matrix.at(0, 1) * matrix.at(1, 0);
			else
			{
				ElementType determinant = 0;
				for (size_t col = 0; col < MATRIX_SIZE; col++)
					determinant += GetMinorSign(ARBITRARY_ROW,col) * matrix.at(ARBITRARY_ROW, col) * _ComputeDeterminant(matrix.SubMartix(ARBITRARY_ROW, col));

				return determinant;
			}
		}

	private:

		union
		{
			ElementType mMatrix[ROWS][COLS];
			ElementType mElements[ROWS * COLS];
		};
	};

	template <typename T, size_t rows,size_t cols, MatrixVectors vectors>
	const MatrixBase<T, rows, cols, vectors> MatrixBase<T, rows, cols, vectors>::Zero = MatrixBase<T, rows, cols, vectors>::CreateZeroMatrix();
	
	template <typename T, size_t rows, size_t cols, MatrixVectors vectors>
	const MatrixBase<T, rows, cols, vectors> MatrixBase<T, rows, cols, vectors>::Identity = MatrixBase<T, rows, cols, vectors>::CreateIdentityMatrix();
}

#endif