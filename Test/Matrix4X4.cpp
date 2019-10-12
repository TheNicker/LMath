#include <catch2/catch.hpp>
#include <LMath/Matrix.h>

using Matrix4X4 = LMath::MatrixBase<double, 4, 4>;

#define CHECK_MATRIX(a, b) \
	 for (size_t i = 0 ; i < decltype(a)::Rows * decltype(a)::Cols ; i++) \
	{\
    CHECK(a.at(i) == Approx(b.at(i))); \
	 };

TEST_CASE("Matrix4x4 Unary  - operator", "[Matrix4x4]")
{
	Matrix4X4 m1(1, 0, 0, 1,
		0, -2, 1, 2,
		2, 1, 0, 1,
		-2, 0, 1, -4);

	Matrix4X4 reference(-1, 0, 0, -1,
		0, 2, -1, -2,
		-2, -1, 0, -1,
		2, 0, -1, 4);

	Matrix4X4 m2 = -m1;

	CHECK_MATRIX(m2, reference);

}

TEST_CASE("Matrix4x4 inverse", "[Matrix4x4]")
{
	Matrix4X4 m1(1, 0, 0, 1,
		 		0, 2, 1, 2, 
				2, 1, 0, 1, 
				2, 0, 1, 4);
	
	Matrix4X4 reference(-2, -0.5, 1, 0.5,
				1 , 0.5 ,0.0, -0.5 
				, -8 , -1 , 2 , 2
				, 3 , 0.5 , -1, -0.5);
	 m1 = m1.Inverse();

	 CHECK_MATRIX(m1, reference);

}


TEST_CASE("Matrix4x4 multiply", "[Matrix4x4]")
{

	
	Matrix4X4 m1(7, 15, 3, 11
				, 5, 1, 7, -8
				,  -1, 10, 2, 12
				,  13, 14, 15, 8);

	Matrix4X4 m2(
		7, -1, 5, 9
		, 2, -5, 0, 1
		, -8, 6, 0, -8
		, 1, 8, 3, 9
	);

	Matrix4X4 reference(
		66, 24, 68, 153
		, -27, -32, 1, -82
		, 9, 59, 31, 93
		, 7, 71, 89, 83);

	Matrix4X4 m4 = m1 * m2;

	CHECK_MATRIX(m4, reference);

}



TEST_CASE("Matrix4x4 explicit cast", "[Matrix4x4]")
{
	Matrix4X4 m1(
	66, 24, 68, 153
	, -27, -32, 1, -82
	, 9, 59, 31, 93
	, 7, 71, 89, 83);


	using Matrix3X3 = LMath::MatrixBase<double, 3, 3>;


	Matrix3X3 reducedReference (
		66, 24, 68, 
		 -27, -32, 1, 
		 9, 59, 31);

	Matrix3X3 reduced = static_cast<Matrix3X3>(m1);

	

	CHECK_MATRIX(reduced, reducedReference);

	using Matrix6X6 = LMath::MatrixBase<double, 6, 6>;

	Matrix6X6 exapndedReference(
		66, 24, 68,0,0,0
	   ,-27, -32, 1, 0, 0,0
		,9, 59, 31, 0, 0, 0
		, 0, 0, 0, 0, 0, 0
		, 0, 0, 0, 0, 0, 0
		, 0, 0, 0, 0, 0, 0
	);


	Matrix6X6 expanded = static_cast<Matrix6X6>(reduced);
	CHECK_MATRIX(expanded, exapndedReference);
};
