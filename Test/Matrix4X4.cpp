#include <catch2/catch.hpp>
#include <LMath/Matrix.h>

using Matrix4X4 = LMath::MatrixBase<double, 4, 4>;

#define CHECK_MATRIX_4(a, b) \
    CHECK(a.at(0,0) == Approx(b.at(0,0))); \
    CHECK(a.at(0,1) == Approx(b.at(0,1))); \
    CHECK(a.at(0,2) == Approx(b.at(0,2))); \
	CHECK(a.at(0,3) == Approx(b.at(0,3))); \
    CHECK(a.at(1,0) == Approx(b.at(1,0))); \
    CHECK(a.at(1,1) == Approx(b.at(1,1))); \
    CHECK(a.at(1,2) == Approx(b.at(1,2))); \
	CHECK(a.at(1,3) == Approx(b.at(1,3))); \
    CHECK(a.at(2,0) == Approx(b.at(2,0))); \
    CHECK(a.at(2,1) == Approx(b.at(2,1))); \
    CHECK(a.at(2,2) == Approx(b.at(2,2))); \
	CHECK(a.at(2,3) == Approx(b.at(2,3))); \
	CHECK(a.at(3,0) == Approx(b.at(3,0))); \
	CHECK(a.at(3,1) == Approx(b.at(3,1))); \
	CHECK(a.at(3,2) == Approx(b.at(3,2))); \
	CHECK(a.at(3,3) == Approx(b.at(3,3))); 

TEST_CASE("Matrix4x4 inverse", "[Matrix3x3]")
{
	Matrix4X4 m1(1, 0, 0, 1,
		 		0, 2, 1, 2, 
				2, 1, 0, 1, 
				2, 0, 1, 4);
	
	Matrix4X4 m2(-2, -0.5, 1, 0.5,
				1 , 0.5 ,0.0, -0.5 
				, -8 , -1 , 2 , 2
				, 3 , 0.5 , -1, -0.5);
	 m1 = m1.Inverse();

	 CHECK_MATRIX_4(m1, m2);

}
