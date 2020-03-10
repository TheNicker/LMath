/**
 *  ============================================================================
 *  MIT License
 *
 *  Copyright (c) 2016 Eric Phillips
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a
 *  copy of this software and associated documentation files (the "Software"),
 *  to deal in the Software without restriction, including without limitation
 *  the rights to use, copy, modify, merge, publish, distribute, sublicense,
 *  and/or sell copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 *  DEALINGS IN THE SOFTWARE.
 *  ============================================================================
 *
 *
 *  This file contains test cases for the Matrix3x3 functions.
 *
 *  Created by Eric Phillips on November 8, 2016.
 */

#include "catch2/catch.hpp"
#include <LMath/Matrix.h>




using Matrix3x3 = LMath::MatrixBase<double, 3, 3>;
using Quaternion = LMath::QuaternionBase<double>;

#define CHECK_MATRIX(a, b) \
	 for (size_t i = 0 ; i < decltype(a)::Rows * decltype(a)::Cols ; i++) \
	{\
    CHECK(a.at(i) == Approx(b.at(i))); \
	 };


#define CHECK_SCALE(a,s0,s1,s2) CHECK( (a.at(0,0) == s0 && a.at(1,1) == s1 && a.at(2,2) == s2))

TEST_CASE("Matrix3x3 plus scalar", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    m1 = m1 + 3;
    Matrix3x3 m2 = Matrix3x3(5, -2, 6, 10, 4, -3, -6, 7, 11);
    CHECK_MATRIX(m1, m2);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    m1 = m1 + 0.2;
    m2 = Matrix3x3(0.44, 0.2082, -0.1, 0.8, -1.2, 0.93, 0.58, 0.296, 0.04);
    CHECK_MATRIX(m1, m2);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    m1 = m1 + 13;
    m2 = Matrix3x3(-14, 96, 45, -140, 66, 96, -51, 36, -33);
    CHECK_MATRIX(m1, m2);
}

TEST_CASE("Matrix3x3 minus scalar", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    m1 = m1 - 3;
    Matrix3x3 m2 = Matrix3x3(-1, -8, 0, 4, -2, -9, -12, 1, 5);
    CHECK_MATRIX(m1, m2);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    m1 = m1 - 0.2;
    m2 = Matrix3x3(0.04, -0.1918, -0.5, 0.4, -1.6, 0.53, 0.18, -0.104, -0.36);
    CHECK_MATRIX(m1, m2);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    m1 = m1 - 13;
    m2 = Matrix3x3(-40, 70, 19, -166, 40, 70, -77, 10, -59);
    CHECK_MATRIX(m1, m2);
}

TEST_CASE("Matrix3x3 times scalar", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    
    
    m1 = m1 * 3.0;
    Matrix3x3 m2 = Matrix3x3(6, -15, 9, 21, 3, -18, -27, 12, 24);
    CHECK_MATRIX(m1, m2);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    m1 = m1 * 0.2;
    m2 = Matrix3x3(0.048, 0.00164, -0.06, 0.12, -0.28, 0.146, 0.076, 0.0192,
        -0.032);
    CHECK_MATRIX(m1, m2);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    m1 = m1 * 13.0;
    m2 = Matrix3x3(-351, 1079, 416, -1989, 689, 1079, -832, 299, -598);
    CHECK_MATRIX(m1, m2);
}

TEST_CASE("Matrix3x3 divided by scalar", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    m1 = m1 / 3;
    Matrix3x3 m2 = Matrix3x3(0.6666666667, -1.6666666667, 1, 2.3333333333,
        0.3333333333, -2, -3, 1.3333333333, 2.6666666667);
    CHECK_MATRIX(m1, m2);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    m1 = m1 / 0.2;
    m2 = Matrix3x3(1.2, 0.041, -1.5, 3, -7, 3.65, 1.9, 0.48, -0.8);
    CHECK_MATRIX(m1, m2);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    m1 = m1 / 13;
    m2 = Matrix3x3(-2.0769230769, 6.3846153846, 2.4615384615, -11.7692307692,
        4.0769230769, 6.3846153846, -4.9230769231, 1.7692307692, -3.5384615385);
    CHECK_MATRIX(m1, m2);
}

TEST_CASE("Matrix3x3 plus Matrix3x3", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    Matrix3x3 m2 = Matrix3x3(-14, 96, 45, -140, 66, 96, -51, 36, -33);
    Matrix3x3 m3 = m1 + m2;
    Matrix3x3 m4 = Matrix3x3(-12, 91, 48, -133, 67, 90, -60, 40, -25);
    CHECK_MATRIX(m3, m4);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    m2 = Matrix3x3(1.7, 0.003, -5.6, 0.44, 0.2082, -0.1, 0.8, -1.2, 0.93);
    m3 = m1 + m2;
    m4 = Matrix3x3(1.94, 0.0112, -5.9, 1.04, -1.1918, 0.63, 1.18, -1.104, 0.77);
    CHECK_MATRIX(m3, m4);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    m2 = Matrix3x3(5, -2, 6, 10, 4, -3, -6, 7, 11);
    m3 = m1 + m2;
    m4 = Matrix3x3(-22, 81, 38, -143, 57, 80, -70, 30, -35);
    CHECK_MATRIX(m3, m4);
}

TEST_CASE("Matrix3x3 minus Matrix3x3", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    Matrix3x3 m2 = Matrix3x3(-14, 96, 45, -140, 66, 96, -51, 36, -33);
    Matrix3x3 m3 = m1 - m2;
    Matrix3x3 m4 = Matrix3x3(16, -101, -42, 147, -65, -102, 42, -32, 41);
    CHECK_MATRIX(m3, m4);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    m2 = Matrix3x3(1.7, 0.003, -5.6, 0.44, 0.2082, -0.1, 0.8, -1.2, 0.93);
    m3 = m1 - m2;
    m4 = Matrix3x3(-1.46, 0.0052, 5.3, 0.16, -1.6082, 0.83, -0.42, 1.296,
        -1.09);
    CHECK_MATRIX(m3, m4);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    m2 = Matrix3x3(5, -2, 6, 10, 4, -3, -6, 7, 11);
    m3 = m1 - m2;
    m4 = Matrix3x3(-32, 85, 26, -163, 49, 86, -58, 16, -57);
    CHECK_MATRIX(m3, m4);
}

TEST_CASE("Matrix3x3 times Matrix3x3", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    Matrix3x3 m2 = Matrix3x3(-14, 96, 45, -140, 66, 96, -51, 36, -33);
    Matrix3x3 m3 = m1 * m2;
    Matrix3x3 m4 = Matrix3x3(519, -30, -489, 68, 522, 609, -842, -312, -285);
    CHECK_MATRIX(m3, m4);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    m2 = Matrix3x3(1.7, 0.003, -5.6, 0.44, 0.2082, -0.1, 0.8, -1.2, 0.93);
    m3 = m1 * m2;
    m4 = Matrix3x3(0.1716080, 0.3624272, -1.6238200, 0.9880000, -1.1656800,
        -2.5411000, 0.5602400, 0.2131272, -2.2864000);
    CHECK_MATRIX(m3, m4);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    m2 = Matrix3x3(5, -2, 6, 10, 4, -3, -6, 7, 11);
    m3 = m1 * m2;
    m4 = Matrix3x3(503, 610, -59, -733, 1099, -164, 186, -102, -959);
    CHECK_MATRIX(m3, m4);
}

TEST_CASE("Matrix3x3 times Vector3", "[Matrix3x3]")
{
	using Vector3 = Matrix3x3::VectorType;
    // Case 1
    Matrix3x3 m = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    Vector3 v = Vector3(0.3535534, -0.1464466, 0.3535534);
    Vector3 r = m * v;
    CHECK(r.X() == Approx(2.5));
    CHECK(r.Y() == Approx(0.2071068));
    CHECK(r.Z() == Approx(-0.9393398));
    // Case 2
    m = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    v = Vector3(0.53, -0.0532, 0.22);
    r = m * v;
    CHECK(r.X() == Approx(0.0607638));
    CHECK(r.Y() == Approx(0.5530800));
    CHECK(r.Z() == Approx(0.1610928));
    // Case 3
    m = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    v = Vector3(13, -1.23, 3.4);
    r = m * v;
    CHECK(r.X() == Approx(-344.29));
    CHECK(r.Y() == Approx(-1771.99));
    CHECK(r.Z() == Approx(-1016.69));
}
//
TEST_CASE("Matrix3x3 equality", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    Matrix3x3 m2 = Matrix3x3(-14, 96, 45, -140, 66, 96, -51, 36, -33);
    CHECK_FALSE(m1 == m2);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    m2 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    CHECK(m1 == m2);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    m2 = Matrix3x3(-27, 83, 33, -153, 53, 83, -64, 23, -46);
    CHECK_FALSE(m1 == m2);
}

TEST_CASE("Matrix3x3 inequality", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    Matrix3x3 m2 = Matrix3x3(-14, 96, 45, -140, 66, 96, -51, 36, -33);
    CHECK(m1 != m2);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    m2 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    CHECK_FALSE(m1 != m2);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    m2 = Matrix3x3(-27, 83, 33, -153, 53, 83, -64, 23, -46);
    CHECK(m1 != m2);
}

TEST_CASE("Matrix3x3 transpose", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
	m1 = m1.Transpose(); 
    Matrix3x3 m2 = Matrix3x3(2, 7, -9, -5, 1, 4, 3, -6, 8);
    CHECK_MATRIX(m1, m2);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
	m1 = m1.Transpose();
    m2 = Matrix3x3(0.24, 0.6, 0.38, 0.0082, -1.4, 0.096, -0.3, 0.73, -0.16);
    CHECK_MATRIX(m1, m2);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
	m1 = m1.Transpose();
    m2 = Matrix3x3(-27, -153, -64, 83, 53, 23, 32, 83, -46);
    CHECK_MATRIX(m1, m2);
}


TEST_CASE("Matrix3x3 scale", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    Matrix3x3 m2 = Matrix3x3(-14, 96, 45, -140, 66, 96, -51, 36, -33);
	Matrix3x3 m3 = m1.Scale(m2);
    Matrix3x3 m4 = Matrix3x3(-28, -480, 135, -980, 66, -576, 459, 144, -264);
    CHECK_MATRIX(m3, m4);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    m2 = Matrix3x3(1.7, 0.003, -5.6, 0.44, 0.2082, -0.1, 0.8, -1.2, 0.93);
    m3 = Matrix3x3(m1.Scale(m2));
    m4 = Matrix3x3(0.408, 0.0000246, 1.68, 0.264, -0.29148, -0.073, 0.304,
        -0.1152, -0.1488);
    CHECK_MATRIX(m3, m4);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    m2 = Matrix3x3(5, -2, 6, 10, 4, -3, -6, 7, 11);
    m3 = m1.Scale(m2);
    m4 = Matrix3x3(-135, -166, 192, -1530, 212, -249, 384, 161, -506);
    CHECK_MATRIX(m3, m4);


	//case 4 Scale matrix
	m1.SetScale(4, 5, 6);
	CHECK_SCALE(m1, 4, 5, 6);
	m1.SetScale({ 7,8,9 });
	CHECK_SCALE(m1, 7,8,9);
	m1 = Matrix3x3::CreateScaleMatrix(10, 11, 12);
	CHECK_SCALE(m1, 10, 11 ,12 );
	m1 = Matrix3x3::CreateScaleMatrix({ 13, 14, 15 });
	CHECK_SCALE(m1, 13, 14, 15);
}



TEST_CASE("Matrix3x3 determinate", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
	double v = m.Determinant();
    CHECK(v == Approx(185));
    // Case 2
    m = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
	v = m.Determinant();
    CHECK(v == Approx(-0.1368773));
    // Case 3
    m = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
	v = m.Determinant();
    CHECK(v == Approx(-911745));
}

TEST_CASE("Matrix3x3 inverse", "[Matrix3x3]")
{

    // Case 1
    Matrix3x3 m1 = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
	m1 = m1.Inverse();
    Matrix3x3 m2 = Matrix3x3(0.172973, 0.2810811, 0.1459459, -0.0108108,
        0.2324324, 0.1783784, 0.2, 0.2, 0.2);
    CHECK_MATRIX(m1, m2);
    // Case 2
    m1 = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
	m1 = m1.Inverse();
    m2 = Matrix3x3(-1.1245106, 0.2008222, 3.0247085, -2.7279903, -0.5523194,
        2.5950245, -4.3075069, 0.1455610, 2.4906975);
    CHECK_MATRIX(m1, m2);
    // Case 3
    m1 = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
	m1 = m1.Inverse();
    m2 = Matrix3x3(0.0047678, -0.0049948, -0.0056957, 0.0135455, -0.0036085,
        0.002912, 0.0001393, 0.0051451, -0.0123587);
    CHECK_MATRIX(m1, m2);
}

TEST_CASE("Matrix3x3 is invertible", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m = Matrix3x3(2, -5, 3, 7, 1, -6, -9, 4, 8);
    CHECK(m.IsInvertible());
    // Case 2
    m = Matrix3x3(0.24, 0.0082, -0.3, 0.6, -1.4, 0.73, 0.38, 0.096, -0.16);
    CHECK(m.IsInvertible());
    // Case 3
    m = Matrix3x3(-27, 83, 32, -153, 53, 83, -64, 23, -46);
    CHECK(m.IsInvertible());
    // Case 4
    m = Matrix3x3(2, 2, 3, 6, 6, 9, 1, 4, 8);
    CHECK_FALSE(m.IsInvertible());
    // Case 5
    m = Matrix3x3(1, 0, 0, -2, 0, 0, 4, 6, 1);
    CHECK_FALSE(m.IsInvertible());
}


TEST_CASE("Matrix3x3 from Quaternion", "[Matrix3x3]")
{
    // Case 1
    Quaternion q (0.3535534, -0.1464466, 0.3535534, 0.8535535);
    Matrix3x3 m1 = Matrix3x3::FromQuaternion(q);
    Matrix3x3 m2 = Matrix3x3(0.707107, -0.707107, 0, 0.5, 0.5, -0.707107, 0.5,
        0.5, 0.707107);
    CHECK_MATRIX(m1, m2);
    // Case 2
    q = Quaternion(0.3919183, 0.3196269, -0.8430416, -0.1830837);
    m1 = Matrix3x3::FromQuaternion(q);
    m2 = Matrix3x3(-0.625761, -0.0581591, -0.777844, 0.55923, -0.728638,
        -0.39541, -0.54377, -0.682425, 0.488477);
    CHECK_MATRIX(m1, m2);
    // Case 3
    q = Quaternion(0, 0.7071068, 0, 0.7071068);
    m1 = Matrix3x3::FromQuaternion(q);
    m2 = Matrix3x3(0, 0, 1, 0, 1, 0, -1, 0, 0);
    CHECK_MATRIX(m1, m2);
    // Case 4
    q = Quaternion(-27, 83, 32, -153);
    m1 = Matrix3x3::FromQuaternion(q);
    m2 = Matrix3x3(0.506224, 0.165673, -0.846339, -0.445353, 0.890612,
        -0.0920408, 0.73851, 0.423513, 0.524633);
    CHECK_MATRIX(m1, m2);
}

TEST_CASE("Matrix3x3 to Quaternion", "[Matrix3x3]")
{
    // Case 1
    Matrix3x3 m = Matrix3x3(0.707107, -0.707107, 0, 0.5, 0.5, -0.707107, 0.5,
        0.5, 0.707107);
	Quaternion q = m.ToQuaternion();
    CHECK(q.X() == Approx(0.3535534));
    CHECK(q.Y() == Approx(-0.1464466));
    CHECK(q.Z() == Approx(0.3535534));
    CHECK(q.W() == Approx(0.8535535));
    // Case 2
    m = Matrix3x3(-0.625761, -0.0581591, -0.777844, 0.55923, -0.728638,
        -0.39541, -0.54377, -0.682425, 0.488477);
	q = m.ToQuaternion();
    CHECK(q.X() == Approx(-0.3919183));
    CHECK(q.Y() == Approx(-0.3196269));
    CHECK(q.Z() == Approx(0.8430416));
    CHECK(q.W() == Approx(0.1830837));
    // Case 3
    m = Matrix3x3(0, 0, 1, 0, 1, 0, -1, 0, 0);
	q = m.ToQuaternion();
    CHECK(q.X() == Approx(0));
    CHECK(q.Y() == Approx(0.7071068));
    CHECK(q.Z() == Approx(0));
    CHECK(q.W() == Approx(0.7071068));
    // Case 4
    m = Matrix3x3(0.506224, 0.165673, -0.846339, -0.445353, 0.890612,
        -0.0920408, 0.73851, 0.423513, 0.524633);
	q = m.ToQuaternion();
    CHECK(q.X() == Approx(0.150814));
    CHECK(q.Y() == Approx(-0.463615));
    CHECK(q.Z() == Approx(-0.178743));
    CHECK(q.W() == Approx(0.854615));
}
