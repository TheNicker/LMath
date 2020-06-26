

#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#include <LMath/Bounds.h>


TEST_CASE("BOUNDS", "[Bounds]") 
{
	// For each section, vector v is anew:
	using namespace LMath;
	using BoundingBox = BoundsBase<double, 3>;
	BoundingBox bd = { { -2.0,-2.0,-2.0 }, { 5, 5, 6 } };
	BoundingBox bd1;
	


	SECTION("ensure boundsBase functionality") {
		REQUIRE(bd.Dimensions == 3);
		REQUIRE(bd.Inflate(-2) == BoundingBox ( { 0.0,0.0,0.0 }, { 3, 3, 4 } ));

		REQUIRE(BoundingBox({ 1.1,2.6,-3.2 }, { 9.7, 15.44, 3.1 }).Round() == BoundingBox({ 1.0,3.0,-3.0 }, { 10.0, 15.0, 3.0 }));
		
	}

	SECTION("test span") 
	{
		REQUIRE(bd.GetSpan(0) == 7);
		REQUIRE(bd.GetSpan(2) == 8);

	}

	SECTION("test intersection") 
	{
		REQUIRE(bd.IsInside({1,3,2} ) == true);
		REQUIRE(bd.IsInside({ 1,3,7 }) == false);
		REQUIRE(bd.IsInside({ 1,3,6 }) == false);
		REQUIRE(bd.IsInsideOrBoundry({ 1,3,6 }) == true);
		REQUIRE((bd.Intersection(BoundingBox({ 0,0,0 }, { 9, 15, 3 })) == BoundingBox({ 0, 0, 0 }, { 5,5,3 })) == true);
		REQUIRE(bd.RelationTo(BoundingBox({ -4.0, -5.0, -6.0 }, { -2.0,5.0,5.0 })) == BoundingBox::Relation::Adjacent);
		
		REQUIRE(bd.IntersectsIncludeBoundry(BoundingBox({ -4.0, -5.0, -6.0 }, { -2.0,5.0,5.0 })) == true);

		REQUIRE(bd.RelationTo(BoundingBox::VectorType( 1,3,6 )) == BoundingBox::Relation::Adjacent);
		REQUIRE(bd.RelationTo(BoundingBox::VectorType(1, 3, 7)) == BoundingBox::Relation::Disjoint);
		REQUIRE(bd.RelationTo(BoundingBox::VectorType(1, 3, 5)) == BoundingBox::Relation::Overlap);
		
		REQUIRE(bd.RelationTo(BoundingBox({ -4.0, -5.0, -6.0 }, { -2.0,5.0,5.0 })) == BoundingBox::Relation::Adjacent);
		

		REQUIRE(BoundingBox({ -2.0, -1, 0.0 }, { 5,5,5 }).Inflate(3).GetSpan(0) == 13.0);
		REQUIRE(BoundingBox({ -2.0, 7, 6.0 }, { 5.0,9.0,203.0 })
          .Join(BoundingBox({ -7.0, 0.0, 5.9 }, { 1.0,12.5,16.0 })) == BoundingBox({ -7.0, 0.0, 5.9 }, { 5,12.5,203.0 }));
	}

#if LMATH_VECTOR_WINDOWS_EXTENSIONS == 1


	SECTION("windows extensions")
	{
		POINT p{ 12,13 };
		using Vec2I = VectorBase<long, 2>;
		REQUIRE( static_cast<POINT>(Vec2I(12,13)) == static_cast<Vec2I>(p));

		REQUIRE( Vec2I(12, 13) == static_cast<Vec2I>(p));
	}
#endif

}


TEST_CASE("VECTOR", "[Vector]")
{
	using namespace LMath;
	using Vec3D = VectorBase<double, 3>;
	using Vec4D = VectorBase<double, 4>;
	SECTION("vector construction")
	{


		REQUIRE(Vec3D::Unit == Vec3D(1.0, 1.0, 1.0));

		REQUIRE(Vec3D::Zero == Vec3D(0.0, 0.0, 0.0));

		Vec4D::ElementType rawData[] = { 4,5,6,7 };

		REQUIRE(Vec4D(rawData) == Vec4D(4.0, 5.0, 6.0, 7.0));

		REQUIRE(Vec3D::Unit + Vec3D::Zero == Vec3D(1, 1, 1));

	}

	SECTION("vector scalar arithmetic")
	{
		Vec3D v1 = { 6.0,6.0,6.0 };
		Vec3D::ElementType scalar = 3.0;
		REQUIRE(v1 * scalar == Vec3D(18.0, 18.0, 18.0));
		REQUIRE(scalar * v1 == Vec3D(18.0, 18.0, 18.0));

		REQUIRE(v1 + scalar == Vec3D(9.0, 9.0, 9.0));
		REQUIRE(scalar + v1 == Vec3D(9.0, 9.0, 9.0));

		REQUIRE(v1 / scalar == Vec3D(2.0, 2.0, 2.0));
		REQUIRE(scalar / v1 == Vec3D(0.5, 0.5, 0.5));

		REQUIRE(v1 - scalar == Vec3D(3.0, 3.0, 3.0));
		REQUIRE(scalar - v1 == Vec3D(-3.0, -3.0, -3.0));

		REQUIRE(1.0 / Vec3D(3.0, 3.0, 3.0) == Vec3D(1.0 / 3.0));


		REQUIRE(-v1 == Vec3D(-6.0, -6.0, -6.0));


	}

	SECTION("vector arithmetic")
	{
		Vec3D v1 = { 6.0,6.0,6.0 };
		Vec3D v2 = { 2.0,2.0,2.0 };

		REQUIRE(v1 * v2 == Vec3D(12.0, 12.0, 12.0));
		REQUIRE(v1 / v2 == Vec3D(3.0, 3.0, 3.0));
		REQUIRE(v1 + v2 == Vec3D(8.0, 8.0, 8.0));
		REQUIRE(v1 - v2 == Vec3D(4.0, 4.0, 4.0));
	}


	SECTION("vector operations")
	{
		Vec3D v1 = { 3.0,4.0,0.0 };
		Vec3D v2 = { 1.0,3.0,1.0 };
		REQUIRE(v1.SquaredLength() == 25.0);
		REQUIRE(v1.Length() == 5.0);
		REQUIRE(v1.Dot(v2) == 15.0);

		REQUIRE(v1.Max(v2) == Vec3D(3.0, 4.0, 1.0));
		REQUIRE(v2.Max(v1) == Vec3D(3.0, 4.0, 1.0));

		REQUIRE(v1.Min(v2) == Vec3D(1.0, 3.0, 0.0));
		REQUIRE(v2.Min(v1) == Vec3D(1.0, 3.0, 0.0));

	}

	SECTION("vector euclid getters")
	{
		Vec4D v1 = { 3.0,4.0, 8.0, 6.0 };
		REQUIRE(v1.X() == 3.0);
		REQUIRE(v1.Y() == 4.0);
		REQUIRE(v1.Z() == 8.0);
		REQUIRE(v1.W() == 6.0);

		v1.X() = 1.3;
		v1.Y() = 2.6;
		v1.Z() = 7.4;
		v1.W() = 33.5;

		REQUIRE(v1.X() == 1.3);
		REQUIRE(v1.Y() == 2.6);
		REQUIRE(v1.Z() == 7.4);
		REQUIRE(v1.W() == 33.5);

		REQUIRE(v1 == Vec4D(1.3, 2.6, 7.4, 33.5));
	}

	SECTION("vector geometric")
	{
		VectorBase<double, 4> vec(2, -3, 4, 1);
		VectorBase<double, 4> normalizedReference = Vec4D(std::sqrt(2.0 / 15.0), -std::sqrt(3.0 / 10.0), 2.0 * std::sqrt(2.0 / 15.0), std::sqrt(1.0 / 30.0));
		VectorBase<double, 4> normalized = vec.Normalized();
		REQUIRE((normalized == normalizedReference));

	}

	SECTION("vector reflection")
	{
		VectorBase<double, 3> vec(0.5, -0.5, 0);
		VectorBase<double, 3> normal(0, 1, 0);
		VectorBase<double, 3> reflected = vec.Reflect(normal);
		REQUIRE((reflected == VectorBase<double, 3>(0.5,0.5,0.0)));
	}

}
	TEST_CASE("Math function", "[math function]")
	{
		using namespace LMath;
		SECTION("Math")
		{
			REQUIRE(NextPowerOfTwo(0u) == 0);
			REQUIRE(NextPowerOfTwo(1u) == 1);
			REQUIRE(NextPowerOfTwo(2u) == 2);
			REQUIRE(NextPowerOfTwo(3u) == 4);
			REQUIRE(NextPowerOfTwo(4u) == 4);
			REQUIRE(NextPowerOfTwo(5u) == 8);
			REQUIRE(NextPowerOfTwo(8u) == 8);
			REQUIRE(NextPowerOfTwo(9u) == 16);
		}
	}
