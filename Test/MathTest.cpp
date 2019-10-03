
#include <catch2/catch.hpp>
#include <LMath/Math.h>

using namespace LMath;
TEST_CASE("Test mathematical modulo", "[Math]")
{
	REQUIRE(Modulo(14u, 4u) == 2);
	REQUIRE(Modulo(-14, 4) == 2);
	REQUIRE(Modulo(7u, 3u) == 1);
	REQUIRE(Modulo(-7, 3) == 2);
	REQUIRE(Modulo(20u, 7u) == 6);
	REQUIRE(Modulo(-20, 7) == 1);
	REQUIRE(Modulo(20.0, 7.0) == 6.0);
	REQUIRE(Modulo(-20.0, 7.0) == 1.0);
}

TEST_CASE("Next power of two ", "[Math]")
{
	REQUIRE(NextPowerOfTwo(0u) == 0);
	REQUIRE(NextPowerOfTwo(1u) == 1);
	REQUIRE(NextPowerOfTwo(2u) == 2);
	REQUIRE(NextPowerOfTwo(3u) == 4);
	REQUIRE(NextPowerOfTwo(4u) == 4);
	REQUIRE(NextPowerOfTwo(5u) == 8);
	REQUIRE(NextPowerOfTwo(8u) == 8);
	REQUIRE(NextPowerOfTwo(9u) == 16);
	REQUIRE(NextPowerOfTwo(111u) == 128);
}

TEST_CASE("Test sign", "[Math]")
{
	REQUIRE(Sign(0) == 0);
	REQUIRE(Sign(1) == 1);
	REQUIRE(Sign(3) == 1);
	REQUIRE(Sign(-1) == -1);
	REQUIRE(Sign(-3) == -1);
}

TEST_CASE("Is power of two", "[Math]")
{
	REQUIRE(IsPowerOfTwo(0u) == true);
	REQUIRE(IsPowerOfTwo(1u) == true);
	REQUIRE(IsPowerOfTwo(2u) == true);
	REQUIRE(IsPowerOfTwo(3u) == false);
	REQUIRE(IsPowerOfTwo(4u) == true);
	REQUIRE(IsPowerOfTwo(5u) == false);
	REQUIRE(IsPowerOfTwo(8u) == true);
	REQUIRE(IsPowerOfTwo(9u) == false);
	REQUIRE(IsPowerOfTwo(111u) == false);
	REQUIRE(IsPowerOfTwo(128u) == true);

}