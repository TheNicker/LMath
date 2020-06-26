
#include <LMath/Vector.h>
#include <LMath/Bounds.h>
#include <LMath/Matrix.h>
#include <iostream>
#include <array>



int main()
{

	using namespace LMath;
	//----------------------------------------------------
	//Declare a set of helper class aliases.
	using Vec3D = VectorBase<double, 3>;
	using Vec2D = VectorBase<double, 2>;
	using Quat = QuaternionBase<double>;
	using Vec4D = VectorBase<double, 4>;
	//----------------------------------------------------


	// Construcut a vector with any argument permutation with smaller vector or numbers.
	Vec4D vec(Vec2D(1, 2), 3, 1);

	/// Rotate the vector 30 degrees yaw and 25 degrees pitch
	vec = Quat::FromEuler(Vec3D(25, 30, 0) * Constans::DegToRad) * vec;

	//----------------------------------------------------

	//vector supports cartesian notation up to four dimensions (x,y,z,w)
	Vec2D v2(1, 2);
	[[maybe_unused]] bool equalsOne = v2.X() == 1; // O.K
	//v2.Z() == 1; // compiler error, v2 has only X and Y components.


	//----------------------------------------------------
	//Any vector to any vector via explicit casting.
	using Vec5I = VectorBase<int, 5>; //5 dimenstion integer vector.

	//Conversion from 3D double vector to 5D integer vector integer to 3
	Vec5I v5i = static_cast<Vec5I>(Vec3D(2.3, 1.1, 7.8));
	std::cout << v5i.ToString();

	//----------------------------------------------------
	// Geometry operations
	Vec3D vec3d(1, 2, 3);
	Vec3D reflectionVector = vec3d.Reflect(Vec3D(0, 1, 0));
	[[maybe_unused]] Vec3D cross = vec3d.Cross(reflectionVector);

	//----------------------------------------------------

	//BoundsBase is based upon VectorBase, same template rules applies.

	//Declare a 2D rectangle type using float.
	using RectF = BoundsBase<float, 2>;

	RectF rect = { {3.f,2.f} , {5.f,6.f} };
	

	//Declare a 3D bounding box type using double.
	using BoundingBoxD = BoundsBase<double, 3>;

	//declate a bounding box instance using two vectors (min and max).
	BoundingBoxD bb = { {-1,-1,-1 },{5,5,5} };

	//test for overlapping, should return true.
	[[maybe_unused]]  bool contains = bb.RelationTo({ { 0, 0, 0 }, { 1, 1, 1 } }) == BoundingBoxD::Relation::Overlap;

}
