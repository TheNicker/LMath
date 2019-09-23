![LMath logo](https://github.com/TheNicker/blob/blob/master/Lmath.png)

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a7f043160b7b4f94b9ec87f83c1608b2)](https://www.codacy.com/manual/TheNicker/LMath?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=TheNicker/LMath&amp;utm_campaign=Badge_Grade)
[![Build status](https://ci.appveyor.com/api/projects/status/yxypaqf9is7w1iwy?svg=true)](https://ci.appveyor.com/project/LiorL/lmath)
![GitHub](https://img.shields.io/github/license/Thenicker/lmath)

**LMath** is a modern and lightweight cross-platform header only **C++ 17** math and geometry library focusing on performance, safety and ease of use.

## LMath is aiming to

* Satisfy any use case with minimum code using templates.
* Have many opted-in features that can be removed at precompile.
* Gradually improve efficiency using its backend for better algorithms and targeting specific platforms.
* Being backward compatible.

## Highlights

* Built from the grounds up using modern techniques for safety and convenience for an intuitive API.
* Features are optional, reducing clutter and compilation time while solving large number of use cases.
* High code reusability minimizing maintenance.
* Compatible with existing systems for easy integration, e.g. use Cartesian notation or any user defined accessor notation.
* Designed for 3D graphics.
* Interface/ implementation separation, allows for instantiation specialization e.g. implementing SIMD instructions for matrix 4X4 multiplication for **any** data type while not changing the Interface.

**Note**: compiles in GCC Clang and MSVC toolset version 14.22 and prior. version 14.23 and onwards is not stable.

## Library overview

|      Class       |                                           Description                                                         |
| ---------------  | ------------------------------------------------------------------------------------------------------------- |
| **Element**      | Base class for all other classes used only for data manipulation                                              |
| **Vector**       | Used for vector arithmetic logic                                                                              |
| **Quaternion**   | A specialization of VectorBase<T,4> with several extra methods                                                |
| **Matrix**       | Provides Matrix N X M of an arbitrary type                                                                    |
| **Bounds**       | Arbitrary dimension bounding limits  e.g. Line for 1D, rectangle for 2D or a box for 3D                       |
| **LMathConfig**  | Configuration file for changing the behvaiour of the library e.g. removing cartesian notation syntactic sugar |


## Vector and Quaternion examples

```c++

//Import LMath namespace
using namespace LMath;
//----------------------------------------------------
//Declare a set of helper class aliases.
using Vec3D = VectorBase<double, 3>;
using Vec2D = VectorBase<double, 2>;
using Quat = QuaternionBase<double>;
using Vec4D = VectorBase<double, 4>;
//----------------------------------------------------


// Construcut a vector with any argument permutation with smaller vectors or numbers.
Vec4D vec(Vec2D(1, 2), 3, 1);

/// Rotate the vector 30 degrees yaw and 25 degrees pitch
vec = Quat::FromEuler(Vec3D(25, 30, 0) * Constans::DegToRad) * vec;

//----------------------------------------------------

//vector supports cartesian notation up to four dimensions (x,y,z,w)
Vec2D v2(1, 2);
v2.X() == 1; // O.K
//v2.Z() == 1; // compiler error, v2 has only X and Y components.


//----------------------------------------------------
//Any vector to any vector implicit casting.
using Vec5I = VectorBase<int, 5>; //5 dimenstion integer vector.

//Conversion from 3D double vector to 5D integer vector integer to 3
Vec5I v5i = static_cast<Vec5I>(Vec3D(2.3, 1.1, 7.8));
std::cout << v5i.ToString();

//----------------------------------------------------
// Geometry operations
Vec3D vec3d(1, 2, 3);
Vec3D reflectionVector = vec3d.Reflect(Vec3D(0, 1, 0));
Vec3D cross = vec3d.Cross(reflectionVector);
```

## Bounds examples

```
//BoundsBase is based upon VectorBase, same template rules applies.

//Declare a 2D rectangle type using float.
using RectF = BoundsBase<float, 2>;

//Declare a 3D bounding box type using double.
using BoundingBoxD = BoundsBase<double, 3>;

//declate a bounding box instance using two vectors (min and max).
BoundingBoxD bb = { {-1,-1,-1 },{5,5,5} };

//test for overlapping, should return true.
bool contains = bb.Contains({ { 0, 0, 0 }, { 1, 1, 1 } }) == BoundingBoxD::IntersectState::Overlap;
```

## Getting involved

<a href="https://discord.gg/6QNHmQR"> <img src="https://discordapp.com/assets/f8389ca1a741a115313bede9ac02e2c0.svg" width="48" height="48" /> Visit our discord development forum</a>

<a href="https://trello.com/b/XQQ2CI8t/lmath"><img src="https://cdn1.iconfinder.com/data/icons/designer-skills/128/trello-64.png" width="48" height="48" /> View or join our trello development board</a>

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
