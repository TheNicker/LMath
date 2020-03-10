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

#ifndef __LMATH_CONFIG_H__
#define __LMATH_CONFIG_H__ 1

//Do not use Cartesian component notation
#define LMATH_VECTOR_CARTESIAN_COMPONENT_DISABLE 0
//Use Cartesian component notation as getter e.g. vector.X()
#define LMATH_VECTOR_CARTESIAN_COMPONENT_GETTER 1
//Use Cartesian component notation as lower case plain variable  e.g. vector.x
#define LMATH_VECTOR_CARTESIAN_COMPONENT_PLAIN_VARIABLE_LOW_CASE 2
//Use Cartesian component notation as big case plain variable  e.g. vector.X
#define LMATH_VECTOR_CARTESIAN_COMPONENT_PLAIN_VARIABLE_BIG_CASE 4



//Default configuration

//Ogre conversion from Vector[N] and to Vector[N]
#ifndef LMATH_VECTOR_OGRE_EXTENSIONS
	#define LMATH_VECTOR_OGRE_EXTENSIONS 0
#endif
//Windows conversion from and to POINT and SIZE
#ifndef LMATH_VECTOR_WINDOWS_EXTENSIONS
	#define LMATH_VECTOR_WINDOWS_EXTENSIONS 0
#endif

#ifndef LMATH_VECTOR_CARTESIAN_COMPONENT
	#define LMATH_VECTOR_CARTESIAN_COMPONENT LMATH_VECTOR_CARTESIAN_COMPONENT_GETTER
#endif

#ifndef LMATH_ALWAYS_INITIALIZE_VARIABLES
	#define LMATH_ALWAYS_INITIALIZE_VARIABLES 0
#endif

#ifndef LMATH_ENABLE_MATRIX4_MUL_IN_VECTOR3
	#define LMATH_ENABLE_MATRIX4_MUL_IN_VECTOR3 1
#endif

#endif