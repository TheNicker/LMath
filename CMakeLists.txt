# CMakeList.txt : CMake project for Test, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project("lmath")
# Add source to this project's executable.
set(CMAKE_CXX_STANDARD 20)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    #enabled warnings flags
    add_compile_options(-Wextra -Wpedantic -Wall -Wpessimizing-move)
    #specifically disable the following warning flags
    add_compile_options(-Wno-c++98-compat -Wno-c++11-compat -Wno-c++14-compat -Wno-c++98-compat-pedantic -Wno-reserved-id-macro -Wno-newline-eof -Wno-global-constructors -Wno-exit-time-destructors -Wno-float-equal -Wno-covered-switch-default -Wno-unreachable-code-break)

  # using Clang
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # using GCC
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  # using Intel C++
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
# using Visual Studio C++
    add_compile_options(/permissive-) # Confrom to standards.
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
endif()


include_directories("./Include")
include_directories("./Test")
add_subdirectory("Example")
add_subdirectory("Test")



# TODO: Add tests and install targets if needed.
