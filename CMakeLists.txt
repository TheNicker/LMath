# CMakeList.txt : CMake project for Test, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project("lmath")
# Add source to this project's executable.
set(CMAKE_CXX_STANDARD 17)
include_directories("./Include")
include_directories("./Test")
add_subdirectory("Example")
add_subdirectory("Test")



# TODO: Add tests and install targets if needed.
