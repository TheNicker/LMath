# CMakeList.txt : CMake project for Test, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project("LMathTestCase")

file(GLOB_RECURSE sourceFiles
    "./*.cpp"
)

# Add source to this project's executable.
add_executable (Test ${sourceFiles} )

# TODO: Add tests and install targets if needed.
