cmake_minimum_required(VERSION 3.7)
project(numerical_analysis)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

# this OpenMP hack was found here https://stackoverflow.com/a/12404666
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(include)

set(SOURCE_FILES test/main.cpp include/FunctionUtils.hpp include/Optimizer.hpp include/VolumousObject.hpp)
add_executable(numerical_analysis ${SOURCE_FILES})