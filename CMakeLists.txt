cmake_minimum_required(VERSION 3.7)
project(numerical_analysis)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

include_directories(include)

set(SOURCE_FILES test/main.cpp include/FunctionUtils.hpp)
add_executable(numerical_analysis ${SOURCE_FILES})
