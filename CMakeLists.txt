cmake_minimum_required(VERSION 3.5)
project(BallTree)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp BallTree.cpp BallTree.h Utility.cpp Utility.h)
add_executable(BallTree ${SOURCE_FILES})