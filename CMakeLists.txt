cmake_minimum_required(VERSION 3.11)
project(note_fonte)

find_package(OpenMP REQUIRED)
find_package(GTest REQUIRED)

set(CMAKE_BUILD_TYPE Release)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

file(GLOB src ${CMAKE_CURRENT_SOURCE_DIR}/src/*.C)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/examples")