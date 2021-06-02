cmake_minimum_required(VERSION 3.13.4)
project(treenome VERSION 2.0.0)

find_package(OpenMP)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY out)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -g -std=c++11")

set(SOURCE_FILES 
  ${PROJECT_SOURCE_DIR}/src/cli.cpp
  ${PROJECT_SOURCE_DIR}/src/GTree.cpp
  ${PROJECT_SOURCE_DIR}/src/InputFile.cpp
  ${PROJECT_SOURCE_DIR}/src/Node.cpp
  ${PROJECT_SOURCE_DIR}/src/SeqRead.cpp
  ${PROJECT_SOURCE_DIR}/src/TreeTop.cpp)

add_executable(treenome ${PROJECT_SOURCE_DIR}/src/main.cpp ${SOURCE_FILES})
