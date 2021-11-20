cmake_minimum_required(VERSION 3.17)
project(Numerical VERSION 1.0 LANGUAGES CXX CUDA C)

set(CMAKE_CXX_STANDARD 17)

include(cmake/InitUCMake.cmake)

Ubpa_InitUCMake()
Ubpa_InitProject()

find_package(Vulkan REQUIRED FATAL_ERROR)
find_package(glfw3 REQUIRED FATAL_ERROR)


message(STATUS "find package: eigen")
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)
  message(STATUS "Building Eigen")
  set(_address "https://gitlab.com/libeigen/eigen.git")
  message(STATUS "fetch: ${_address}")
  FetchContent_Declare(
    eigen3
    GIT_REPOSITORY ${_address}
  )
  FetchContent_MakeAvailable(eigen3)
  message(STATUS "Eigen3 build done")
endif()

message(STATUS "find package: open_mesh")
find_package(OpenMesh QUIET)
if(NOT OpenMesh_FOUND)
  message(STATUS "Building OpenMesh")
  set(_address "https://www.graphics.rwth-aachen.de:9000/OpenMesh/OpenMesh.git")
  message(STATUS "fetch: ${_address}")
  FetchContent_Declare(
    OpenMesh
    GIT_REPOSITORY ${_address}
  )
  FetchContent_MakeAvailable(OpenMesh)
  message(STATUS "OpenMesh build done")
endif()

message(STATUS "find package: assimp")
find_package(Assimp QUIET)
if(NOT Assimp_FOUND)
  message(STATUS "Building Assimp")
  set(_address "https://github.com/assimp/assimp.git")
  message(STATUS "fetch: ${_address}")
  FetchContent_Declare(
    assimp
    GIT_REPOSITORY ${_address}
  )
  FetchContent_MakeAvailable(assimp)
  message(STATUS "Assimp build done")
endif()

set(DONUT_SHADERS_OUTPUT_DIR "${PROJECT_SOURCE_DIR}/bin/shaders/framework")

message(STATUS "find package: donut")
find_package(donut QUIET)
if(NOT DONUT_FOUND)
  message(STATUS "Building donut")
  set(_address "https://github.com/NVIDIAGameWorks/donut.git")
  message(STATUS "fetch: ${_address}")
  FetchContent_Declare(
    donut
    GIT_REPOSITORY ${_address}
    GIT_TAG main
  )

  FetchContent_MakeAvailable(donut)
  message(STATUS "donut build done")
endif()

message("eigen3 is at ${Eigen3_SOURCE_DIR}")
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${Eigen3_SOURCE_DIR})
include_directories(${Vulkan_INCLUDE_DIRS})

Ubpa_AddSubDirsRec(src)
