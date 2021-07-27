cmake_minimum_required(VERSION 3.17)
project(Numerical VERSION 1.0 LANGUAGES CXX CUDA C)

set(CMAKE_CXX_STANDARD 17)

include(cmake/InitUCMake.cmake)

Ubpa_InitUCMake()
Ubpa_InitProject()


message(STATUS "find package: glm")
find_package(glm QUIET)
if(NOT GLM_FOUND)
  message(STATUS "Building glm")
  set(_address "https://github.com/Jerry-Shen0527/glm.git")
  message(STATUS "fetch: ${_address}")
  FetchContent_Declare(
    glm
    GIT_REPOSITORY ${_address}
  )
  FetchContent_MakeAvailable(glm)
  message(STATUS "glm build done")
endif()

message("glm is at ${glm_SOURCE_DIR}")
include_directories(${glm_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)

Ubpa_AddSubDirsRec(src)
