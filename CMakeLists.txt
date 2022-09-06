cmake_minimum_required(VERSION 3.17)
project(Numerical VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

option(${PROJECT_NAME}_ENABLE_UNIT_TESTING on)

if(${PROJECT_NAME}_ENABLE_UNIT_TESTING)
  enable_testing()
  message(STATUS "Build unit tests for the project. Tests should always be found in the test folder\n")
  add_subdirectory(ext/googletest)
  set(${CMAKE_PROJECT_NAME}_USE_GTEST on)
  add_subdirectory(src/tests)
endif()
