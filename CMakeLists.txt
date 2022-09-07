cmake_minimum_required(VERSION 3.17)
project(Numerical VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

option(${PROJECT_NAME}_ENABLE_UNIT_TESTING on)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)


file(GLOB Numerical_sources ${CMAKE_CURRENT_SOURCE_DIR}/src/Numerical/*.cpp)
message("sources : ${jit_sources}")

add_library(${PROJECT_NAME} ${Numerical_sources})

add_subdirectory(ext/drjit)
target_link_libraries(${PROJECT_NAME} drjit drjit-core)

set(OUTPUT_DIR 
  RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}"
  RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}"
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}"
)

set_target_properties(
  nanothread
  PROPERTIES
  ${OUTPUT_DIR}
)

set_target_properties(
  drjit-core
  PROPERTIES
  ${OUTPUT_DIR}
)



if(${PROJECT_NAME}_ENABLE_UNIT_TESTING)
  enable_testing()
  message(STATUS "Build unit tests for the project. Tests should always be found in the test folder\n")
  add_subdirectory(ext/googletest)
  set(${CMAKE_PROJECT_NAME}_USE_GTEST on)
  add_subdirectory(tests)
endif()
