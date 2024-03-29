cmake_minimum_required(VERSION 3.18)

#
# Project details
#
 
project(
  ${CMAKE_PROJECT_NAME}_Tests
  LANGUAGES CXX
)

message("Adding tests under ${CMAKE_PROJECT_NAME}_Tests...")
file(GLOB test_sources ${CMAKE_CURRENT_SOURCE_DIR}/ "*.cpp")

foreach(file ${test_sources})
  string(REGEX REPLACE "(.*/)([a-zA-Z0-9_ ]+)(\.cpp)" "\\2" test_name ${file}) 
  message("file name: ${test_name}")
  add_executable(${test_name}_test ${file})

  #
  # Set the compiler standard
  #
  
  target_compile_features(${test_name}_test PUBLIC cxx_std_17)

  #
  # Setup code coverage if enabled
  #

  if (${CMAKE_PROJECT_NAME}_ENABLE_CODE_COVERAGE)
    target_compile_options(${CMAKE_PROJECT_NAME} PUBLIC -O0 -g -fprofile-arcs -ftest-coverage)
    target_link_options(${CMAKE_PROJECT_NAME} PUBLIC -fprofile-arcs -ftest-coverage)
    message("Code coverage is enabled and provided with GCC.")
  endif()

  #
  # Load the desired unit testing framework
  #
  # Currently supported: GoogleTest (and GoogleMock), Catch2.

  message("project lib for test: ${CMAKE_PROJECT_NAME}_LIB")

  set(${CMAKE_PROJECT_NAME}_TEST_LIB ${CMAKE_PROJECT_NAME})
  if(${CMAKE_PROJECT_NAME}_USE_GTEST)

    if(${CMAKE_PROJECT_NAME}_USE_GOOGLE_MOCK)
      set(GOOGLE_MOCK_LIBRARIES GTest::gmock GTest::gmock_main)
    endif()

    target_link_libraries(
      ${test_name}_test
      PUBLIC
        GTest::gtest_main
        ${GOOGLE_MOCK_LIBRARIES}
        ${${CMAKE_PROJECT_NAME}_TEST_LIB}
    )
    message("${test_name}_Tests linked with ${${CMAKE_PROJECT_NAME}_TEST_LIB}")
  else()
    message(FATAL_ERROR "Unknown testing library. Please setup your desired unit testing library by using `target_link_libraries`.")  
  endif()

  #
  # Add the unit tests
  #

  add_test(
    NAME
      ${test_name}
    COMMAND
      ${test_name}_test
  )
endforeach()

message("Finished adding unit tests for ${CMAKE_PROJECT_NAME}.")
