cmake_minimum_required(VERSION 3.30)
project(conductivity_evaluators LANGUAGES CXX)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

find_package(OpenCL REQUIRED)
add_definitions(-DCL_TARGET_OPENCL_VERSION=300)

add_library(conductivity_evaluators STATIC
	src/GraphicsUtils.cpp
	src/ParallelHeatSimulation.cpp
	src/SequentialHeatSimulation.cpp
	src/Simulation.cpp
)

target_include_directories(conductivity_evaluators
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(conductivity_evaluators
	PUBLIC
		OpenCL::OpenCL
)

target_compile_features(conductivity_evaluators PUBLIC cxx_std_20)

# Copy kernel file to build directory
add_custom_command(
  TARGET conductivity_evaluators  # Tie the command to the "MyApp" target
  POST_BUILD    # Run after the target is built
  COMMAND ${CMAKE_COMMAND} -E copy_directory  # CMake's cross-platform copy command
		${CMAKE_CURRENT_SOURCE_DIR}/kernels/       # Source directory
    ${CMAKE_CURRENT_BINARY_DIR}/kernels/       # Destination directory
  COMMENT "Copying config files to build directory"  # Optional: Log message
)

option(CONDUCTIVITY_EVALUATORS_BUILD_TESTS "Build conductivity_evaluators tests" OFF)
if(CONDUCTIVITY_EVALUATORS_BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()
