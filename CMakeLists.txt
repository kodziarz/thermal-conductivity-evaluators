cmake_minimum_required(VERSION 3.30)

project(thermal_conductivity_evaluators)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenCL REQUIRED)

add_library(thermal_conductivity_evaluators STATIC
	src/GraphicsUtils.cpp
	src/ParallelHeatSimulation.cpp
	src/SequentialHeatSimulation.cpp
	src/Simulation.cpp
)

target_include_directories(thermal_conductivity_evaluators
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(thermal_conductivity_evaluators
	PUBLIC
		OpenCL::OpenCL
)

target_compile_features(thermal_conductivity_evaluators PUBLIC cxx_std_20)

target_compile_definitions(${PROJECT_NAME}
 	PUBLIC
		CL_TARGET_OPENCL_VERSION=300)


# Tests
option(BUILD_TESTS "Build tests" OFF)

if(BUILD_TESTS)

	include(FetchContent)

	FetchContent_Declare(
		googletest
		URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
	)

	FetchContent_MakeAvailable(googletest)

    enable_testing()
	add_subdirectory(tests)
endif()
