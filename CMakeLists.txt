cmake_minimum_required(VERSION 3.30)
project(conductivity_evaluators LANGUAGES CXX)

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

option(CONDUCTIVITY_EVALUATORS_BUILD_TESTS "Build conductivity_evaluators tests" OFF)
if(CONDUCTIVITY_EVALUATORS_BUILD_TESTS)
	enable_testing()
	add_subdirectory(tests)
endif()
