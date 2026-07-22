cmake_minimum_required(VERSION 3.30)

project(thermal_conductivity_evaluators)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenCL REQUIRED)

# embed kernels into source code

set(GENERATED_DIR ${CMAKE_BINARY_DIR}/generated)

file(MAKE_DIRECTORY ${GENERATED_DIR})

function(embed_kernel INPUT_FILE)

    get_filename_component(FILE_NAME ${INPUT_FILE} NAME_WE)

    set(OUTPUT_FILE
        ${GENERATED_DIR}/${FILE_NAME}.hpp
    )

    add_custom_command(
        OUTPUT ${OUTPUT_FILE}
        COMMAND
            ${CMAKE_COMMAND}
            -DINPUT=${INPUT_FILE}
            -DOUTPUT=${OUTPUT_FILE}
            -DVARIABLE=${FILE_NAME}
            -P ${PROJECT_SOURCE_DIR}/cmake/embed_kernel.cmake
        DEPENDS ${INPUT_FILE}
    )

    list(APPEND GENERATED_KERNEL_HEADERS ${OUTPUT_FILE})
    set(GENERATED_KERNEL_HEADERS
        ${GENERATED_KERNEL_HEADERS}
        PARENT_SCOPE
    )

endfunction()

file(GLOB KERNEL_FILES
    CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/kernels/*.cl"
)

foreach(KERNEL_FILE IN LISTS KERNEL_FILES)
    embed_kernel(${KERNEL_FILE})
endforeach()

# The main configuration

add_library(${PROJECT_NAME} STATIC
	src/GraphicsUtils.cpp
	src/ParallelHeatSimulation.cpp
	src/SequentialHeatSimulation.cpp
	src/Simulation.cpp

	${GENERATED_KERNEL_HEADERS}
)

target_include_directories(${PROJECT_NAME}
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/include
	PRIVATE
        ${GENERATED_DIR}
)

target_link_libraries(${PROJECT_NAME}
	PUBLIC
		OpenCL::OpenCL
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

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
