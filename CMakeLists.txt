add_library(conductivity_evaluators STATIC
	src/GraphicsUtils.cpp
	src/ParallelHeatSimulation.cpp
	src/SequentialHeatSimulation.cpp
	src/Simulation.cpp
)

target_include_directories(conductivity_evaluators
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/include
		${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(conductivity_evaluators
	PUBLIC
		OpenCL::OpenCL
)

target_compile_features(conductivity_evaluators PUBLIC cxx_std_20)
