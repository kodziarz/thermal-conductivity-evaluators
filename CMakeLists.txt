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
