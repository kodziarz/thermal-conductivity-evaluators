include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)

set(CONDUCTIVITY_EVALUATORS_TEST_COMMON_INCLUDES
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
)

add_compile_options(-O0 -g)

add_executable(test_simulation test_simulation.cpp ../src/SequentialHeatSimulation.cpp ../src/ParallelHeatSimulation.cpp ../src/GraphicsUtils.cpp)
target_include_directories(test_simulation PRIVATE ${CONDUCTIVITY_EVALUATORS_TEST_COMMON_INCLUDES})
target_link_libraries(test_simulation PRIVATE conductivity_evaluators  GTest::gtest_main)
add_test(NAME test_simulation COMMAND test_simulation)

add_executable(test_parallel_simulation test_parallel_simulation.cpp)
target_include_directories(test_parallel_simulation PRIVATE ${CONDUCTIVITY_EVALUATORS_TEST_COMMON_INCLUDES})
target_link_libraries(test_parallel_simulation PRIVATE conductivity_evaluators GTest::gtest_main)
add_test(NAME test_parallel_simulation COMMAND test_parallel_simulation)

# add_executable(test_simulation_comparison test_simulation_comparison.cpp)
# target_include_directories(test_simulation_comparison PRIVATE ${CONDUCTIVITY_EVALUATORS_TEST_COMMON_INCLUDES})
# target_link_libraries(test_simulation_comparison PRIVATE conductivity_evaluators GTest::gtest)
# add_test(NAME test_simulation_comparison COMMAND test_simulation_comparison)
