add_executable(thermal_conductivity_tests
    test_parallel_simulation.cpp
    test_simulation_comparison.cpp
    test_simulation.cpp
)

target_link_libraries(thermal_conductivity_tests
    PRIVATE
        thermal_conductivity_evaluators
        GTest::gtest_main
)

include(GoogleTest)

gtest_discover_tests(thermal_conductivity_tests)