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

add_executable(dump_results
    dump_results.cpp
)

target_link_libraries(dump_results
    PRIVATE
        thermal_conductivity_evaluators
)

include(GoogleTest)

gtest_discover_tests(thermal_conductivity_tests)