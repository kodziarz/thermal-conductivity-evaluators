#pragma once

#include <chrono>

namespace conductivity_evaluators
{
    using simulation_value_t = double;
    using simulation_steps_index_t = int;
    using timestamp = std::chrono::time_point<std::chrono::high_resolution_clock>;

    struct CellProperties
    {
        simulation_value_t k;
        simulation_value_t invC;
        simulation_value_t qGen;
    };
}

using conductivity_evaluators::simulation_steps_index_t;
using conductivity_evaluators::simulation_value_t;
using conductivity_evaluators::timestamp;
using conductivity_evaluators::CellProperties;
