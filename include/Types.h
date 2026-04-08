#pragma once

#include <chrono>

namespace conductivity_evaluators
{
    using simulation_value_t = double;
    using simulation_steps_index_t = int;
    using cell_type_t = int;
    using timestamp = std::chrono::time_point<std::chrono::high_resolution_clock>;

    enum Cell
    {
        CONDUCTOR = 0,
        GENERATOR = 1,
        ADIABATIC = 2,
        DRAIN = 3
    };
}

using conductivity_evaluators::cell_type_t;
using conductivity_evaluators::simulation_steps_index_t;
using conductivity_evaluators::simulation_value_t;
using conductivity_evaluators::timestamp;
