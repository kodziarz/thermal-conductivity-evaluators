#include <chrono>

namespace conductivity_evaluators
{
    using simulation_value_t = double;
    using simulation_steps_index_t = int;
    using cell_type_t = int;
    using timestamp = std::chrono::time_point<std::chrono::high_resolution_clock>;
}

using conductivity_evaluators::cell_type_t;
using conductivity_evaluators::simulation_steps_index_t;
using conductivity_evaluators::simulation_value_t;
using conductivity_evaluators::timestamp;
