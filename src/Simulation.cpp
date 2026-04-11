#include "Simulation.h"

namespace conductivity_evaluators
{

    void Simulation::setSimulationParams(const SimulationParams &params)
    {
        simulationSteps = params.simulationSteps;
        drainTemperature = params.drainTemperature;
        delta_time = params.delta_time;
        DRAIN_ALPHA = params.DRAIN_ALPHA;
        CONDUCTOR_ALPHA = params.CONDUCTOR_ALPHA;
        GENERATOR_ALPHA = params.GENERATOR_ALPHA;
        CONDUCTOR_BETA = params.CONDUCTOR_BETA;
        GENERATOR_BETA = params.GENERATOR_BETA;
        ETA = params.ETA;
        RESULTANT_POWER_TOL = params.RESULTANT_POWER_TOL;

        if (params.startTemperatures != NULL)
        {
            std::memcpy((void *)startTemperatures, params.startTemperatures, boardHeight * boardWidth * boardThickness * sizeof(simulation_value_t));
        }
    }

}