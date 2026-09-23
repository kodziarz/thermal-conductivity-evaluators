#include "Simulation.h"

namespace conductivity_evaluators
{

    void Simulation::setSimulationParams(const SimulationParams &params)
    {
        simulationSteps = params.simulationSteps;
        drainTemperature = params.drainTemperature;
        delta_time = params.delta_time;
        ETA = params.ETA;

        if (params.startTemperatures != NULL)
        {
            for (int i = 0; i < boardHeight * boardWidth; i++)
            {
                startTemperatures[i] = params.startTemperatures[i];
            }
        }

        int boardSize = boardHeight * boardWidth;
        if (params.k_values != NULL)
        {
            std::memcpy(k_values, params.k_values, boardSize * sizeof(simulation_value_t));
        }
        if (params.invC_values != NULL)
        {
            std::memcpy(invC_values, params.invC_values, boardSize * sizeof(simulation_value_t));
        }
        if (params.qGen_values != NULL)
        {
            std::memcpy(qGen_values, params.qGen_values, boardSize * sizeof(simulation_value_t));
        }
    }

}
