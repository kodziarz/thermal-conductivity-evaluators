#ifndef SEQUNTIALHEATSIMULATION_H
#define SEQUNTIALHEATSIMULATION_H

#include <cstring>
#include <string>
#include <vector>
#include "Simulation.h"
#include "Fenotype.h"

namespace conductivity_evaluators
{

    /*
    Possible optimizations:
    - better array structure
    - array with alphas and betas instead of if in the loop
    - remove double calculations on borders (from both neighbors perspectives)
    */

    const std::string dir = "../heat";

    class SequentialHeatSimulation : public Simulation
    {
        void static save_heat_data(int height, int width, simulation_value_t **board, simulation_value_t maxTemp, simulation_steps_index_t timestep);

    public:
        SequentialHeatSimulation(int boardHeight, int boardWidth, const SimulationParams &params = SimulationParams{}) : Simulation(boardHeight, boardWidth, params) {};

        simulation_value_t evaluateFenotype(const cell_type_t *fenotype, simulation_value_t **returnedFinalTemperatures = NULL, simulation_steps_index_t *returnedEquilibriumMoment = NULL);
        std::vector<simulation_value_t> evaluateGeneration(const std::vector<cell_type_t> &fenotypes, simulation_value_t *minFinalTemperatures = NULL, simulation_steps_index_t *lastEquilibriumMoment = NULL) override;

    protected:
        simulation_value_t calculateMutualAlpha(cell_type_t neighbor_type, cell_type_t considered_cell_type);
    };

}

using conductivity_evaluators::SequentialHeatSimulation;

#endif // SEQUNTIALHEATSIMULATION_H
