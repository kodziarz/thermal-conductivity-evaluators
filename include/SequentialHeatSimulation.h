#ifndef SEQUNTIALHEATSIMULATION_H
#define SEQUNTIALHEATSIMULATION_H

#include <cstring>
#include <string>
#include <vector>
#include "Simulation.h"

namespace conductivity_evaluators
{

    const std::string dir = "../heat";

    class SequentialHeatSimulation : public Simulation
    {
        void static save_heat_data(int height, int width, simulation_value_t **board, simulation_value_t maxTemp, simulation_steps_index_t timestep);

    public:
        SequentialHeatSimulation(int boardHeight, int boardWidth, const SimulationParams &params = SimulationParams{}) : Simulation(boardHeight, boardWidth, params) {};

        simulation_value_t evaluateSystemLayout(
            const simulation_value_t *k,
            const simulation_value_t *invC,
            const simulation_value_t *qGen,
            simulation_value_t **returnedFinalTemperatures = NULL,
            simulation_steps_index_t *returnedEquilibriumMoment = NULL);

        std::vector<simulation_value_t> evaluateGeneration(
            const simulation_value_t *k_values,
            const simulation_value_t *invC_values,
            const simulation_value_t *qGen_values,
            int individualsNumber,
            simulation_value_t *minFinalTemperatures = NULL,
            simulation_steps_index_t *lastEquilibriumMoment = NULL) override;
    };

}

using conductivity_evaluators::SequentialHeatSimulation;

#endif // SEQUNTIALHEATSIMULATION_H
