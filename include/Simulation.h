#ifndef SIMULATION_H
#define SIMULATION_H
#include <chrono>
#include <algorithm>
#include <cstring>
#include "GraphicsUtils.h"
#include <vector>

#define SIMULATION_STEPS 150'000

namespace conductivity_evaluators
{

    namespace DefaultSimulationParams
    {
        const simulation_value_t ETA = 0.0000001;
        const int numTimeSteps = SIMULATION_STEPS;

        const simulation_value_t startTemperature = 300.0;
        const simulation_value_t drainTemperature = 280.0;

        const simulation_value_t delta_time = 0.01;
    }

    typedef struct
    {
        simulation_steps_index_t simulationSteps = SIMULATION_STEPS;
        const simulation_value_t *startTemperatures = NULL;
        simulation_value_t drainTemperature = DefaultSimulationParams::drainTemperature;
        simulation_value_t delta_time = DefaultSimulationParams::delta_time;
        simulation_value_t ETA = DefaultSimulationParams::ETA;
        const simulation_value_t *k_values = NULL;
        const simulation_value_t *invC_values = NULL;
        const simulation_value_t *qGen_values = NULL;
    } SimulationParams;

    class Simulation
    {
    public:
        Simulation(int boardHeight, int boardWidth,
                   const SimulationParams &params = SimulationParams{}) : boardHeight(boardHeight), boardWidth(boardWidth), individualHeight(boardHeight - 2), individualWidth(boardWidth - 2), simulationSteps(params.simulationSteps), drainTemperature(params.drainTemperature), delta_time(params.delta_time), ETA(params.ETA)
        {
            int boardSize = boardWidth * boardHeight;
            startTemperatures = new simulation_value_t[boardSize]{0};
            if (params.startTemperatures == nullptr)
                std::fill_n(startTemperatures, boardSize, params.drainTemperature);
            else
                std::memcpy(startTemperatures, params.startTemperatures, boardSize * sizeof(simulation_value_t));

            k_values = new simulation_value_t[boardSize]{0};
            invC_values = new simulation_value_t[boardSize]{0};
            qGen_values = new simulation_value_t[boardSize]{0};

            if (params.k_values != nullptr)
                std::memcpy(k_values, params.k_values, boardSize * sizeof(simulation_value_t));
            if (params.invC_values != nullptr)
                std::memcpy(invC_values, params.invC_values, boardSize * sizeof(simulation_value_t));
            if (params.qGen_values != nullptr)
                std::memcpy(qGen_values, params.qGen_values, boardSize * sizeof(simulation_value_t));
        };

        ~Simulation()
        {
            delete[] startTemperatures;
            delete[] k_values;
            delete[] invC_values;
            delete[] qGen_values;
        }

        virtual std::vector<simulation_value_t> evaluateGeneration(
            const simulation_value_t *k_values,
            const simulation_value_t *invC_values,
            const simulation_value_t *qGen_values,
            int individualsNumber,
            simulation_value_t *minFinalTemperatures = NULL,
            simulation_steps_index_t *lastEquilibriumMoment = NULL) = 0;

        virtual void setSimulationParams(const SimulationParams &params);

    protected:
        const int boardHeight, boardWidth, individualWidth, individualHeight;
        simulation_steps_index_t simulationSteps;
        simulation_value_t *startTemperatures, *k_values, *invC_values, *qGen_values;
        simulation_value_t drainTemperature, delta_time, ETA;
    };

}

using conductivity_evaluators::Simulation;
using conductivity_evaluators::SimulationParams;

#endif // SIMULATION_H
