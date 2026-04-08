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
        // Unrealistic Parameters (But they bring more better results)
        const simulation_value_t ETA = 0.0000001;
        const int numTimeSteps = SIMULATION_STEPS;

        const simulation_value_t startTemperature = 300.0;
        const simulation_value_t drainTemperature = 280.0;

        const simulation_value_t GENERATOR_ALPHA = 0.1;
        const simulation_value_t GENERATOR_BETA = 10;

        const simulation_value_t CONDUCTOR_ALPHA = 10;
        const simulation_value_t CONDUCTOR_BETA = 0;

        const simulation_value_t delta_time = 0.01;

        const simulation_value_t DRAIN_ALPHA = CONDUCTOR_ALPHA;

    }

    typedef struct
    {
        simulation_steps_index_t simulationSteps = SIMULATION_STEPS;
        const simulation_value_t *startTemperatures = NULL;
        simulation_value_t drainTemperature = DefaultSimulationParams::drainTemperature;
        simulation_value_t delta_time = DefaultSimulationParams::delta_time;
        simulation_value_t DRAIN_ALPHA = DefaultSimulationParams::DRAIN_ALPHA;
        simulation_value_t CONDUCTOR_ALPHA = DefaultSimulationParams::CONDUCTOR_ALPHA;
        simulation_value_t GENERATOR_ALPHA = DefaultSimulationParams::GENERATOR_ALPHA;
        simulation_value_t CONDUCTOR_BETA = DefaultSimulationParams::CONDUCTOR_BETA;
        simulation_value_t GENERATOR_BETA = DefaultSimulationParams::GENERATOR_BETA;
        simulation_value_t ETA = DefaultSimulationParams::ETA;
    } SimulationParams;

    class Simulation
    {
    public:
        Simulation(int boardHeight, int boardWidth,
                   const SimulationParams &params = SimulationParams{}) : boardHeight(boardHeight), boardWidth(boardWidth), individualHeight(boardHeight - 2), individualWidth(boardWidth - 2), simulationSteps(params.simulationSteps), drainTemperature(params.drainTemperature), delta_time(params.delta_time), DRAIN_ALPHA(params.DRAIN_ALPHA), CONDUCTOR_ALPHA(params.CONDUCTOR_ALPHA), GENERATOR_ALPHA(params.GENERATOR_ALPHA), CONDUCTOR_BETA(params.CONDUCTOR_BETA), GENERATOR_BETA(params.GENERATOR_BETA), ETA(params.ETA)
        {
            startTemperatures = new simulation_value_t[boardWidth * boardHeight]{0};
            if (params.startTemperatures == nullptr)
                std::fill_n(startTemperatures, boardHeight * boardWidth, params.drainTemperature);
            else
                std::memcpy(startTemperatures, params.startTemperatures, boardHeight * boardWidth * sizeof(simulation_value_t));
        };

        ~Simulation()
        {
            delete[] startTemperatures;
        }

        virtual std::vector<simulation_value_t> evaluateGeneration(const std::vector<cell_type_t> &fenotypes, simulation_value_t *minFinalTemperatures = NULL, simulation_steps_index_t *lastEquilibriumMoment = NULL) = 0;

        virtual void setSimulationParams(const SimulationParams &params);

    protected:
        const int boardHeight, boardWidth, individualWidth, individualHeight;
        simulation_steps_index_t simulationSteps;
        simulation_value_t *startTemperatures, drainTemperature, delta_time, DRAIN_ALPHA, CONDUCTOR_ALPHA, GENERATOR_ALPHA, CONDUCTOR_BETA, GENERATOR_BETA, ETA;
    };

}

using conductivity_evaluators::Simulation;
using conductivity_evaluators::SimulationParams;

#endif // SIMULATION_H
