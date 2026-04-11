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
        const simulation_value_t RESULTANT_POWER_TOL = 1e-4;
        const int numTimeSteps = SIMULATION_STEPS;

        const simulation_value_t startTemperature = 300.0;
        const simulation_value_t drainTemperature = 280.0;

        const simulation_value_t GENERATOR_ALPHA = 2;
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
        simulation_value_t RESULTANT_POWER_TOL = DefaultSimulationParams::RESULTANT_POWER_TOL;
    } SimulationParams;

    class Simulation
    {
    public:
        Simulation(int boardHeight, int boardWidth, int thickness,
                   const SimulationParams &params = SimulationParams{}) : boardHeight(boardHeight), boardWidth(boardWidth), boardThickness(thickness), individualHeight(boardHeight - 2), individualWidth(boardWidth - 2), individualThickness(thickness - 2), simulationSteps(params.simulationSteps), drainTemperature(params.drainTemperature), delta_time(params.delta_time), DRAIN_ALPHA(params.DRAIN_ALPHA), CONDUCTOR_ALPHA(params.CONDUCTOR_ALPHA), GENERATOR_ALPHA(params.GENERATOR_ALPHA), CONDUCTOR_BETA(params.CONDUCTOR_BETA), GENERATOR_BETA(params.GENERATOR_BETA), ETA(params.ETA), RESULTANT_POWER_TOL(params.RESULTANT_POWER_TOL)
        {
            startTemperatures = new simulation_value_t[boardWidth * boardHeight * boardThickness]{0};
            if (params.startTemperatures == nullptr)
                std::fill_n(startTemperatures, boardHeight * boardWidth, params.drainTemperature);
            else
                std::memcpy(startTemperatures, params.startTemperatures, boardHeight * boardWidth * sizeof(simulation_value_t));
        };

        ~Simulation()
        {
            delete[] startTemperatures;
        }

        virtual std::vector<simulation_value_t> evaluateGeneration(const std::vector<cell_type_t> &fenotypes, simulation_value_t *finalTemperatureDistributions = NULL, simulation_steps_index_t *equilibriumMoments = NULL) = 0;

        virtual void setSimulationParams(const SimulationParams &params);

        virtual void getSimulationParams(SimulationParams &params)
        {
            params.simulationSteps = simulationSteps;
            params.drainTemperature = drainTemperature;
            params.delta_time = delta_time;
            params.DRAIN_ALPHA = DRAIN_ALPHA;
            params.CONDUCTOR_ALPHA = CONDUCTOR_ALPHA;
            params.GENERATOR_ALPHA = GENERATOR_ALPHA;
            params.CONDUCTOR_BETA = CONDUCTOR_BETA;
            params.GENERATOR_BETA = GENERATOR_BETA;
            params.ETA = ETA;
            params.RESULTANT_POWER_TOL = RESULTANT_POWER_TOL;

            // Deep copy of temperature differences
            if (params.startTemperatures == NULL)
            {
                params.startTemperatures = new simulation_value_t[boardHeight * boardWidth * boardThickness];
            }
            std::memcpy((void *)params.startTemperatures, startTemperatures, boardHeight * boardWidth * boardThickness * sizeof(simulation_value_t));
        };

    protected:
        const int boardHeight, boardWidth, boardThickness, individualWidth, individualHeight, individualThickness;
        simulation_steps_index_t simulationSteps;
        simulation_value_t *startTemperatures, drainTemperature, delta_time, DRAIN_ALPHA, CONDUCTOR_ALPHA, GENERATOR_ALPHA, CONDUCTOR_BETA, GENERATOR_BETA, ETA, RESULTANT_POWER_TOL;
    };

}

using conductivity_evaluators::Simulation;
using conductivity_evaluators::SimulationParams;

#endif // SIMULATION_H
