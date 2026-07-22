#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <chrono>
#include "SequentialHeatSimulation.h"

#if defined(__has_include)
#if __has_include("env.h")
#include "env.h"
#endif
#else
#if defined(__GNUC__) || defined(__clang__)
#warning "__has_include is not supported"
#elif defined(_MSC_VER)
#pragma message("__has_include is not supported")
#endif
#endif

namespace conductivity_evaluators
{

    void SequentialHeatSimulation::save_heat_data(int height, int width, simulation_value_t **board, simulation_value_t maxTemp, simulation_steps_index_t timestep)
    {
        std::string filename = "../heat/timestep" + std::to_string(timestep) + ".csv";
        std::ofstream out(filename);

        if (!out)
        {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        out << "maxTemp," << maxTemp << "\n";

        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                out << board[i][j];
                if (j + 1 < width)
                    out << ",";
            }
            out << "\n";
        }
    }

    std::vector<simulation_value_t> SequentialHeatSimulation::evaluateGeneration(const std::vector<cell_type_t> &systemLayouts, simulation_value_t *minFinalTemperatures, simulation_steps_index_t *lastEquilibriumMoment)
    {
#ifdef BENCHMARK
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
#endif
        std::vector<simulation_value_t> results;
        if (minFinalTemperatures != NULL)
        {
            // initialize temperatures
            for (int i = 0; i < boardHeight * boardWidth; i++)
            {
                minFinalTemperatures[i] = std::numeric_limits<simulation_value_t>::infinity();
            }
        }

        if (lastEquilibriumMoment != NULL)
        {
            *lastEquilibriumMoment = 0;
        }

        for (const cell_type_t *systemLayout = systemLayouts.data(); systemLayout < systemLayouts.data() + systemLayouts.size(); systemLayout += boardHeight * boardWidth)
        {
            simulation_value_t *finalTemperatures;
            simulation_steps_index_t equilibriumMoment = 0;
            results.push_back(evaluateSystemLayout(systemLayout, &finalTemperatures, &equilibriumMoment));

            if (minFinalTemperatures != NULL)
            {
                for (int i = 0; i < boardHeight * boardWidth; i++)
                {
                    minFinalTemperatures[i] = std::min(minFinalTemperatures[i], finalTemperatures[i]);
                }
            }
            delete[] finalTemperatures;

            if (lastEquilibriumMoment != NULL)
            {
                *lastEquilibriumMoment = std::max(*lastEquilibriumMoment, equilibriumMoment);
            }
        }
#ifdef BENCHMARK
        std::chrono::time_point<std::chrono::high_resolution_clock> stop = std::chrono::high_resolution_clock::now();
        std::cout << "BENCHMARK: generation evaluation took "
                  << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
                  << " microseconds" << std::endl;
#endif
        return results;
    }

    simulation_value_t SequentialHeatSimulation::evaluateSystemLayout(const cell_type_t *systemLayout, simulation_value_t **returnedFinalTemperatures, simulation_steps_index_t *returnedEquilibriumMoment)
    {
        // simulation_value_t *inputTs = new simulation_value_t[boardHeight * boardWidth];
        std::vector<simulation_value_t> inputTs = std::vector<simulation_value_t>(boardHeight * boardWidth);
        simulation_value_t *outputTs = new simulation_value_t[boardHeight * boardWidth];

        // initialize temperature arrays to the starting temperature
        for (int idx = 0; idx < boardHeight * boardWidth; ++idx)
        {
            inputTs[idx] = startTemperatures[idx];
            outputTs[idx] = startTemperatures[idx];
        }

        simulation_value_t maxT = inputTs[0];

        for (int i = 0; i < simulationSteps; i++)
        {
            for (int row = 1; row < boardHeight - 1; row++)
            {
                for (int col = 1; col < boardWidth - 1; col++)
                {
                    const int cellIndex = row * boardWidth + col;
                    const simulation_value_t currentT = inputTs[cellIndex];
                    simulation_value_t flow = 0;

                    simulation_value_t beta = GENERATOR_BETA;
                    if (systemLayout[cellIndex] != Cell::GENERATOR)
                    {
                        beta = CONDUCTOR_BETA;
                    }
                    // order of neighbors: bottom, up, right, left
                    flow += (inputTs[cellIndex + boardWidth] - currentT) *
                            calculateMutualAlpha(systemLayout[cellIndex], systemLayout[cellIndex + boardWidth]);

                    flow += (inputTs[cellIndex - boardWidth] - currentT) *
                            calculateMutualAlpha(systemLayout[cellIndex], systemLayout[cellIndex - boardWidth]);

                    flow += (inputTs[cellIndex + 1] - currentT) *
                            calculateMutualAlpha(systemLayout[cellIndex], systemLayout[cellIndex + 1]);

                    flow += (inputTs[cellIndex - 1] - currentT) *
                            calculateMutualAlpha(systemLayout[cellIndex], systemLayout[cellIndex - 1]);

                    simulation_value_t temperatureIncrease = delta_time * (flow + beta);
                    outputTs[cellIndex] += temperatureIncrease;

                    simulation_value_t errorDenominator = std::max(
                        std::fabs(inputTs[cellIndex]),
                        std::fabs(outputTs[cellIndex]));
                    if (returnedEquilibriumMoment != NULL && errorDenominator != 0 && std::fabs(temperatureIncrease / errorDenominator) >= ETA)
                    {
                        *returnedEquilibriumMoment = i;
                    }
                    maxT = std::max(maxT, outputTs[cellIndex]);
                }
            }
            std::memcpy(inputTs.data(), outputTs, boardHeight * boardWidth * sizeof(simulation_value_t));
        }

        if (returnedFinalTemperatures != NULL)
        {
            *returnedFinalTemperatures = outputTs;
        }
        else
        {
            delete[] outputTs;
        }
        // delete[] inputTs;

        return -maxT;
    }

    simulation_value_t SequentialHeatSimulation::calculateMutualAlpha(cell_type_t considered_cell_type, cell_type_t neighbor_type)
    {
        simulation_value_t considered_cell_alpha;
        switch (considered_cell_type)
        {
        case Cell::DRAIN:
        case Cell::ADIABATIC:
            return 0;
            break;
        case Cell::GENERATOR:
            considered_cell_alpha = GENERATOR_ALPHA;
            break;
        case Cell::CONDUCTOR:
            considered_cell_alpha = CONDUCTOR_ALPHA;
            break;
        default:
            throw std::runtime_error(std::string("calculateMutualAlpha(): unexpected considered_cell_alpha type value was received: ") + std::to_string(considered_cell_type));
            break;
        }
        switch (neighbor_type)
        {
        case Cell::ADIABATIC:
            return 0;
            break;
        case Cell::DRAIN:
            return (considered_cell_alpha + DRAIN_ALPHA) / 2;
            break;
        case Cell::GENERATOR:
            return (considered_cell_alpha + GENERATOR_ALPHA) / 2;
            break;
        case Cell::CONDUCTOR:
            return (considered_cell_alpha + CONDUCTOR_ALPHA) / 2;
            break;
        default:
            throw std::runtime_error(std::string("calculateMutualAlpha(): unexpected neighbor_type type value was received: ") + std::to_string(neighbor_type));
            break;
        }
    }

}