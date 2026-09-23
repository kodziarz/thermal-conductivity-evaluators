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

    std::vector<simulation_value_t> SequentialHeatSimulation::evaluateGeneration(
        const simulation_value_t *k_values,
        const simulation_value_t *invC_values,
        const simulation_value_t *qGen_values,
        int individualsNumber,
        simulation_value_t *minFinalTemperatures,
        simulation_steps_index_t *lastEquilibriumMoment)
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

        int boardSize = boardHeight * boardWidth;
        for (int ind = 0; ind < individualsNumber; ind++)
        {
            const simulation_value_t *indK = k_values + ind * boardSize;
            const simulation_value_t *indInvC = invC_values + ind * boardSize;
            const simulation_value_t *indQGen = qGen_values + ind * boardSize;

            simulation_value_t *finalTemperatures;
            simulation_steps_index_t equilibriumMoment = 0;
            results.push_back(evaluateSystemLayout(indK, indInvC, indQGen, &finalTemperatures, &equilibriumMoment));

            if (minFinalTemperatures != NULL)
            {
                for (int i = 0; i < boardSize; i++)
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

    simulation_value_t SequentialHeatSimulation::evaluateSystemLayout(
        const simulation_value_t *k,
        const simulation_value_t *invC,
        const simulation_value_t *qGen,
        simulation_value_t **returnedFinalTemperatures,
        simulation_steps_index_t *returnedEquilibriumMoment)
    {
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

                    simulation_value_t alpha_current = k[cellIndex] * invC[cellIndex];
                    simulation_value_t beta = qGen[cellIndex] * invC[cellIndex];

                    // order of neighbors: bottom, up, right, left
                    simulation_value_t alpha_neighbor;
                    alpha_neighbor = k[cellIndex + boardWidth] * invC[cellIndex + boardWidth];
                    flow += (inputTs[cellIndex + boardWidth] - currentT) *
                            ((k[cellIndex] == 0 || k[cellIndex + boardWidth] == 0) ? 0 : (alpha_current + alpha_neighbor) / 2);

                    alpha_neighbor = k[cellIndex - boardWidth] * invC[cellIndex - boardWidth];
                    flow += (inputTs[cellIndex - boardWidth] - currentT) *
                            ((k[cellIndex] == 0 || k[cellIndex - boardWidth] == 0) ? 0 : (alpha_current + alpha_neighbor) / 2);

                    alpha_neighbor = k[cellIndex + 1] * invC[cellIndex + 1];
                    flow += (inputTs[cellIndex + 1] - currentT) *
                            ((k[cellIndex] == 0 || k[cellIndex + 1] == 0) ? 0 : (alpha_current + alpha_neighbor) / 2);

                    alpha_neighbor = k[cellIndex - 1] * invC[cellIndex - 1];
                    flow += (inputTs[cellIndex - 1] - currentT) *
                            ((k[cellIndex] == 0 || k[cellIndex - 1] == 0) ? 0 : (alpha_current + alpha_neighbor) / 2);

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

            // Reset DRAIN cells (invC == 0) to drainTemperature
            for (int idx = 0; idx < boardHeight * boardWidth; idx++)
            {
                if (invC[idx] == 0)
                {
                    outputTs[idx] = drainTemperature;
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

        return -maxT;
    }

}
