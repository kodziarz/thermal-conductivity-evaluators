#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <chrono>
#include "SequentialHeatSimulation.h"
#include "env.h"

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

    std::vector<simulation_value_t> SequentialHeatSimulation::evaluateGeneration(const std::vector<cell_type_t> &fenotypes, simulation_value_t *finalTemperatureDistributions, simulation_steps_index_t *equilibriumMoments)
    {
#ifdef BENCHMARK
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
#endif
        int boardSize = boardHeight * boardWidth * boardThickness;

        std::vector<simulation_value_t> results;
        if (finalTemperatureDistributions != NULL)
        {
            std::fill_n(finalTemperatureDistributions, boardSize, -1);
        }

        if (equilibriumMoments != NULL)
        {
            *std::fill_n(equilibriumMoments, fenotypes.size() / (boardSize), 0);
        }

        int fenotype_i = 0;
        for (const cell_type_t *fenotype = fenotypes.data(); fenotype < fenotypes.data() + fenotypes.size(); fenotype += boardSize)
        {
            simulation_value_t *finalTemperatures;
            simulation_steps_index_t equilibriumMoment = 0;
            results.push_back(evaluateFenotype(fenotype, &finalTemperatures, &equilibriumMoment));

            if (finalTemperatures != NULL)
            {
                std::memcpy(finalTemperatureDistributions + fenotype_i + boardSize, finalTemperatures, boardSize * sizeof(simulation_value_t));
            }
            delete[] finalTemperatures;

            if (equilibriumMoments != NULL)
            {
                equilibriumMoments[fenotype_i] = equilibriumMoment;
            }
            fenotype_i++;
        }
#ifdef BENCHMARK
        std::chrono::time_point<std::chrono::high_resolution_clock> stop = std::chrono::high_resolution_clock::now();
        std::cout << "BENCHMARK: generation evaluation took "
                  << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
                  << " microseconds" << std::endl;
#endif
        return results;
    }

    simulation_value_t SequentialHeatSimulation::evaluateFenotype(const cell_type_t *fenotype, simulation_value_t **returnedFinalTemperatures, simulation_steps_index_t *returnedEquilibriumMoment)
    {
        int depthStratumSize = boardHeight * boardWidth;
        int boardSize = boardHeight * boardWidth * boardThickness;

        std::vector<simulation_value_t> inputTs = std::vector<simulation_value_t>(boardSize);
        simulation_value_t *outputTs = new simulation_value_t[boardSize];

        // initialize temperature arrays to the starting temperature
        for (int idx = 0; idx < boardSize; ++idx)
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
                    for (int depthStratum = 1; depthStratum < boardThickness - 1; depthStratum++)
                    {

                        const int cellIndex = depthStratum * depthStratumSize + row * boardWidth + col;
                        const simulation_value_t currentT = inputTs[cellIndex];
                        simulation_value_t flow = 0;

                        simulation_value_t beta = GENERATOR_BETA;
                        if (fenotype[cellIndex] != Cell::GENERATOR)
                        {
                            beta = CONDUCTOR_BETA;
                        }
                        // order of neighbors: bottom, up, right, left, below, above
                        flow += (inputTs[cellIndex + boardWidth] - currentT) *
                                calculateMutualAlpha(fenotype[cellIndex], fenotype[cellIndex + boardWidth]);

                        flow += (inputTs[cellIndex - boardWidth] - currentT) *
                                calculateMutualAlpha(fenotype[cellIndex], fenotype[cellIndex - boardWidth]);

                        flow += (inputTs[cellIndex + 1] - currentT) *
                                calculateMutualAlpha(fenotype[cellIndex], fenotype[cellIndex + 1]);

                        flow += (inputTs[cellIndex - 1] - currentT) *
                                calculateMutualAlpha(fenotype[cellIndex], fenotype[cellIndex - 1]);

                        flow += (inputTs[cellIndex + depthStratumSize] - currentT) *
                                calculateMutualAlpha(fenotype[cellIndex], fenotype[cellIndex + depthStratumSize]);

                        flow += (inputTs[cellIndex - depthStratumSize] - currentT) *
                                calculateMutualAlpha(fenotype[cellIndex], fenotype[cellIndex - depthStratumSize]);

                        simulation_value_t resultantPower = flow + beta;

                        simulation_value_t temperatureIncrease = delta_time * resultantPower;
                        outputTs[cellIndex] += temperatureIncrease;

                        if (returnedEquilibriumMoment != NULL && resultantPower < RESULTANT_POWER_TOL)
                        {
                            *returnedEquilibriumMoment = i;
                        }
                        maxT = std::max(maxT, outputTs[cellIndex]);
                    }
                }
            }
            std::memcpy(inputTs.data(), outputTs, boardSize * sizeof(simulation_value_t));
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