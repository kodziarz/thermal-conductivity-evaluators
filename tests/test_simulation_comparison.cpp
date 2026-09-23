#include <gtest/gtest.h>
#include <cmath>
#include "SequentialHeatSimulation.h"
#include "ParallelHeatSimulation.h"
#include "SystemLayout.h"
#include <algorithm>
#include <fstream>

#define TEST_REL_ERROR 0.00001

#define ASSERT_IN_REL_ERROR(expected, checked, eta)             \
    {                                                           \
        ASSERT_GE(upperToleranceBound(expected, eta), checked); \
        ASSERT_LE(lowerToleranceBound(expected, eta), checked); \
    }

constexpr simulation_value_t lowerToleranceBound(simulation_value_t expectedValue, simulation_value_t eta)
{
    return expectedValue - std::fabs(expectedValue) * eta;
}

constexpr simulation_value_t upperToleranceBound(simulation_value_t expectedValue, simulation_value_t eta)
{
    return expectedValue + std::fabs(expectedValue) * eta;
}

static void createPerCellArrays(
    const int *layout, int size,
    simulation_value_t genK, simulation_value_t genInvC, simulation_value_t genQGen,
    simulation_value_t condK, simulation_value_t condInvC, simulation_value_t condQGen,
    simulation_value_t adiaK, simulation_value_t adiaInvC, simulation_value_t adiaQGen,
    simulation_value_t drainK, simulation_value_t drainInvC, simulation_value_t drainQGen,
    std::vector<simulation_value_t> &k_out,
    std::vector<simulation_value_t> &invC_out,
    std::vector<simulation_value_t> &qGen_out)
{
    k_out.resize(size);
    invC_out.resize(size);
    qGen_out.resize(size);

    SystemLayout::MaterialProperties genProps = {genK, genInvC, genQGen};
    SystemLayout::MaterialProperties condProps = {condK, condInvC, condQGen};
    SystemLayout::MaterialProperties adiaProps = {adiaK, adiaInvC, adiaQGen};
    SystemLayout::MaterialProperties drainProps = {drainK, drainInvC, drainQGen};

    SystemLayout::cellTypeLayoutToProperties(layout, size, condProps, genProps, adiaProps, drainProps,
                                            k_out.data(), invC_out.data(), qGen_out.data());
}

#pragma region 6x6BottomDrain
TEST(HeatSimulationTestComparison, EvaluateGenerationWith6x6GeneratorFullBottomDrainSystemLayout)
{

    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEta = TEST_REL_ERROR / 1000;

    // test

    SystemLayout_t fen = SystemLayout::createGeneratorSystemLayout(boardHeight, boardWidth);
    std::vector<int> systemLayouts(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    std::vector<simulation_value_t> kValues, invCValues, qGenValues;
    createPerCellArrays(systemLayouts.data(), boardHeight * boardWidth,
                        10, 1.0, 1,
                        100, 1.0, 0,
                        0, 1.0, 0,
                        100, 0, 0,
                        kValues, invCValues, qGenValues);

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];
    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    SimulationParams testParams = {
        .simulationSteps = 0,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .ETA = usedEta,
        .k_values = kValues.data(),
        .invC_values = invCValues.data(),
        .qGen_values = qGenValues.data()};

    std::ofstream seqFile("seqTsRegister.txt");
    std::ofstream parFile("parTsRegister.txt");

    // Check if the file was opened successfully
    if (seqFile.is_open() && parFile.is_open())
    {

        SequentialHeatSimulation sequentialSimulator(boardHeight, boardWidth, testParams);
        simulation_value_t *sequentialMinTemperatures = new simulation_value_t[boardHeight * boardWidth];
        simulation_steps_index_t sequentialEquilibriumStep;

        ParallelHeatSimulation parallelSimulator(boardHeight, boardWidth, testParams);
        simulation_value_t *parallelMinTemperatures = new simulation_value_t[boardHeight * boardWidth];
        simulation_steps_index_t parallelEquilibriumStep;

        for (int i = 0; i < 1'333; i++)
        {
            testParams.simulationSteps = i;
            sequentialSimulator.setSimulationParams(testParams);
            parallelSimulator.setSimulationParams(testParams);

            auto sequentialResult = sequentialSimulator.evaluateGeneration(kValues.data(), invCValues.data(), qGenValues.data(), 1, sequentialMinTemperatures, &sequentialEquilibriumStep);

            auto parallelResult = parallelSimulator.evaluateGeneration(kValues.data(), invCValues.data(), qGenValues.data(), 1, parallelMinTemperatures, &parallelEquilibriumStep);

            // assert

            ASSERT_EQ(sequentialResult.size(), 1u);
            ASSERT_EQ(parallelResult.size(), 1u);

            seqFile << "Step: " << i << std::endl;
            parFile << "Step: " << i << std::endl;

            for (int row = 0; row < boardHeight; row++)
            {
                for (int column = 0; column < boardWidth; column++)
                {
                    seqFile << std::fixed << std::setprecision(13) << sequentialMinTemperatures[row * boardWidth + column] << " ";
                    parFile << std::fixed << std::setprecision(13) << parallelMinTemperatures[row * boardWidth + column] << " ";
                }
                seqFile << std::endl;
                parFile << std::endl;
            }
            seqFile << std::endl;
            parFile << std::endl;
        }

        delete[] sequentialMinTemperatures;
        delete[] parallelMinTemperatures;

        // Close the file
        seqFile.close();
        parFile.close();

        std::cout << "Simulations' registers written!" << std::endl;
    }
    else
    {
        // If the file couldn't be opened, show an error message
        std::cerr << "Error opening the file!" << std::endl;
    }
}
#pragma endregion 6x6BottomDrain
