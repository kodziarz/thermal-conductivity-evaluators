#include <gtest/gtest.h>
#include <cmath>
#include "SequentialHeatSimulation.h"
#include "ParallelHeatSimulation.h"
#include "Fenotype.h"
#include <algorithm>

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

#pragma region 6x6BottomDrain
TEST(SequentialHeatSimulationTest, EvaluateGenerationWith6x6GeneratorFullBottomDrainFenotype)
{

    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEta = TEST_REL_ERROR / 1000;
    // const simulation_value_t expectedFitness = -280.672727;
    // const simulation_steps_index_t expectedEquilibriumMoment = 1'331;
    // const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
    //     280,
    //     280,
    //     280,
    //     280,
    //     280,
    //     280,
    //     280,
    //     280.672727,
    //     280.672727,
    //     280.672727,
    //     280.672727,
    //     280,
    //     280,
    //     280.572727,
    //     280.572727,
    //     280.572727,
    //     280.572727,
    //     280,
    //     280,
    //     280.372727,
    //     280.372727,
    //     280.372727,
    //     280.372727,
    //     280,
    //     280,
    //     280.072727,
    //     280.072727,
    //     280.072727,
    //     280.072727,
    //     280,
    //     280,
    //     280,
    //     280,
    //     280,
    //     280,
    //     280,
    // };

    // test

    Fenotype_t fen = Fenotype::createGeneratorFenotype(boardHeight, boardWidth);
    std::vector<cell_type_t> fenotypes(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];
    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    SimulationParams testParams = {
        .simulationSteps = 0,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 100,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 10,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    std::ofstream seqFile("seqTsRegister.txt");
    std::ofstream parFile("parTsRegister.txt");

    // Check if the file was opened successfully
    if (seqFile.is_open() && parFile.is_open())
    {

        for (int i = 0; i < 10; i++)
        {
            testParams.simulationSteps = i;

            SequentialHeatSimulation sequentialSimulator(boardHeight, boardWidth, testParams);
            simulation_value_t *sequentialMinTemperatures = new simulation_value_t[boardHeight * boardWidth];
            simulation_steps_index_t sequentialEquilibriumStep;

            auto sequentialResult = sequentialSimulator.evaluateGeneration(fenotypes, sequentialMinTemperatures, &sequentialEquilibriumStep);

            ParallelHeatSimulation parallelSimulator(boardHeight, boardWidth, testParams);
            simulation_value_t *parallelMinTemperatures = new simulation_value_t[boardHeight * boardWidth];
            simulation_steps_index_t parallelEquilibriumStep;

            auto parallelResult = parallelSimulator.evaluateGeneration(fenotypes, parallelMinTemperatures, &parallelEquilibriumStep);

            // assert

            ASSERT_EQ(sequentialResult.size(), 1u);
            ASSERT_EQ(parallelResult.size(), 1u);

            for (int row = 0; row < boardHeight; row++)
            {
                for (int column = 0; column < boardWidth; column++)
                {
                    seqFile << std::fixed << std::setprecision(6) << sequentialMinTemperatures[row * boardWidth + column] << " ";
                    parFile << std::fixed << std::setprecision(6) << parallelMinTemperatures[row * boardWidth + column] << " ";
                }
                seqFile << std::endl;
                parFile << std::endl;
            }
            seqFile << std::endl;
            parFile << std::endl;

            delete[] sequentialMinTemperatures;
            delete[] parallelMinTemperatures;
        }

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

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
