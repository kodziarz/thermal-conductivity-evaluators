#include <gtest/gtest.h>
#include <cmath>
#include "SequentialHeatSimulation.h"
#include "Fenotype.h"
#include <algorithm>

#define TEST_ERROR 0.0001

#pragma region 6x6BottomDrain
TEST(SequentialHeatSimulationTest, EvaluateGenerationWith6x6GeneratorFullBottomDrainFenotype)
{

    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEps = TEST_ERROR / 1000;
    const simulation_value_t expectedFitness = -280.672727;
    const simulation_steps_index_t expectedEquilibriumMoment = 1'949;
    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280.672727,
        280.672727,
        280.672727,
        280.672727,
        280,
        280,
        280.572727,
        280.572727,
        280.572727,
        280.572727,
        280,
        280,
        280.372727,
        280.372727,
        280.372727,
        280.372727,
        280,
        280,
        280.072727,
        280.072727,
        280.072727,
        280.072727,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
    };

    // test

    Fenotype_t fen = Fenotype::createGeneratorFenotype(boardHeight, boardWidth);
    std::vector<cell_type_t> fenotypes(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];
    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 100,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 10,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .EPS = usedEps};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_NEAR(result[0], expectedFitness, TEST_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_NEAR(minTemperatures[row * boardWidth + column], expectedMinTemperatures[row * boardWidth + column], TEST_ERROR);
        }
    }
    ASSERT_EQ(equilibriumStep, expectedEquilibriumMoment);
    delete[] minTemperatures;
}
#pragma endregion 6x6BottomDrain

#pragma region 6x6LowDrainAl
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithGeneratorFenotypeLowDrainAlpha)
{
    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEps = TEST_ERROR / 1000;
    const simulation_value_t expectedFintess = -281.000000;
    const simulation_steps_index_t expectedEquilibriumMoment = 2'903;
    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        281.000000,
        281.000000,
        281.000000,
        281.000000,
        280,
        280,
        280.900000,
        280.900000,
        280.900000,
        280.900000,
        280,
        280,
        280.700000,
        280.700000,
        280.700000,
        280.700000,
        280,
        280,
        280.400000,
        280.400000,
        280.400000,
        280.400000,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
    };

    // test

    Fenotype_t fen = Fenotype::createGeneratorFenotype(boardHeight, boardWidth);
    std::vector<cell_type_t> fenotypes(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];
    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 10,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 10,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .EPS = usedEps};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_NEAR(result[0], expectedFintess, TEST_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_NEAR(minTemperatures[row * boardWidth + column], expectedMinTemperatures[row * boardWidth + column], TEST_ERROR);
        }
    }
    ASSERT_EQ(equilibriumStep, expectedEquilibriumMoment);
    delete[] minTemperatures;
}
#pragma endregion 6x6LowDrainAl

#pragma region 6x6StartTs
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithGeneratorFenotypeAndStartTemperatures)
{
    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEps = TEST_ERROR / 1000;
    const simulation_value_t expectedFintess = -281.000000;
    const simulation_steps_index_t expectedEquilibriumMoment = 2390;
    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        281.000000,
        281.000000,
        281.000000,
        281.000000,
        280,
        280,
        280.900000,
        280.900000,
        280.900000,
        280.900000,
        280,
        280,
        280.700000,
        280.700000,
        280.700000,
        280.700000,
        280,
        280,
        280.400000,
        280.400000,
        280.400000,
        280.400000,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
    };

    // test

    Fenotype_t fen = Fenotype::createGeneratorFenotype(boardHeight, boardWidth);
    std::vector<cell_type_t> fenotypes(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            if (row == 0 || column == 0 || row == boardWidth - 1 || column == boardHeight - 1)
            {
                startTemperatures[row * boardHeight + column] = 280;
            }
            else
            {
                startTemperatures[row * boardHeight + column] = 280.7;
            }
        }
    }

    // std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280.1);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 10,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 10,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .EPS = usedEps};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_NEAR(result[0], expectedFintess, TEST_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_NEAR(minTemperatures[row * boardWidth + column], expectedMinTemperatures[row * boardWidth + column], TEST_ERROR);
        }
    }
    ASSERT_EQ(equilibriumStep, expectedEquilibriumMoment);
    delete[] minTemperatures;
}
#pragma endregion 6x6StartTs

#pragma region 6x6LeftStrip
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithLeftConductorStripFenotype)
{
    // prepare
    int boardHeight = 8;
    int boardWidth = 8;

    const simulation_value_t usedEps = TEST_ERROR / 1000;
    const simulation_value_t expectedFintess = -281.000000;
    const simulation_steps_index_t expectedEquilibriumMoment = 500'000;
    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        281.000000,
        281.000000,
        281.000000,
        281.000000,
        280,
        280,
        280.900000,
        280.900000,
        280.900000,
        280.900000,
        280,
        280,
        280.700000,
        280.700000,
        280.700000,
        280.700000,
        280,
        280,
        280.400000,
        280.400000,
        280.400000,
        280.400000,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
    };

    // test

    Fenotype_t fen = Fenotype::createLeftConductorStripFenotype(boardHeight, boardWidth, 2);
    std::vector<cell_type_t> fenotypes(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];

    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 100,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 10,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .EPS = usedEps};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            std::cout << std::fixed << std::setprecision(6) << minTemperatures[row * boardWidth + column] << ",\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_NEAR(result[0], expectedFintess, TEST_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_NEAR(minTemperatures[row * boardWidth + column], expectedMinTemperatures[row * boardWidth + column], TEST_ERROR);
        }
    }
    ASSERT_EQ(equilibriumStep, expectedEquilibriumMoment);
    delete[] minTemperatures;
}
#pragma endregion 6x6LeftStrip

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
