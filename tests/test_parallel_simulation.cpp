#include <gtest/gtest.h>
#include <cmath>
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
TEST(ParallelHeatSimulationTest, EvaluateGenerationWith6x6GeneratorFullBottomDrainFenotype)
{

    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEta = TEST_REL_ERROR / 1000;
    const simulation_value_t expectedFitness = -280.672727;
    const simulation_steps_index_t expectedEquilibriumMoment = 1'331;
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
        .ETA = usedEta};

    ParallelHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_IN_REL_ERROR(result[0], expectedFitness, TEST_REL_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(minTemperatures[row * boardWidth + column], expectedMinTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_EQ(equilibriumStep, expectedEquilibriumMoment);
    delete[] minTemperatures;
}
#pragma endregion 6x6BottomDrain

#pragma region 6x6Multiple
TEST(ParallelHeatSimulationTest, EvaluateGenerationWithMultiple6x6GeneratorFullBottomDrainFenotype)
{

    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEta = TEST_REL_ERROR / 1000;
    const simulation_value_t expectedFitness = -280.672727;
    const simulation_steps_index_t expectedEquilibriumMoment = 1'331;
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

    for (int i = 0; i < 9; i++)
    {
        fenotypes.insert(fenotypes.end(), fen, fen + boardHeight * boardWidth);
    }
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
        .ETA = usedEta};

    ParallelHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 10u);
    for (int index = 0; index < 10; index++)
    {
        ASSERT_IN_REL_ERROR(result[index], expectedFitness, TEST_REL_ERROR);
    }

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(minTemperatures[row * boardWidth + column], expectedMinTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_EQ(equilibriumStep, expectedEquilibriumMoment);
    delete[] minTemperatures;
}
#pragma endregion 6x6Multiple

#pragma region 6x6LowDrainAl
TEST(ParallelHeatSimulationTest, EvaluateGenerationWithGeneratorFenotypeLowDrainAlpha)
{
    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEta = TEST_REL_ERROR / 1000;
    const simulation_value_t expectedFitness = -281.000000;
    const simulation_steps_index_t expectedEquilibriumMoment = 1'983;
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
        .ETA = usedEta};

    ParallelHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_IN_REL_ERROR(result[0], expectedFitness, TEST_REL_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(minTemperatures[row * boardWidth + column], expectedMinTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_EQ(equilibriumStep, expectedEquilibriumMoment);
    delete[] minTemperatures;
}
#pragma endregion 6x6LowDrainAl

#pragma region 6x6StartTs
TEST(ParallelHeatSimulationTest, EvaluateGenerationWithGeneratorFenotypeAndStartTemperatures)
{
    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEta = TEST_REL_ERROR / 1000;
    const simulation_value_t expectedFitness = -281.000000;
    const simulation_steps_index_t expectedEquilibriumMoment = 1'470;
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
        .ETA = usedEta};

    ParallelHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_IN_REL_ERROR(result[0], expectedFitness, TEST_REL_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(minTemperatures[row * boardWidth + column], expectedMinTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_EQ(equilibriumStep, expectedEquilibriumMoment);
    delete[] minTemperatures;
}
#pragma endregion 6x6StartTs

#pragma region 16x8LeftStrip
TEST(ParallelHeatSimulationTest, EvaluateGenerationWithLeftConductorStripFenotype)
{
    // prepare
    int boardHeight = 16;
    int boardWidth = 8;

    const simulation_value_t usedEta = TEST_REL_ERROR / 10'000;
    const simulation_value_t expectedFitness = -287.85321037376912;
    const simulation_steps_index_t expectedEquilibriumMoment = 19'543;
    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        281.909756,
        281.929283,
        282.006648,
        284.933180,
        286.882217,
        287.856089,
        280,
        280,
        281.890229,
        281.909740,
        281.987047,
        284.910675,
        286.857382,
        287.829961,
        280,
        280,
        281.851192,
        281.870667,
        281.947848,
        284.865090,
        286.806675,
        287.776412,
        280,
        280,
        281.792680,
        281.812093,
        281.889054,
        284.795162,
        286.727818,
        287.692599,
        280,
        280,
        281.714755,
        281.734068,
        281.810670,
        284.698685,
        286.616834,
        287.573569,
        280,
        280,
        281.617517,
        281.636670,
        281.712706,
        284.572073,
        286.467265,
        287.411273,
        280,
        280,
        281.501125,
        281.520028,
        281.595175,
        284.409635,
        286.268881,
        287.192985,
        280,
        280,
        281.365831,
        281.384339,
        281.458098,
        284.202413,
        286.005637,
        286.898802,
        280,
        280,
        281.212029,
        281.229911,
        281.301501,
        283.936282,
        285.652452,
        286.497783,
        280,
        280,
        281.040345,
        281.057212,
        281.125406,
        283.588762,
        285.170107,
        285.942096,
        280,
        280,
        280.851793,
        280.866941,
        280.929785,
        283.123255,
        284.497116,
        285.158397,
        280,
        280,
        280.648095,
        280.660081,
        280.714343,
        282.477356,
        283.536705,
        284.035980,
        280,
        280,
        280.432410,
        280.437806,
        280.476118,
        281.535122,
        282.136368,
        282.412837,
        280,
        280,
        280.211329,
        280.201578,
        280.113655,
        280.050646,
        280.060807,
        280.066165,
        280,
        280,
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
        .GENERATOR_ALPHA = 1,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    ParallelHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_IN_REL_ERROR(result[0], expectedFitness, TEST_REL_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(minTemperatures[row * boardWidth + column], expectedMinTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_EQ(equilibriumStep, expectedEquilibriumMoment);
    delete[] minTemperatures;
}
#pragma endregion 16x8LeftStrip

#pragma region 8x8_Isotropy
TEST(ParallelHeatSimulationTest, Evaluate4RotatedGenerationsWithLeftConductorStripFenotype)
{
    // prepare
    int boardLength = 8;

    const simulation_value_t usedEta = TEST_REL_ERROR / 100'000;
    const simulation_value_t expectedFitness = -285.17430098531764;
    const simulation_steps_index_t expectedEquilibriumMoment = 14'958;
    const simulation_value_t expectedUpwardsMinTemperatures[boardLength * boardLength] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280.296288,
        280.313184,
        280.381221,
        282.836493,
        284.408762,
        285.174301,
        280,
        280,
        280.279392,
        280.295721,
        280.361825,
        282.719496,
        284.215493,
        284.939840,
        280,
        280,
        280.246167,
        280.261205,
        280.323019,
        282.464174,
        283.793873,
        284.429725,
        280,
        280,
        280.197903,
        280.210511,
        280.264666,
        282.020307,
        283.066101,
        283.555462,
        280,
        280,
        280.137030,
        280.145079,
        280.185465,
        281.286289,
        281.894760,
        282.170561,
        280,
        280,
        280.068110,
        280.067299,
        280.044939,
        280.044623,
        280.056091,
        280.061460,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
    };

    const simulation_value_t *expectedLeftwardsMinTemperatures = Fenotype::rotateLeftSquareFenotypeBy90Deg(boardLength, expectedUpwardsMinTemperatures);
    const simulation_value_t *expectedDownwardsMinTemperatures = Fenotype::rotateLeftSquareFenotypeBy90Deg(boardLength, expectedLeftwardsMinTemperatures);
    const simulation_value_t *expectedRightwardsMinTemperatures = Fenotype::rotateLeftSquareFenotypeBy90Deg(boardLength, expectedDownwardsMinTemperatures);

    const simulation_value_t *expectedMinTemperatures[4] = {
        expectedUpwardsMinTemperatures,
        expectedLeftwardsMinTemperatures,
        expectedDownwardsMinTemperatures,
        expectedRightwardsMinTemperatures,
    };

    // test

    Fenotype_t upFen = Fenotype::createLeftConductorStripFenotype(boardLength, boardLength, 2);
    std::vector<cell_type_t> upwardsFenotypes(upFen, upFen + boardLength * boardLength);

    Fenotype_t leftFen = Fenotype::rotateLeftSquareFenotypeBy90Deg(boardLength, upFen);
    std::vector<cell_type_t> leftwardsFenotypes(leftFen, leftFen + boardLength * boardLength);

    Fenotype_t downFen = Fenotype::rotateLeftSquareFenotypeBy90Deg(boardLength, leftFen);
    std::vector<cell_type_t> downwardsFenotypes(downFen, downFen + boardLength * boardLength);

    Fenotype_t rightFen = Fenotype::rotateLeftSquareFenotypeBy90Deg(boardLength, downFen);
    std::vector<cell_type_t> rightwardsFenotypes(rightFen, rightFen + boardLength * boardLength);

    std::vector<cell_type_t> fenotypes[4] = {
        upwardsFenotypes,
        leftwardsFenotypes,
        downwardsFenotypes,
        rightwardsFenotypes,
    };

    delete[] upFen;
    delete[] leftFen;
    delete[] downFen;
    delete[] rightFen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardLength * boardLength];

    std::fill_n(startTemperatures, boardLength * boardLength, (simulation_value_t)280);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 100,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 1,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    ParallelHeatSimulation simulator(boardLength, boardLength, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardLength * boardLength];
    simulation_steps_index_t equilibriumStep;

    for (int index = 0; index < 4; index++)
    {

        auto result = simulator.evaluateGeneration(fenotypes[index], minTemperatures, &equilibriumStep);

        // assert

        ASSERT_EQ(result.size(), 1u);
        ASSERT_IN_REL_ERROR(result[0], expectedFitness, TEST_REL_ERROR);

        for (int row = 0; row < boardLength; row++)
        {
            for (int column = 0; column < boardLength; column++)
            {
                ASSERT_IN_REL_ERROR(minTemperatures[row * boardLength + column], expectedMinTemperatures[index][row * boardLength + column], TEST_REL_ERROR);
            }
        }
        ASSERT_EQ(equilibriumStep, expectedEquilibriumMoment);
    }
    delete[] minTemperatures;
}
#pragma endregion 8x8_Isotropy

#pragma region 32x32LeftStrip
TEST(ParallelHeatSimulationTest, EvaluateGenerationWithBigLeftConductorStripFenotype)
{
    // prepare
    int boardHeight = 32;
    int boardWidth = 32;

    const simulation_value_t usedEta = TEST_REL_ERROR / 100'000;
    const simulation_value_t expectedFitness = -535.19777089204274;
    const simulation_steps_index_t expectedEquilibriumMoment = 534'733;
    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        316.990313,
        317.081414,
        317.442448,
        334.766361,
        351.205802,
        366.784597,
        381.526456,
        395.454874,
        408.593011,
        420.963601,
        432.588848,
        443.490340,
        453.688961,
        463.204818,
        472.057173,
        480.264380,
        487.843829,
        494.811905,
        501.183939,
        506.974183,
        512.195771,
        516.860704,
        520.979822,
        524.562795,
        527.618103,
        530.153032,
        532.173665,
        533.684872,
        534.690310,
        535.192419,
        280,
        280,
        316.899211,
        316.990194,
        317.350759,
        334.650859,
        351.066487,
        366.621577,
        381.339955,
        395.245218,
        408.360631,
        420.709025,
        432.312695,
        443.193309,
        453.371830,
        462.868434,
        471.702443,
        479.892264,
        487.455337,
        494.408085,
        500.765876,
        506.542989,
        511.752583,
        516.406679,
        520.516133,
        524.090626,
        527.138652,
        529.667503,
        531.683266,
        533.190817,
        534.193817,
        534.694706,
        280,
        280,
        316.717127,
        316.807872,
        317.167496,
        334.419859,
        350.787745,
        366.295317,
        380.966623,
        394.825475,
        407.895345,
        420.199256,
        431.759685,
        442.598469,
        452.736722,
        462.194758,
        470.992020,
        479.147025,
        486.677304,
        493.599363,
        499.928637,
        505.679464,
        510.865050,
        515.497455,
        519.587567,
        523.145094,
        526.178546,
        528.695233,
        530.701253,
        532.201491,
        533.199612,
        533.698061,
        280,
        280,
        316.444298,
        316.534684,
        316.892892,
        334.073362,
        350.369353,
        365.805371,
        380.405799,
        394.194778,
        407.196089,
        419.433050,
        430.928410,
        441.704258,
        451.781938,
        461.181967,
        469.923975,
        478.026637,
        485.507625,
        492.383564,
        498.669991,
        504.381328,
        509.530853,
        514.130682,
        518.191751,
        521.723802,
        524.735375,
        527.233802,
        529.225198,
        530.714457,
        531.705255,
        532.200041,
        280,
        280,
        316.081083,
        316.170989,
        316.527297,
        333.611373,
        349.810969,
        365.151059,
        379.656480,
        393.351812,
        406.261257,
        418.408525,
        429.816736,
        440.508313,
        450.504907,
        459.827310,
        468.495394,
        476.528048,
        483.943127,
        490.757414,
        496.986578,
        502.645152,
        507.746506,
        512.302828,
        516.325114,
        519.823152,
        522.805519,
        525.279573,
        527.251450,
        528.726059,
        529.707085,
        530.196983,
        280,
        280,
        315.627963,
        315.717263,
        316.071180,
        333.033891,
        349.112129,
        364.331461,
        378.717304,
        392.294796,
        405.088671,
        417.123139,
        428.421781,
        439.007447,
        448.902169,
        458.127083,
        466.702362,
        474.647157,
        481.979553,
        488.716522,
        494.873897,
        500.466344,
        505.507342,
        510.009167,
        513.982884,
        517.438336,
        520.384142,
        522.827691,
        524.775140,
        526.231416,
        527.200214,
        527.683998,
        280,
        280,
        315.085544,
        315.174111,
        315.525131,
        332.340907,
        348.272231,
        363.345396,
        377.586531,
        391.021459,
        403.675561,
        415.573657,
        426.739889,
        437.197619,
        446.969341,
        456.076599,
        464.539928,
        472.378789,
        479.611533,
        486.255358,
        492.326284,
        497.839130,
        502.807499,
        507.243767,
        511.159076,
        514.563328,
        517.465186,
        519.872073,
        521.790170,
        523.224420,
        524.178529,
        524.654966,
        280,
        280,
        314.454557,
        314.542260,
        314.889867,
        331.532402,
        347.290526,
        362.191405,
        376.262019,
        389.529007,
        402.018526,
        413.756117,
        424.766584,
        435.073893,
        444.701075,
        453.670153,
        462.002074,
        469.716659,
        476.832559,
        483.367226,
        489.336888,
        494.756535,
        499.639904,
        503.999478,
        507.846479,
        511.190871,
        514.041363,
        516.405409,
        518.289212,
        519.697732,
        520.634683,
        521.102539,
        280,
        280,
        313.735867,
        313.822570,
        314.166228,
        330.608335,
        346.166100,
        360.867722,
        374.741183,
        387.814085,
        400.113488,
        411.665775,
        422.496521,
        432.630383,
        442.091011,
        450.900968,
        459.081668,
        466.653330,
        473.634940,
        480.044227,
        485.897643,
        491.210357,
        495.996249,
        500.267909,
        504.036641,
        507.312469,
        510.104143,
        512.419146,
        514.263701,
        515.642776,
        516.560094,
        517.018134,
        280,
        280,
        312.930475,
        313.016037,
        313.355188,
        329.568636,
        344.897851,
        359.372240,
        373.020957,
        385.872718,
        397.955632,
        409.297049,
        419.923421,
        429.860196,
        439.131714,
        447.761142,
        455.770409,
        463.180169,
        470.009767,
        476.277225,
        481.999231,
        487.191138,
        491.866967,
        496.039412,
        499.719855,
        502.918373,
        505.643749,
        507.903488,
        509.703825,
        511.049736,
        511.944944,
        512.391929,
        280,
        280,
        312.039521,
        312.123795,
        312.457860,
        328.413192,
        343.484461,
        357.702472,
        371.097735,
        383.700254,
        395.539338,
        406.643438,
        417.040000,
        426.755350,
        435.814601,
        444.241576,
        452.058764,
        459.287281,
        465.946850,
        472.055797,
        477.631046,
        482.688130,
        487.241204,
        491.303058,
        494.885140,
        497.997564,
        500.649141,
        502.847384,
        504.598530,
        505.907554,
        506.778174,
        507.212866,
        280,
        280,
        311.064292,
        311.147124,
        311.475494,
        327.141836,
        341.924359,
        355.855492,
        368.967305,
        381.291279,
        392.858090,
        403.697434,
        413.837866,
        423.306689,
        432.129852,
        440.331893,
        447.935892,
        454.963448,
        461.434670,
        467.368186,
        472.781148,
        477.689260,
        482.106792,
        486.046614,
        489.520219,
        492.537747,
        495.108011,
        497.238522,
        498.935507,
        500.203924,
        501.047482,
        501.468645,
        280,
        280,
        310.006232,
        310.087459,
        310.409494,
        325.754320,
        340.215678,
        353.827869,
        366.624760,
        378.639518,
        389.904371,
        400.450408,
        410.307417,
        419.503767,
        428.066313,
        436.020346,
        443.389563,
        450.196052,
        456.460308,
        462.201243,
        467.436222,
        472.181092,
        476.450218,
        480.256520,
        483.611511,
        486.525330,
        489.006774,
        491.063330,
        492.701194,
        493.925300,
        494.739330,
        495.145734,
        280,
        280,
        308.866944,
        308.946395,
        309.261415,
        324.250294,
        338.356195,
        351.615583,
        364.064392,
        375.737714,
        386.669526,
        396.892473,
        406.437700,
        415.334725,
        423.611369,
        431.293706,
        438.406055,
        444.970992,
        451.009370,
        456.540366,
        461.581520,
        466.148788,
        470.256590,
        473.917862,
        477.144104,
        479.945419,
        482.330561,
        484.306965,
        485.880778,
        487.056889,
        487.838945,
        488.229368,
        280,
        280,
        307.648205,
        307.725696,
        308.032979,
        322.629266,
        336.343252,
        349.213910,
        361.279554,
        372.577468,
        383.143600,
        393.012321,
        402.216251,
        410.786140,
        418.750812,
        426.137138,
        432.970052,
        439.272585,
        445.065916,
        450.369437,
        455.200815,
        459.576063,
        463.509609,
        467.014357,
        470.101746,
        472.781808,
        475.063215,
        476.953321,
        478.458198,
        479.582667,
        480.330325,
        480.703559,
        280,
        280,
        306.351975,
        306.427310,
        306.726084,
        320.890561,
        334.173662,
        346.617285,
        358.262486,
        369.149050,
        379.315138,
        388.797020,
        397.628906,
        405.842843,
        413.468677,
        420.534064,
        427.064517,
        433.083471,
        438.612370,
        443.670751,
        448.276343,
        452.445149,
        456.191539,
        459.528324,
        462.466833,
        465.016973,
        467.187294,
        468.985030,
        470.416149,
        471.485383,
        472.196258,
        472.551113,
        280,
        280,
        304.980411,
        305.053378,
        305.342810,
        319.033250,
        331.843577,
        343.819112,
        355.004092,
        365.441154,
        375.170931,
        384.231771,
        392.659573,
        400.487716,
        407.747060,
        414.466001,
        420.670562,
        426.384500,
        431.629431,
        436.424951,
        440.788755,
        444.736754,
        448.283179,
        451.440678,
        454.220399,
        456.632073,
        458.684072,
        460.383474,
        461.736105,
        462.746579,
        463.418332,
        463.753641,
        280,
        280,
        303.535880,
        303.606250,
        303.885438,
        317.056068,
        329.348310,
        340.811524,
        351.493651,
        361.440584,
        370.695707,
        379.299611,
        387.289956,
        394.701451,
        401.565914,
        407.912392,
        413.767307,
        419.154618,
        424.095990,
        428.610956,
        432.717066,
        436.430030,
        439.763846,
        442.730911,
        445.342119,
        447.606953,
        449.533557,
        451.128800,
        452.398329,
        453.346609,
        453.976964,
        454.291593,
        280,
        280,
        302.020980,
        302.088501,
        302.356460,
        314.957290,
        326.682093,
        337.585050,
        347.718438,
        357.131864,
        365.871747,
        373.981060,
        381.499244,
        388.462275,
        394.902816,
        400.850415,
        406.331727,
        411.370750,
        415.989036,
        420.205900,
        424.038609,
        427.502545,
        430.611359,
        433.377095,
        435.810310,
        437.920163,
        439.714506,
        441.199945,
        442.381905,
        443.264671,
        443.851427,
        444.144280,
        280,
        280,
        300.438559,
        300.502954,
        300.758598,
        312.734556,
        323.837744,
        334.128169,
        343.663219,
        352.496720,
        360.678399,
        368.253683,
        375.263733,
        381.745645,
        387.732721,
        393.254785,
        398.338505,
        403.007690,
        407.283577,
        411.185077,
        414.729005,
        417.930266,
        420.802034,
        423.355891,
        425.601953,
        427.548978,
        429.204452,
        430.574665,
        431.664773,
        432.478842,
        433.019891,
        433.289919,
        280,
        280,
        298.791743,
        298.852702,
        299.094823,
        310.384606,
        320.806175,
        330.426688,
        339.309578,
        347.513431,
        355.091484,
        362.091580,
        368.556408,
        374.523901,
        380.027690,
        385.097559,
        389.759878,
        394.037994,
        397.952572,
        401.521899,
        404.762141,
        407.687558,
        410.310701,
        412.642563,
        414.692717,
        416.469428,
        417.979745,
        419.229579,
        420.223767,
        420.966121,
        421.459468,
        421.705677,
        280,
        280,
        297.083967,
        297.141139,
        297.368373,
        307.902882,
        317.575680,
        326.462853,
        334.635000,
        342.155972,
        349.082561,
        355.464782,
        361.346459,
        366.765908,
        371.756629,
        376.347935,
        380.565511,
        384.431893,
        387.966880,
        391.187873,
        394.110167,
        396.747194,
        399.110720,
        401.211017,
        403.057001,
        404.656351,
        406.015600,
        407.140219,
        408.034677,
        408.702490,
        409.146266,
        409.367725,
        280,
        280,
        295.319020,
        295.371993,
        295.582781,
        305.282881,
        314.130825,
        322.214063,
        329.611620,
        336.392923,
        342.618035,
        348.338564,
        353.598775,
        358.436683,
        362.885026,
        366.972088,
        370.722387,
        374.157242,
        377.295238,
        380.152603,
        382.743523,
        385.080394,
        387.174033,
        389.033849,
        390.667987,
        392.083443,
        393.286157,
        394.281093,
        395.072303,
        395.662970,
        396.055453,
        396.251307,
        280,
        280,
        293.501099,
        293.549373,
        293.741888,
        302.515046,
        310.450691,
        317.650971,
        324.204513,
        330.186088,
        335.658121,
        340.672694,
        345.273426,
        349.497060,
        353.374744,
        356.933046,
        360.194751,
        363.179497,
        365.904274,
        368.383829,
        370.630981,
        372.656881,
        374.471225,
        376.082418,
        377.497717,
        378.723338,
        379.764553,
        380.625757,
        381.310535,
        381.821700,
        382.161333,
        382.330808,
        280,
        280,
        291.634905,
        291.677807,
        291.849850,
        299.584734,
        306.505933,
        312.734632,
        318.369393,
        323.488814,
        328.155689,
        332.420691,
        336.325203,
        339.903417,
        343.183878,
        346.190636,
        348.944112,
        351.461759,
        353.758577,
        355.847503,
        357.739735,
        359.444972,
        360.971616,
        362.326932,
        363.517175,
        364.547693,
        365.423012,
        366.146903,
        366.722436,
        367.152017,
        367.437427,
        367.579841,
        280,
        280,
        289.725809,
        289.762262,
        289.911067,
        296.468115,
        302.253684,
        307.412244,
        312.049628,
        316.244103,
        320.055151,
        323.529200,
        326.703303,
        329.607553,
        332.266744,
        334.701539,
        336.929333,
        338.964886,
        340.820805,
        342.507909,
        344.035524,
        345.411696,
        346.643377,
        347.736562,
        348.696401,
        349.527290,
        350.232944,
        350.816455,
        351.280333,
        351.626552,
        351.856565,
        351.971334,
        280,
        280,
        287.780259,
        287.808024,
        287.929876,
        293.122982,
        297.628453,
        301.611044,
        305.172783,
        308.382835,
        311.291626,
        313.937672,
        316.351274,
        318.556771,
        320.574026,
        322.419466,
        324.106821,
        325.647676,
        327.051875,
        328.327836,
        329.482785,
        330.522943,
        331.453668,
        332.279572,
        333.004612,
        333.632159,
        334.165057,
        334.605674,
        334.955929,
        335.217329,
        335.390986,
        335.477633,
        280,
        280,
        285.806946,
        285.820015,
        285.909126,
        289.465488,
        292.526108,
        295.230702,
        297.647634,
        299.822837,
        301.790858,
        303.578601,
        305.207365,
        306.694245,
        308.053142,
        309.295497,
        310.430830,
        311.467140,
        312.411205,
        313.268797,
        314.044861,
        314.743646,
        315.368804,
        315.923473,
        316.410343,
        316.831701,
        317.189480,
        317.485281,
        317.720408,
        317.895878,
        318.012446,
        318.070607,
        280,
        280,
        283.820563,
        283.800073,
        283.832142,
        285.303739,
        286.779794,
        288.138027,
        289.364222,
        290.470028,
        291.470377,
        292.378517,
        293.205349,
        293.959714,
        294.648812,
        295.278561,
        295.853873,
        296.378866,
        296.857021,
        297.291300,
        297.684233,
        298.037993,
        298.354445,
        298.635191,
        298.881601,
        299.094843,
        299.275898,
        299.425583,
        299.544561,
        299.633349,
        299.692332,
        299.721761,
        280,
        280,
        281.854670,
        281.743447,
        280.903046,
        280.137534,
        280.151303,
        280.177393,
        280.201202,
        280.222679,
        280.242108,
        280.259746,
        280.275805,
        280.290455,
        280.303838,
        280.316068,
        280.327241,
        280.337436,
        280.346722,
        280.355155,
        280.362786,
        280.369656,
        280.375801,
        280.381253,
        280.386038,
        280.390179,
        280.393695,
        280.396602,
        280.398912,
        280.400636,
        280.401782,
        280.402353,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
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
        .GENERATOR_ALPHA = 1,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    ParallelHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_IN_REL_ERROR(result[0], expectedFitness, TEST_REL_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(minTemperatures[row * boardWidth + column], expectedMinTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_IN_REL_ERROR(equilibriumStep, expectedEquilibriumMoment, 0.00001);
    delete[] minTemperatures;
}
#pragma endregion 32x32LeftStrip

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
