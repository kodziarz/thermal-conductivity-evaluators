#include <gtest/gtest.h>
#include <cmath>
#include "SequentialHeatSimulation.h"
#include "Fenotype.h"
#include <algorithm>

#define TEST_REL_ERROR 0.00001
#define TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR 0.00001

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
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 100,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 10,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_IN_REL_ERROR(expectedFitness, result[0], TEST_REL_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(expectedMinTemperatures[row * boardWidth + column], minTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_IN_REL_ERROR(expectedEquilibriumMoment, equilibriumStep, TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR);
    delete[] minTemperatures;
}
#pragma endregion 6x6BottomDrain

#pragma region 6x6Multiple
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithMultiple6x6GeneratorFullBottomDrainFenotype)
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
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 100,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 10,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 10u);
    for (int index = 0; index < 10; index++)
    {
        ASSERT_IN_REL_ERROR(expectedFitness, result[index], TEST_REL_ERROR);
    }

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(expectedMinTemperatures[row * boardWidth + column], minTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_IN_REL_ERROR(expectedEquilibriumMoment, equilibriumStep, TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR);

    delete[] minTemperatures;
}
#pragma endregion 6x6Multiple

#pragma region 6x6LowDrainAl
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithGeneratorFenotypeLowDrainAlpha)
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
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 10,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 10,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_IN_REL_ERROR(expectedFitness, result[0], TEST_REL_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(expectedMinTemperatures[row * boardWidth + column], minTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_IN_REL_ERROR(expectedEquilibriumMoment, equilibriumStep, TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR);

    delete[] minTemperatures;
}
#pragma endregion 6x6LowDrainAl

#pragma region 6x6StartTs
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithGeneratorFenotypeAndStartTemperatures)
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
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 10,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 10,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_IN_REL_ERROR(expectedFitness, result[0], TEST_REL_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(expectedMinTemperatures[row * boardWidth + column], minTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_IN_REL_ERROR(expectedEquilibriumMoment, equilibriumStep, TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR);
    delete[] minTemperatures;
}
#pragma endregion 6x6StartTs

#pragma region 16x8LeftStrip
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithLeftConductorStripFenotype)
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
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 100,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 1,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_IN_REL_ERROR(expectedFitness, result[0], TEST_REL_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(expectedMinTemperatures[row * boardWidth + column], minTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_IN_REL_ERROR(expectedEquilibriumMoment, equilibriumStep, TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR);
    delete[] minTemperatures;
}
#pragma endregion 16x8LeftStrip

#pragma region 8x8_Isotropy
TEST(SequentialHeatSimulationTest, Evaluate4RotatedGenerationsWithLeftConductorStripFenotype)
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
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 100,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 1,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    SequentialHeatSimulation simulator(boardLength, boardLength, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardLength * boardLength];
    simulation_steps_index_t equilibriumStep;

    for (int index = 0; index < 4; index++)
    {

        auto result = simulator.evaluateGeneration(fenotypes[index], minTemperatures, &equilibriumStep);

        // assert

        ASSERT_EQ(result.size(), 1u);
        ASSERT_IN_REL_ERROR(expectedFitness, result[0], TEST_REL_ERROR);

        for (int row = 0; row < boardLength; row++)
        {
            for (int column = 0; column < boardLength; column++)
            {
                ASSERT_IN_REL_ERROR(expectedMinTemperatures[index][row * boardLength + column], minTemperatures[row * boardLength + column], TEST_REL_ERROR);
            }
        }
        ASSERT_IN_REL_ERROR(expectedEquilibriumMoment, equilibriumStep, TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR);
    }
    delete[] minTemperatures;
}
#pragma endregion 8x8_Isotropy

#pragma region 32x32LeftStrip
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithBigLeftConductorStripFenotype)
{
    // prepare
    int boardHeight = 32;
    int boardWidth = 32;

    const simulation_value_t usedEta = TEST_REL_ERROR / 100'000;
    const simulation_value_t expectedFitness = -535.22383168287854;
    const simulation_steps_index_t expectedEquilibriumMoment = 669'885;
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
        316.9934670787159,
        317.0845770903038,
        317.4456440440507,
        334.7712233296016,
        351.2123174373631,
        366.7927452926414,
        381.5362151777811,
        395.4662141215931,
        408.6059008203153,
        420.9780032242480,
        432.6047217507854,
        443.5076388881478,
        453.7076357472531,
        463.2248159120056,
        472.0784367407891,
        480.2868480920986,
        487.8674382899695,
        494.8365870142345,
        501.2096246978929,
        507.0007979385419,
        512.2232403816897,
        516.8889485085188,
        521.0087617558902,
        524.5923464094078,
        527.6481827376886,
        530.1835548752965,
        532.2045430094621,
        533.7160174805124,
        534.7216344654748,
        535.2238329772831,
        280,
        280,
        316.9023570671281,
        316.9933482902496,
        317.3539459227288,
        334.6557085074193,
        351.0729836898748,
        366.6297032628273,
        381.3496861191561,
        395.2565263667491,
        408.3734851151711,
        420.7233871017241,
        432.3285231400455,
        443.2105591664903,
        453.3904524417100,
        462.8883752480788,
        471.7236462183720,
        479.9146692456604,
        487.4788797636984,
        494.4326980549642,
        500.7914891410438,
        506.5695287361851,
        511.7799746981414,
        516.4348433881091,
        520.5449903499115,
        524.1200947348152,
        527.1686469285319,
        529.6979388789094,
        531.7140566727478,
        533.2218749667705,
        534.2250529387995,
        534.7260314892619,
        280,
        280,
        316.7202558324188,
        316.8110089089153,
        317.1706656569202,
        334.4246810875003,
        350.7942055519180,
        366.3033979496842,
        380.9762996693142,
        394.8367201111425,
        407.9081261719622,
        420.2135369275120,
        431.7754245412675,
        442.6156221961431,
        452.7552396051217,
        462.2145864203318,
        471.0131036390642,
        479.1693029085957,
        486.7007134643228,
        493.6238363010032,
        499.9541050752744,
        505.7058531671554,
        510.8922862867234,
        515.5254599959977,
        519.6162615209953,
        523.1743952515802,
        526.2083713628848,
        528.7254970392320,
        530.7318698360195,
        532.2323727751929,
        533.2306708338613,
        533.7292085518736,
        280,
        280,
        316.4474015212134,
        316.5377959463349,
        316.8960357347906,
        334.0781446337721,
        350.3757594806413,
        365.8133833147244,
        380.4153944973218,
        394.2059282366105,
        407.2087625340895,
        419.4472098951726,
        430.9440159014547,
        441.7212654717782,
        451.8002973624064,
        461.2016271891665,
        469.9448790090614,
        478.0487252854587,
        485.5308348841170,
        492.4078286095744,
        498.6952416920351,
        504.4074925705809,
        509.5578572857413,
        514.1584487882960,
        518.2202004866500,
        521.7528533877960,
        524.7649462323656,
        527.2638080792844,
        529.2555528570759,
        530.7450754642905,
        531.7360490697495,
        532.2309233326680,
        280,
        280,
        316.0841527848866,
        316.1740663157062,
        316.5304062307010,
        333.6161022321843,
        349.8173044221791,
        365.1589813312979,
        379.6659667686852,
        393.3628358039545,
        406.2737858326790,
        418.4225242177060,
        429.8321636976857,
        440.5251264271940,
        450.5230571836631,
        459.8467459649707,
        468.5160599226606,
        476.5498843401834,
        483.9660721772353,
        490.7814015612657,
        497.0115405128467,
        502.6710181375336,
        507.7732014975071,
        512.3302773849277,
        516.3532382496539,
        519.8518715807587,
        522.8347520996672,
        525.3092361886347,
        527.2814580488797,
        528.7563271553142,
        529.7375266483485,
        530.2275123765514,
        280,
        280,
        315.6309905177401,
        315.7202985588249,
        316.0742464323778,
        333.0385536421136,
        349.1183746446213,
        364.3392708196499,
        378.7266554422142,
        392.3056623779096,
        405.1010207750321,
        417.1369374453533,
        428.4369882444732,
        439.0240193557344,
        448.9200589801850,
        458.1462395644966,
        466.7227303765311,
        474.6686799755018,
        482.0021679234983,
        488.7401649455295,
        494.8985006606757,
        500.4918379693416,
        505.5336531819678,
        510.0362210043867,
        514.0106035464122,
        517.4666425860790,
        520.4129543970801,
        522.8569265268781,
        524.8047159946641,
        526.2612484599085,
        527.2302179919493,
        527.7140871488082,
        280,
        280,
        315.0885202095091,
        315.1770951668844,
        315.5281470387508,
        332.3454912592993,
        348.2783696945712,
        363.3530718605135,
        377.5957218026593,
        391.0321374905030,
        403.6876974442527,
        415.5872165442681,
        426.7548324792047,
        437.2139037711706,
        446.9869198169375,
        456.0954229364038,
        464.5599420435695,
        472.3999372619025,
        479.6337545958496,
        486.2785896368015,
        492.3504592151083,
        497.8641798973308,
        502.8333522567781,
        507.2703499043808,
        511.1863123456620,
        514.5911408201977,
        517.4934963758620,
        519.9007995273037,
        521.8192309431605,
        523.2537326978768,
        524.2080097109025,
        524.6845310780943,
        280,
        280,
        314.4574749439027,
        314.5451855370270,
        314.8928229538369,
        331.5368946617897,
        347.2965410138792,
        362.1989251252210,
        376.2710224174540,
        389.5394683372481,
        402.0304149672738,
        413.7693988083282,
        424.7812213569920,
        435.0898434328913,
        444.7182935800757,
        453.6885903207159,
        462.0216775995450,
        469.7373724327934,
        476.8543235613192,
        483.3899797908416,
        489.3605666657485,
        494.7810701482182,
        499.6652260435748,
        504.0255140108386,
        507.8731551117994,
        511.2181119733206,
        514.0690907589992,
        516.4335442634755,
        518.3176755529679,
        519.7264416777061,
        520.6635570758599,
        521.1314963747427,
        280,
        280,
        313.7387190851722,
        313.8254296048050,
        314.1691167098792,
        330.6127234201720,
        346.1719745739634,
        360.8750652090847,
        374.7499744047349,
        387.8242984738091,
        400.1250952793325,
        411.6787423648453,
        422.5108107076289,
        432.6459550234120,
        442.1078207498434,
        450.9189671669433,
        459.1008056012053,
        466.6735513085111,
        473.6561874259106,
        480.0664392996204,
        485.9207575089488,
        491.2343079863422,
        496.0209677586006,
        500.2933249837415,
        504.0626821175185,
        507.3390612024279,
        510.1312104234718,
        512.4466112147638,
        514.2914853276717,
        515.6708013842789,
        516.5882805402529,
        517.0464009704421,
        280,
        280,
        312.9332527068090,
        313.0188222041535,
        313.3580025618835,
        329.5729077350842,
        344.9035686527459,
        359.3793867324669,
        373.0295115186390,
        385.8826558739684,
        397.9669253114679,
        409.3096646641579,
        419.9373240853392,
        429.8753452033698,
        439.1480672290277,
        447.7786519961033,
        455.7890263299261,
        463.1998397742393,
        470.0304355342961,
        476.2988324729034,
        482.0217160842075,
        487.2144365297241,
        491.8910120208671,
        496.0641360481468,
        499.7451871722472,
        502.9442402955428,
        505.6700785178383,
        507.9302048445691,
        509.7308531588088,
        511.0769979916179,
        511.9723627305633,
        512.4194259964634,
        280,
        280,
        312.0422168311011,
        312.1264982201933,
        312.4605913060568,
        328.4173363055637,
        343.4900055694978,
        357.7094015494408,
        371.1060290634333,
        383.7098881920050,
        395.5502854284791,
        406.6556668950457,
        417.0534757662664,
        426.7700344757852,
        435.8304509668796,
        444.2585472586015,
        452.0768079482605,
        459.3063459243281,
        465.9668824642355,
        472.0767389736008,
        477.6528378253766,
        482.7107100276030,
        487.2645077471199,
        491.3270200158545,
        494.9096902279174,
        498.0226342898000,
        500.6746585079114,
        502.8732764870077,
        504.6247244715187,
        505.9339746929625,
        506.8047463940612,
        507.2395142885175,
        280,
        280,
        311.0668995663015,
        311.1497386169641,
        311.4781358868433,
        327.1458406116442,
        341.9297157702695,
        355.8621848323935,
        368.9753149936959,
        381.3005824021865,
        392.8686613154635,
        403.7092417213454,
        413.8508776089621,
        423.3208659667102,
        432.1451549041891,
        440.3482781232480,
        447.9533122802833,
        454.9818535106812,
        461.4540094248213,
        467.3884031319920,
        472.8021862162144,
        477.7110580083144,
        482.1292889242782,
        486.0697460403572,
        489.5439194338911,
        492.5619481279514,
        495.1326447371415,
        497.2635181241736,
        498.9607935474378,
        500.2294299147943,
        501.0731338643434,
        501.4943704751701,
        280,
        280,
        310.0087432508391,
        310.0899774431088,
        310.4120378717416,
        325.7581744839285,
        340.2208320675709,
        353.8343070161965,
        366.6324636768173,
        378.6484651076292,
        389.9145357098903,
        400.4617610659768,
        410.3199269815925,
        419.5173968779709,
        428.0810245600039,
        436.0360980500032,
        443.4063095390289,
        450.2137464133963,
        456.4788985924807,
        462.2206779134356,
        467.4564458992791,
        472.2020468652815,
        476.4718439014442,
        480.2787557875284,
        483.6342933394612,
        486.5485940510964,
        489.0304541886528,
        491.0873577252306,
        492.7255016794021,
        493.9498175545741,
        494.7639886734897,
        495.1704632727909,
        280,
        280,
        308.8693527431072,
        308.9488099450635,
        309.2638548904237,
        324.2539873847857,
        338.3611309999174,
        351.6217474880328,
        364.0717675897948,
        375.7462786416700,
        386.6792553505393,
        396.9033398511451,
        406.4496723735267,
        415.3477700036432,
        423.6254484079375,
        431.3087799778173,
        438.4220814125180,
        444.9879240114797,
        451.0271606183736,
        456.5589640300951,
        461.6008726022888,
        466.1688396521926,
        470.2772840288009,
        473.9391398689738,
        477.1659040854524,
        479.9676805484431,
        482.3532202412657,
        484.3299569088172,
        485.9040378904890,
        487.0803499507332,
        487.8625400023737,
        488.2530306698361,
        280,
        280,
        307.6505050334192,
        307.7280019515681,
        308.0353091554390,
        322.6327891649013,
        336.3479570593088,
        349.2197843462508,
        361.2865805527065,
        372.5856265187641,
        383.1528671994991,
        393.0226706146022,
        402.2276526577922,
        410.7985623552041,
        418.7642190903519,
        426.1514920408957,
        432.9853121218314,
        439.2887076017162,
        445.0828558395257,
        450.3871449863863,
        455.2192408276928,
        459.5951551125036,
        463.5293126926970,
        467.0346155742180,
        470.1225025850505,
        472.8030038160808,
        475.0847893192729,
        476.9752117784066,
        478.4803430231267,
        479.6050043556192,
        480.3527907155590,
        480.7260887344667,
        280,
        280,
        306.3541604055825,
        306.4295007382667,
        306.7282972064849,
        320.8939030600934,
        334.1781237261940,
        346.6228522849837,
        358.2691437560475,
        369.1567796812282,
        379.3239163141382,
        388.8068227500200,
        397.6397052879020,
        405.8546076690955,
        413.4813735574367,
        420.5476569736597,
        427.0789674322811,
        433.0987384341133,
        438.6284101517120,
        443.6875192483171,
        448.2937906096966,
        452.4632272775359,
        456.2101960553699,
        459.5475071502549,
        462.4864868645553,
        465.0370428116611,
        467.2077214414489,
        469.0057578625325,
        470.4371180681152,
        471.5065337331808,
        472.2175297698994,
        472.5724448181281,
        280,
        280,
        304.9824754450617,
        305.0554476411995,
        305.3449010489499,
        319.0364021428029,
        331.8477825004185,
        343.8243573114711,
        355.0103625053000,
        365.4484321360101,
        375.1791956258530,
        384.2409987834848,
        392.6697380747593,
        400.4987894759055,
        407.7590104967057,
        414.4787948640916,
        420.6841621996027,
        426.3988685508292,
        431.6445270849771,
        436.4407312455587,
        440.8051750853260,
        444.7537673326718,
        448.3007371010958,
        451.4587301069808,
        454.2388949113586,
        456.6509591246638,
        458.7032957724336,
        460.4029801622634,
        461.7558376537251,
        462.7664827391940,
        463.4383498128338,
        463.7737159501313,
        280,
        280,
        303.5378182884031,
        303.6081927693564,
        303.8874008889670,
        317.0590219617593,
        329.3522468212344,
        340.8164319552104,
        351.4995168176998,
        361.4473907317059,
        370.7034352698261,
        379.3082386833545,
        387.2994587517922,
        394.7118016631277,
        401.5770840894555,
        407.9243497864647,
        413.7800179512751,
        419.1680464847020,
        424.1100983918938,
        428.6257035636999,
        432.7324111534621,
        436.4459298668147,
        439.7802549094459,
        442.7477812653104,
        445.3594035493385,
        447.6246030033059,
        449.5515223614623,
        451.1470293604667,
        452.4167696454320,
        453.3652097571403,
        453.9956707922147,
        454.3103532195362,
        280,
        280,
        302.0227866507914,
        302.0903122780633,
        302.3582896965383,
        314.9600379940426,
        326.6857508675775,
        337.5896068704650,
        347.7238820786112,
        357.1381787033162,
        365.8789160384382,
        373.9890619283622,
        381.5080565859747,
        388.4718743354052,
        394.9131744115899,
        400.8615022411029,
        406.3435133343971,
        411.3832010448763,
        416.0021164342626,
        420.2195734639700,
        424.0528360980928,
        427.5172860717644,
        430.6265714046481,
        433.3927364955617,
        435.8263350174645,
        437.9365269778446,
        439.7311613097382,
        441.2168452728132,
        442.3990018105001,
        443.2819158518246,
        443.8687703794526,
        444.1616729163667,
        280,
        280,
        300.4402293859079,
        300.5046288177122,
        300.7602898396059,
        312.7370894503049,
        323.8411117845960,
        334.1323625804890,
        343.6682259229924,
        352.5025259645378,
        360.6849882522878,
        368.2610364057287,
        375.2718313283865,
        381.7544646809760,
        387.7422369804435,
        393.2649714320261,
        398.3493321004005,
        403.0191279262098,
        407.2955928363766,
        411.1976377598912,
        414.7420737032430,
        417.9438069175871,
        420.8160081419056,
        423.3702582949089,
        425.6166730471986,
        427.5640085809548,
        429.2197506269180,
        430.5901886106330,
        431.6804764720157,
        432.4946814602910,
        433.0358219574896,
        433.3058951501968,
        280,
        280,
        298.7932726892200,
        298.8542359731092,
        299.0963719776158,
        310.3869181829849,
        320.8092442400218,
        330.4305057439312,
        339.3141330683599,
        347.5187109795831,
        355.0974746004750,
        362.0982641139255,
        368.5637676409140,
        374.5319160797160,
        380.0363373972294,
        385.1068144062049,
        389.7697157090284,
        394.0483857232521,
        397.9634892252091,
        401.5333110360412,
        404.7740140374675,
        407.6998597535017,
        410.3233959505448,
        412.6556154950465,
        414.7060902955511,
        416.4830836719433,
        417.9936440064313,
        419.2436820708706,
        420.2380340067240,
        420.9805115599191,
        421.4739408401032,
        421.7201905768194,
        280,
        280,
        297.0853527086431,
        297.1425277301194,
        297.3697761378529,
        307.9049670640065,
        317.5784412485847,
        326.4662830868824,
        334.6390896269612,
        342.1607102849885,
        349.0879350561319,
        355.4707778086129,
        361.3530590416752,
        366.7730945997921,
        371.7643821226003,
        376.3562330865834,
        380.5743306063034,
        384.4412100326086,
        387.9766673052322,
        391.1981031216632,
        394.1208116571502,
        396.7582221084737,
        399.1221004117918,
        401.2227174392475,
        403.0689889680825,
        404.6685918049024,
        406.0280596560596,
        407.1528616597604,
        408.0474659241699,
        408.7153899326404,
        409.1592392662678,
        409.3807357402417,
        280,
        280,
        295.3202577065899,
        295.3732340627008,
        295.5840339624822,
        305.2847326866130,
        314.1332706034377,
        322.2170957280632,
        329.6152320676422,
        336.3971054773060,
        342.6227775304798,
        348.3438530227471,
        353.6045961174104,
        358.4430211552169,
        362.8918634068436,
        366.9794052112720,
        370.7301635970406,
        374.1654564956942,
        377.3038668414951,
        380.1616224882767,
        382.7529073610470,
        385.0901166115175,
        387.1840661489672,
        389.0441648821355,
        390.6785563326951,
        392.0942349235904,
        393.2971411532105,
        394.2922389880077,
        395.0835780976212,
        395.6743429802711,
        396.0668905521524,
        396.2627773777041,
        280,
        280,
        293.5021863484260,
        293.5504628020033,
        293.7429880019498,
        302.5166591165349,
        310.4528127504995,
        317.6535971542997,
        324.2076374382668,
        330.1897020261419,
        335.6622165657627,
        340.6772606345137,
        345.2784512500310,
        349.5025304968499,
        353.3806451383138,
        356.9393607546674,
        360.2014620749403,
        363.1865855116797,
        365.9117210768247,
        368.3916126289491,
        370.6390786872910,
        372.6652708276294,
        374.4798826904716,
        376.0913196076798,
        377.5068365570201,
        378.7326504036199,
        379.7740310452506,
        380.6353750415051,
        381.3202644981022,
        381.8315133387366,
        382.1712025844329,
        382.3407058407845,
        280,
        280,
        291.6358385366847,
        291.6787427689105,
        291.8507935241352,
        299.5861030270870,
        306.5077241277349,
        312.7368427003792,
        318.3720185049930,
        323.4918486232604,
        328.1591260719435,
        332.4245216995422,
        336.3294177513785,
        339.9080044438666,
        343.1888258949226,
        346.1959305941720,
        348.9497384364018,
        351.4677023992956,
        353.7648193252220,
        355.8540282634513,
        357.7465239315856,
        359.4520053212850,
        360.9788741776572,
        362.3343943011393,
        363.5248198841331,
        364.5554990886659,
        365.4309575827142,
        366.1549656347073,
        366.7305915145933,
        367.1602432921878,
        367.4457006061056,
        367.5881375602638,
        280,
        280,
        289.7265864927178,
        289.7630413366547,
        289.9118526822262,
        296.4692353399523,
        302.2551380329837,
        307.4140310144984,
        312.0517452580750,
        316.2465478899726,
        320.0579173992193,
        323.5322823403613,
        326.7066936121027,
        329.6112436323434,
        332.2707234033666,
        334.7057972907245,
        336.9338586772277,
        338.9696663239074,
        340.8258255613449,
        342.5131571680770,
        344.0409834543437,
        345.4173523482995,
        346.6492143977801,
        347.7425635351343,
        348.7025495897544,
        349.5335684842436,
        350.2393345622805,
        350.8229384000637,
        351.2868926334233,
        351.6331677093629,
        351.8632189875853,
        351.9780062339487,
        280,
        280,
        287.7808796048139,
        287.8086450188225,
        287.9305021339589,
        293.1238476175218,
        297.6295616497584,
        301.6123980665653,
        305.1743836228457,
        308.3846802793452,
        311.2937132946091,
        313.9399966505906,
        316.3538307243371,
        318.5595530700661,
        320.5770267955038,
        322.4226764881604,
        324.1102326579053,
        325.6512786577898,
        327.0556594282017,
        328.3317913931964,
        329.4869003694413,
        330.5272062198173,
        331.4580675300580,
        332.2840958518919,
        333.0092464555351,
        333.6368906963022,
        334.1698737821289,
        334.6105607698720,
        334.9608729097020,
        335.2223159242835,
        335.3960014009543,
        335.4826621540332,
        280,
        280,
        285.8074073029015,
        285.8204762718551,
        285.9095904165244,
        289.4660913464272,
        292.5268628819721,
        295.2316159791682,
        297.6487108874065,
        299.8240763099629,
        301.7922588492908,
        303.5801602430644,
        305.2090795645984,
        306.6961111280896,
        308.0551542204318,
        309.2976492085177,
        310.4331168084528,
        311.4695562211720,
        312.4137421005038,
        313.2714486070940,
        314.0476204104361,
        314.7465046314989,
        315.3717536507708,
        315.9265058868686,
        316.4134496842203,
        316.8348740633296,
        317.1927091000892,
        317.4885579876219,
        317.7237223112576,
        317.8992216771434,
        318.0158085379436,
        318.0739788272249,
        280,
        280,
        283.8208660320356,
        283.8003738507832,
        283.8324420749964,
        285.3040644696999,
        286.7801825525440,
        288.1384920807384,
        289.3647676376585,
        290.4706552238186,
        291.4710855495363,
        292.3793059077872,
        293.2062161629120,
        293.9606576572714,
        294.6498297496256,
        295.2796493170351,
        295.8550291462259,
        296.3800873179511,
        296.8583041455571,
        297.2926405242492,
        297.6856280337198,
        298.0394382449809,
        298.3559365546673,
        298.6367243606010,
        298.8831723311575,
        299.0964467727161,
        299.2775305672861,
        299.4272397692827,
        299.5462366705867,
        299.6350399351153,
        299.6940322464799,
        299.7234657897262,
        280,
        280,
        281.8548169424222,
        281.7435847952311,
        280.9031166615437,
        280.1375419048413,
        280.1513107777751,
        280.1774021535923,
        280.2012123586798,
        280.2226913981261,
        280.2421222172579,
        280.2597616756456,
        280.2758215220003,
        280.2904735884680,
        280.3038578037732,
        280.3160891637808,
        280.3272631414740,
        280.3374597588589,
        280.3467466395335,
        280.3551813106355,
        280.3628129552224,
        280.3696837600471,
        280.3758299623257,
        280.3812826697200,
        280.3860685071023,
        280.3902101291005,
        280.3937266270659,
        280.3966338516456,
        280.3989446667006,
        280.4006691462608,
        280.4018147231442,
        280.4023862954831,
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
        280};

    // test

    Fenotype_t fen = Fenotype::createLeftConductorStripFenotype(boardHeight, boardWidth, 2);
    std::vector<cell_type_t> fenotypes(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];

    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .DRAIN_ALPHA = 100,
        .CONDUCTOR_ALPHA = 100,
        .GENERATOR_ALPHA = 1,
        .CONDUCTOR_BETA = 0,
        .GENERATOR_BETA = 1,
        .ETA = usedEta};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(fenotypes, minTemperatures, &equilibriumStep);

    // assert

    ASSERT_EQ(result.size(), 1u);
    ASSERT_IN_REL_ERROR(expectedFitness, result[0], TEST_REL_ERROR);

    for (int row = 0; row < boardHeight; row++)
    {
        for (int column = 0; column < boardWidth; column++)
        {
            ASSERT_IN_REL_ERROR(expectedMinTemperatures[row * boardWidth + column], minTemperatures[row * boardWidth + column], TEST_REL_ERROR);
        }
    }
    ASSERT_IN_REL_ERROR(expectedEquilibriumMoment, equilibriumStep, TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR);
    delete[] minTemperatures;
}
#pragma endregion 32x32LeftStrip

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
