#include <gtest/gtest.h>
#include <cmath>
#include "SequentialHeatSimulation.h"
#include "SystemLayout.h"
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
TEST(SequentialHeatSimulationTest, EvaluateGenerationWith6x6GeneratorFullBottomDrainSystemLayout)
{

    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEta = TEST_REL_ERROR / 1000;
    const simulation_value_t expectedFitness = -281.398877;
    const simulation_steps_index_t expectedEquilibriumMoment = 2'855;

    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        281.398877,
        281.398877,
        281.398877,
        281.398877,
        280,
        280,
        281.298970,
        281.298970,
        281.298970,
        281.298970,
        280,
        280,
        281.099150,
        281.099150,
        281.099150,
        281.099150,
        280,
        280,
        280.799400,
        280.799400,
        280.799400,
        280.799400,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
    };

    // test

    SystemLayout_t fen = SystemLayout::createGeneratorSystemLayout(boardHeight, boardWidth);
    std::vector<int> systemLayouts(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];
    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    std::vector<simulation_value_t> kValues, invCValues, qGenValues;
    createPerCellArrays(systemLayouts.data(), boardHeight * boardWidth,
                        10, 1.0, 1,
                        100, 1.0, 0,
                        0, 1.0, 0,
                        100, 0, 0,
                        kValues, invCValues, qGenValues);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .ETA = usedEta,
        .k_values = kValues.data(),
        .invC_values = invCValues.data(),
        .qGen_values = qGenValues.data()};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(kValues.data(), invCValues.data(), qGenValues.data(), 1, minTemperatures, &equilibriumStep);

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
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithMultiple6x6GeneratorFullBottomDrainSystemLayout)
{

    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEta = TEST_REL_ERROR / 1000;
    const simulation_value_t expectedFitness = -281.398877;
    const simulation_steps_index_t expectedEquilibriumMoment = 2'855;

    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        281.398877,
        281.398877,
        281.398877,
        281.398877,
        280,
        280,
        281.298970,
        281.298970,
        281.298970,
        281.298970,
        280,
        280,
        281.099150,
        281.099150,
        281.099150,
        281.099150,
        280,
        280,
        280.799400,
        280.799400,
        280.799400,
        280.799400,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
    };

    // test

    SystemLayout_t fen = SystemLayout::createGeneratorSystemLayout(boardHeight, boardWidth);
    std::vector<int> systemLayouts(fen, fen + boardHeight * boardWidth);

    for (int i = 0; i < 9; i++)
    {
        systemLayouts.insert(systemLayouts.end(), fen, fen + boardHeight * boardWidth);
    }
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];
    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    int individualsNumber = 10;
    int singleSize = boardHeight * boardWidth;
    std::vector<simulation_value_t> kValues(singleSize * individualsNumber), invCValues(singleSize * individualsNumber), qGenValues(singleSize * individualsNumber);
    for (int i = 0; i < individualsNumber; i++)
    {
        std::vector<simulation_value_t> singleK, singleInvC, singleQGen;
        createPerCellArrays(systemLayouts.data() + i * singleSize, singleSize,
                            10, 1.0, 1,
                            100, 1.0, 0,
                            0, 1.0, 0,
                            100, 0, 0,
                            singleK, singleInvC, singleQGen);
        std::copy(singleK.begin(), singleK.end(), kValues.begin() + i * singleSize);
        std::copy(singleInvC.begin(), singleInvC.end(), invCValues.begin() + i * singleSize);
        std::copy(singleQGen.begin(), singleQGen.end(), qGenValues.begin() + i * singleSize);
    }

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .ETA = usedEta,
        .k_values = kValues.data(),
        .invC_values = invCValues.data(),
        .qGen_values = qGenValues.data()};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(kValues.data(), invCValues.data(), qGenValues.data(), individualsNumber, minTemperatures, &equilibriumStep);

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
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithGeneratorSystemLayoutLowDrainAlpha)
{
    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEta = TEST_REL_ERROR / 1000;
    const simulation_value_t expectedFitness = -281.398877;
    const simulation_steps_index_t expectedEquilibriumMoment = 2'855;

    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        281.398877,
        281.398877,
        281.398877,
        281.398877,
        280,
        280,
        281.298970,
        281.298970,
        281.298970,
        281.298970,
        280,
        280,
        281.099150,
        281.099150,
        281.099150,
        281.099150,
        280,
        280,
        280.799400,
        280.799400,
        280.799400,
        280.799400,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
    };

    // test

    SystemLayout_t fen = SystemLayout::createGeneratorSystemLayout(boardHeight, boardWidth);
    std::vector<int> systemLayouts(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];
    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    std::vector<simulation_value_t> kValues, invCValues, qGenValues;
    createPerCellArrays(systemLayouts.data(), boardHeight * boardWidth,
                        10, 1.0, 1,
                        100, 1.0, 0,
                        0, 1.0, 0,
                        10, 0, 0,
                        kValues, invCValues, qGenValues);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .ETA = usedEta,
        .k_values = kValues.data(),
        .invC_values = invCValues.data(),
        .qGen_values = qGenValues.data()};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(kValues.data(), invCValues.data(), qGenValues.data(), 1, minTemperatures, &equilibriumStep);

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
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithGeneratorSystemLayoutAndStartTemperatures)
{
    // prepare
    int boardHeight = 6;
    int boardWidth = 6;

    const simulation_value_t usedEta = TEST_REL_ERROR / 1000;
    const simulation_value_t expectedFitness = -281.398879;
    const simulation_steps_index_t expectedEquilibriumMoment = 2'506;
    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280,
        280,
        280,
        280,
        280,
        280,
        280,
        281.398879,
        281.398879,
        281.398879,
        281.398879,
        280,
        280,
        281.298972,
        281.298972,
        281.298972,
        281.298972,
        280,
        280,
        281.099151,
        281.099151,
        281.099151,
        281.099151,
        280,
        280,
        280.799401,
        280.799401,
        280.799401,
        280.799401,
        280,
        280,
        280,
        280,
        280,
        280,
        280,
    };

    // test

    SystemLayout_t fen = SystemLayout::createGeneratorSystemLayout(boardHeight, boardWidth);
    std::vector<int> systemLayouts(fen, fen + boardHeight * boardWidth);
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

    std::vector<simulation_value_t> kValues, invCValues, qGenValues;
    createPerCellArrays(systemLayouts.data(), boardHeight * boardWidth,
                        10, 1.0, 1,
                        100, 1.0, 0,
                        0, 1.0, 0,
                        10, 0, 0,
                        kValues, invCValues, qGenValues);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .ETA = usedEta,
        .k_values = kValues.data(),
        .invC_values = invCValues.data(),
        .qGen_values = qGenValues.data()};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(kValues.data(), invCValues.data(), qGenValues.data(), 1, minTemperatures, &equilibriumStep);

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
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithLeftConductorStripSystemLayout)
{
    // prepare
    int boardHeight = 16;
    int boardWidth = 8;

    const simulation_value_t usedEta = TEST_REL_ERROR / 10'000;
    const simulation_value_t expectedFitness = -288.218395;
    const simulation_steps_index_t expectedEquilibriumMoment = 20'710;
    const simulation_value_t expectedMinTemperatures[boardHeight * boardWidth] = {
        280, 280, 280, 280, 280, 280, 280, 280,
        280, 282.254585, 282.274154, 282.351670, 285.285886, 287.241095, 288.218395, 280,
        280, 282.235015, 282.254577, 282.332068, 285.264948, 287.219086, 288.195792, 280,
        280, 282.195885, 282.215430, 282.292863, 285.222807, 287.174590, 288.149989, 280,
        280, 282.137209, 282.156726, 282.234057, 285.158882, 287.106559, 288.079679, 280,
        280, 282.059016, 282.078487, 282.155653, 285.072158, 287.013162, 287.982578, 280,
        280, 281.961352, 281.980750, 282.057655, 284.960987, 286.891425, 287.854980, 280,
        280, 281.844291, 281.863574, 281.940071, 284.822756, 286.736643, 287.691018, 280,
        280, 281.707946, 281.727050, 281.802909, 284.653368, 286.541438, 287.481506, 280,
        280, 281.552498, 281.571322, 281.646187, 284.446409, 286.294295, 287.212132, 280,
        280, 281.378226, 281.396611, 281.469925, 284.191824, 285.977254, 286.860657, 280,
        280, 281.185569, 281.203261, 281.274154, 283.873739, 285.562288, 286.392640, 280,
        280, 280.975220, 280.991803, 281.058898, 283.466716, 285.005558, 285.755023, 280,
        280, 280.748287, 280.763046, 280.824128, 282.928690, 284.238240, 284.866907, 280,
        280, 280.506596, 280.518202, 280.569415, 282.185692, 283.151829, 283.607487, 280,
        280, 280, 280, 280, 280, 280, 280, 280,
    };

    // test

    SystemLayout_t fen = SystemLayout::createLeftConductorStripSystemLayout(boardHeight, boardWidth, 2);
    std::vector<int> systemLayouts(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];

    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    std::vector<simulation_value_t> kValues, invCValues, qGenValues;
    createPerCellArrays(systemLayouts.data(), boardHeight * boardWidth,
                        1, 1.0, 1,
                        100, 1.0, 0,
                        0, 1.0, 0,
                        100, 0, 0,
                        kValues, invCValues, qGenValues);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .ETA = usedEta,
        .k_values = kValues.data(),
        .invC_values = invCValues.data(),
        .qGen_values = qGenValues.data()};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(kValues.data(), invCValues.data(), qGenValues.data(), 1, minTemperatures, &equilibriumStep);

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
TEST(SequentialHeatSimulationTest, Evaluate4RotatedGenerationsWithLeftConductorStripSystemLayout)
{
    // prepare
    int boardLength = 8;

    const simulation_value_t usedEta = TEST_REL_ERROR / 100'000;
    const simulation_value_t expectedFitness = -285.912264;
    const simulation_steps_index_t expectedEquilibriumMoment = 17'525;
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
        280.457066,
        280.475390,
        280.548474,
        283.258701,
        285.034440,
        285.912264,
        280,
        280,
        280.438742,
        280.456807,
        280.528980,
        283.193194,
        284.932363,
        285.790097,
        280,
        280,
        280.402353,
        280.419841,
        280.490005,
        283.052738,
        284.711727,
        285.525674,
        280,
        280,
        280.348476,
        280.364930,
        280.431563,
        282.816030,
        284.336140,
        285.075206,
        280,
        280,
        280.278146,
        280.292824,
        280.353619,
        282.443683,
        283.741602,
        284.363812,
        280,
        280,
        280.193136,
        280.204695,
        280.255740,
        281.863485,
        282.822777,
        283.274634,
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

    const simulation_value_t *expectedLeftwardsMinTemperatures = SystemLayout::rotateLeftSquareSystemLayoutBy90Deg(boardLength, expectedUpwardsMinTemperatures);
    const simulation_value_t *expectedDownwardsMinTemperatures = SystemLayout::rotateLeftSquareSystemLayoutBy90Deg(boardLength, expectedLeftwardsMinTemperatures);
    const simulation_value_t *expectedRightwardsMinTemperatures = SystemLayout::rotateLeftSquareSystemLayoutBy90Deg(boardLength, expectedDownwardsMinTemperatures);

    const simulation_value_t *expectedMinTemperatures[4] = {
        expectedUpwardsMinTemperatures,
        expectedLeftwardsMinTemperatures,
        expectedDownwardsMinTemperatures,
        expectedRightwardsMinTemperatures,
    };

    // test

    SystemLayout_t upFen = SystemLayout::createLeftConductorStripSystemLayout(boardLength, boardLength, 2);
    std::vector<int> upwardsSystemLayouts(upFen, upFen + boardLength * boardLength);

    SystemLayout_t leftFen = SystemLayout::rotateLeftSquareSystemLayoutBy90Deg(boardLength, upFen);
    std::vector<int> leftwardsSystemLayouts(leftFen, leftFen + boardLength * boardLength);

    SystemLayout_t downFen = SystemLayout::rotateLeftSquareSystemLayoutBy90Deg(boardLength, leftFen);
    std::vector<int> downwardsSystemLayouts(downFen, downFen + boardLength * boardLength);

    SystemLayout_t rightFen = SystemLayout::rotateLeftSquareSystemLayoutBy90Deg(boardLength, downFen);
    std::vector<int> rightwardsSystemLayouts(rightFen, rightFen + boardLength * boardLength);

    std::vector<int> systemLayouts[4] = {
        upwardsSystemLayouts,
        leftwardsSystemLayouts,
        downwardsSystemLayouts,
        rightwardsSystemLayouts,
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
        .ETA = usedEta,
        .k_values = nullptr,
        .invC_values = nullptr,
        .qGen_values = nullptr};

    SequentialHeatSimulation simulator(boardLength, boardLength, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardLength * boardLength];
    simulation_steps_index_t equilibriumStep;

    for (int index = 0; index < 4; index++)
    {
        std::vector<simulation_value_t> kValues, invCValues, qGenValues;
        createPerCellArrays(systemLayouts[index].data(), boardLength * boardLength,
                            1, 1.0, 1,
                            100, 1.0, 0,
                            0, 1.0, 0,
                            100, 0, 0,
                            kValues, invCValues, qGenValues);

        auto result = simulator.evaluateGeneration(kValues.data(), invCValues.data(), qGenValues.data(), 1, minTemperatures, &equilibriumStep);

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
TEST(SequentialHeatSimulationTest, EvaluateGenerationWithBigLeftConductorStripSystemLayout)
{
    // prepare
    int boardHeight = 32;
    int boardWidth = 32;

        const simulation_value_t usedEta = TEST_REL_ERROR / 100'000;
    const simulation_value_t expectedFitness = -552.445548;
    const simulation_steps_index_t expectedEquilibriumMoment = 714'610;
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
        321.78378648014319,
        321.87839662362717,
        322.25329772934418,
        340.28093448711326,
        357.42471254148927,
        373.70565014495429,
        389.14467940678537,
        403.76255377796838,
        417.5797579441047,
        430.61642094773458,
        442.8922332565225,
        454.42636837307833,
        465.23740945208965,
        475.34328125711147,
        484.76118765856171,
        493.50755475119053,
        501.59797955748394,
        509.04718418605188,
        515.86897523282539,
        522.07620814865186,
        527.68075624962353,
        532.69348401544494,
        537.12422430515392,
        540.98175911701696,
        544.27380352870568,
        547.0069924732295,
        549.18687003380194,
        550.8178809753307,
        551.90336427011596,
        552.44554841938134,
        280,
        280,
        321.68917635666583,
        321.78368172877924,
        322.15816883511746,
        340.16479618452428,
        357.28755694715875,
        373.54756339093365,
        388.9658401375325,
        403.56323074239589,
        417.36030676620874,
        430.37728018078724,
        442.63391984149621,
        454.14947263140283,
        464.94258974607021,
        475.03125844880248,
        484.43273949038269,
        493.16351026008516,
        501.23921362003273,
        508.67461227577263,
        515.48354845388178,
        521.678908593237,
        527.27259271004084,
        532.27548806687889,
        536.69744676162759,
        540.54726685123592,
        543.83267663644085,
        546.56032175452128,
        548.73575475647726,
        550.36342688093544,
        551.44668177822734,
        551.9877509831756,
        280,
        280,
        321.50006088103191,
        321.59435623734942,
        321.96801348136097,
        339.93252745528122,
        357.01315961661294,
        373.23121122623974,
        388.60789283559222,
        403.16422903044776,
        416.92096583809882,
        429.89848168468973,
        442.11670266665806,
        453.5950227603318,
        464.3522294445753,
        474.40643506045285,
        483.77501408594384,
        492.47454636841519,
        500.52076623708734,
        507.92851731497586,
        514.71171276627888,
        520.88330065123944,
        526.45523401559512,
        531.43844531450998,
        535.84282475965176,
        539.67720218080149,
        542.94933200800688,
        545.66588100477691,
        547.83241841526615,
        549.45340822716639,
        550.53220329561464,
        551.07104112065906,
        280,
        280,
        321.21665006893784,
        321.31062921391708,
        321.68303698600613,
        339.58414351032593,
        356.60134676305478,
        372.75623392961114,
        388.0702967447188,
        402.56483341470414,
        416.2608534735561,
        429.17898652778604,
        441.3393957028477,
        452.76169644329104,
        463.4648811493002,
        473.46724996274725,
        482.7863478541845,
        491.43890801455012,
        499.44080142625927,
        506.80699238055655,
        513.55149962265409,
        519.68736274308151,
        525.22661339166791,
        530.18025086769342,
        534.55822163374648,
        538.36940230915707,
        541.62158571909549,
        544.3214696052562,
        546.47464764134156,
        548.0856024396187,
        549.15770028213115,
        549.69318736046478,
        280,
        280,
        320.83926013157338,
        320.9328154303156,
        321.30354845151749,
        339.11966578640329,
        356.05185389154786,
        372.12208581575919,
        387.35223255211747,
        401.76396106883806,
        415.37863565919423,
        428.21722366114057,
        440.30020642655421,
        451.64749622929344,
        462.27835960250843,
        472.21134739966493,
        481.46423169004919,
        490.05394943503768,
        497.99655275094347,
        505.30716545026394,
        511.9999454661401,
        518.08805270379321,
        523.5836218254899,
        528.49773945939216,
        532.84042532432852,
        536.62061677895281,
        539.84615633189196,
        542.52378168668986,
        544.65911793932105,
        546.25667159485192,
        547.31982612202808,
        547.85083881897049,
        280,
        280,
        320.3683149149777,
        320.46133678930784,
        320.82996212750834,
        338.53912021206958,
        355.36432105776242,
        371.32802767261114,
        386.45259227457717,
        400.76014924126792,
        414.27251190308959,
        427.01107435772656,
        438.99671927299113,
        450.2497324124023,
        460.78972437876962,
        470.63555983944866,
        479.80529428407937,
        488.31611817960646,
        496.18430823307006,
        503.42518535193932,
        510.05307880413352,
        516.08129602299778,
        521.52209747290499,
        526.38667598476172,
        530.68513998324102,
        534.42650005506528,
        537.61865834639241,
        540.26840032426355,
        542.38138848955509,
        543.96215768497632,
        545.01411169985681,
        545.53952093276837,
        280,
        280,
        319.80434784331658,
        319.89672424615048,
        320.26279924780914,
        337.8425347595969,
        354.53828626292886,
        370.37311626482631,
        385.36996525584698,
        399.55153822729454,
        412.9401957297319,
        425.55785081527893,
        437.425872939719,
        448.56499961032506,
        458.99525627206344,
        468.73588464627039,
        477.80527948583187,
        486.22093349847552,
        493.99939001973405,
        501.15620289022922,
        507.70590290599063,
        513.66197016130241,
        519.03681158567679,
        523.84174298418134,
        528.08697491783039,
        531.78160180293366,
        534.93359366076788,
        537.54979000808373,
        539.63589544196645,
        541.19647653742265,
        542.23495974131959,
        542.75363001115477,
        280,
        280,
        319.14800438779207,
        319.23962024900618,
        319.60269036005764,
        337.02993615454949,
        353.57317671948755,
        369.25619051828693,
        384.10261979435387,
        398.13584909167093,
        411.37888923620329,
        423.854268329989,
        435.58393096660939,
        446.5891465086022,
        456.8904269021981,
        466.50745416539758,
        475.45901738898556,
        483.76295884654519,
        491.43612862294583,
        498.49434703987146,
        504.95237407694987,
        510.82388495087361,
        516.12145101452188,
        520.85652516540949,
        525.03943100033814,
        528.67935501445425,
        531.784341212698,
        534.3612875758464,
        536.415943898921,
        537.95291059473573,
        538.97563812850967,
        539.4864268202906,
        280,
        280,
        318.40004508968269,
        318.49078172566936,
        318.85037817133957,
        336.101345566943,
        352.46829762466876,
        367.97585386117328,
        382.648479749633,
        396.51035540284721,
        409.58525092567066,
        421.89641025224097,
        433.4664448342017,
        444.31723807238285,
        454.46986092402057,
        463.94449870067285,
        472.7603887189901,
        480.93576818799193,
        488.48783151452005,
        495.43269607845212,
        501.7853754622559,
        507.55975910438178,
        512.76859737121003,
        517.42349109674774,
        521.53488471335004,
        525.1120621824856,
        528.16314502578427,
        530.69509184883498,
        532.71369884037711,
        534.22360081587397,
        535.22827245606527,
        535.73002946794224,
        280,
        280,
        317.56134917382803,
        317.65108365154197,
        318.00672093955399,
        335.05677304708655,
        351.22281795692356,
        366.53045202376666,
        381.00509526501708,
        394.67184800745287,
        407.55535579522586,
        419.67968470403457,
        431.06820860964569,
        441.7435093417642,
        451.7272900685864,
        461.04030174240319,
        469.70228201647518,
        477.73190572791879,
        485.14674582853331,
        491.96324352513386,
        498.19668634803287,
        503.86119288399362,
        508.96970297181377,
        513.53397224999605,
        517.56457005452683,
        521.07087978111201,
        524.06110094335986,
        526.54225227169002,
        528.52017530432806,
        529.99953802029665,
        530.98383815449813,
        531.47540591738982,
        280,
        280,
        316.63291879806798,
        316.72152324253403,
        317.072696447216,
        333.89621038999519,
        349.83575265242695,
        364.91804537746566,
        379.16960647767212,
        392.61659158345947,
        405.2846463617945,
        417.19877175915587,
        428.38320391935849,
        438.8613097142802,
        448.65549807539844,
        457.78714667685506,
        466.27654302370541,
        474.14283864849102,
        481.40401490587954,
        488.07685875947936,
        494.17694695293005,
        499.71863702426111,
        504.71506373566433,
        509.17813963111359,
        513.11855858686113,
        516.54580137365417,
        519.46814239728099,
        521.89265692122376,
        523.82522819956125,
        525.27055405927467,
        526.23215256936214,
        526.71236652108189,
        280,
        280,
        315.61588399517308,
        315.70322482700118,
        316.04940660259865,
        332.61962200689118,
        348.30594031122843,
        363.13637460450423,
        377.13869874387257,
        390.34027134239471,
        402.76787294464799,
        414.44755944786652,
        425.40453373104071,
        435.6630363745964,
        445.2462553882051,
        454.17625407765604,
        462.47391560135088,
        470.15890239082063,
        477.24962841507772,
        483.76324222181495,
        489.71561875211455,
        495.12135806408725,
        499.99378928437039,
        504.34497831070593,
        508.18573799629593,
        511.52563974539618,
        514.37302563245521,
        516.73502032036549,
        518.61754219612612,
        520.02531326505698,
        520.9618674496312,
        521.42955702852305,
        280,
        280,
        314.5115083772609,
        314.59744546347133,
        314.93808272275368,
        331.22693323950369,
        346.63201530412698,
        361.18281810638859,
        374.90854745813681,
        387.83792777725455,
        399.99902106260987,
        411.4190665309618,
        422.12434307471148,
        432.14005525299916,
        441.49024228485598,
        450.19770854932921,
        458.28397343583043,
        465.76923800878484,
        472.67236580869451,
        479.01087515105667,
        484.80094044928148,
        490.05740032868118,
        494.79377057648338,
        499.02226025866867,
        502.75378960872519,
        505.9980085443301,
        508.76331488949063,
        511.05687157012346,
        512.8846222113259,
        514.25130469743442,
        515.16046236498141,
        515.61445258789922,
        280,
        280,
        313.32119568938685,
        313.40558138650238,
        313.7400915644111,
        329.71801535587008,
        344.81237277129151,
        359.05433904190875,
        372.47474994813552,
        385.10387673558847,
        396.97122321840811,
        408.10534947306786,
        418.5337244122706,
        428.28260757875893,
        437.37695889890671,
        445.84037397283521,
        453.69504175455955,
        460.96172113865265,
        467.65973293658249,
        473.80696390721954,
        479.41987982088381,
        484.51354491874315,
        489.1016455305421,
        493.19651600060445,
        496.80916542500006,
        499.94930401175822,
        502.62536813805195,
        504.84454339484898,
        506.61278508469866,
        507.93483577752295,
        508.81423963835647,
        509.25335332577015,
        280,
        280,
        312.04649732004327,
        312.12917538248917,
        312.45694218250236,
        328.09266618972021,
        342.84512447604612,
        356.74741917722127,
        369.83224112410858,
        382.13161128467266,
        393.67665159249037,
        404.49739040786648,
        414.62260486805849,
        424.07969974393382,
        432.8946203772677,
        441.09179590734516,
        448.69410826442646,
        455.72288219527394,
        462.19789175027262,
        468.13737906615654,
        473.55808180938288,
        478.47526621796356,
        482.90276323696605,
        486.85300575079805,
        490.33706535690493,
        493.3646874955939,
        495.94432405176923,
        498.08316278298582,
        499.78715311288869,
        501.06102796879213,
        501.9083214454605,
        502.33138215198221,
        280,
        280,
        310.68912090325824,
        310.76992522195621,
        311.09029371261761,
        326.35058498990531,
        340.7280427319281,
        354.2579757449667,
        366.97518846634301,
        378.91368075588127,
        390.1063872035337,
        400.58496210130801,
        410.37961195250017,
        419.51897381699717,
        428.03003522362076,
        435.93808985566614,
        443.26672259248494,
        450.03781754435522,
        456.27158321637154,
        461.98658967835752,
        467.19981344975764,
        471.92668662869528,
        476.18114754213468,
        479.9756908398229,
        483.32141548980303,
        486.22806956207205,
        488.70409102046182,
        490.75664399502477,
        492.3916501923797,
        493.61381523301634,
        494.42664979422653,
        494.83248549518362,
        280,
        280,
        309.25094018210103,
        309.32969330669761,
        309.64196419372036,
        324.49133946921319,
        338.45848854858735,
        351.58125611611553,
        363.89686042220103,
        375.44154090900997,
        386.25025984908683,
        396.35646495509178,
        405.79191374928081,
        414.58655566649253,
        422.76846473538302,
        430.36381414006121,
        437.3968836725482,
        443.89009164126287,
        449.86404383468687,
        455.33759336950078,
        460.32790648782185,
        464.85053049654579,
        468.91946100950315,
        472.54720644526486,
        475.74484835790736,
        478.52209665441939,
        480.88733910445825,
        482.84768479960843,
        484.40900139163926,
        485.57594605086945,
        486.3519901517235,
        486.73943772501042,
        280,
        280,
        307.73400634995841,
        307.81051773173346,
        308.11394057154808,
        322.51432218157987,
        336.0333185676605,
        348.71170308531339,
        360.58946017056996,
        371.70536720752597,
        382.09665153981769,
        391.79872992955018,
        400.84502881326603,
        409.26687731788155,
        417.09346140866251,
        424.35182631653862,
        431.06691483624115,
        437.26163050925402,
        442.95691655796765,
        448.17184334741205,
        452.92369890183886,
        457.22807849364233,
        461.09897052467909,
        464.54883685062765,
        467.58868639332712,
        470.22814138622761,
        472.4754959448714,
        474.33776688353032,
        475.82073684026557,
        476.92898984927308,
        477.66593952953218,
        478.03385005620578,
        280,
        280,
        306.14056114890371,
        306.21462501733214,
        306.50839005182002,
        320.41869204285609,
        333.44876299781339,
        345.64278064008204,
        357.04391372204282,
        367.69382055654802,
        377.63225409773031,
        386.89677989953861,
        395.52260029527224,
        403.5424699544011,
        410.98668434998245,
        417.88312246018268,
        424.25732689809689,
        430.13260750278823,
        435.53015746755585,
        440.46917388799216,
        444.96697698043636,
        449.03912410040363,
        452.69951611243994,
        455.96049469609181,
        458.83292989469754,
        461.32629769695382,
        463.44874774700719,
        465.20716145620031,
        466.60720087607547,
        467.65334871568308,
        468.34893986687666,
        468.69618475341713,
        280,
        280,
        304.47305209150323,
        304.54444484105505,
        304.82767300555724,
        318.20329474840099,
        330.70026312907152,
        342.36674571706016,
        353.24959704800665,
        363.39375128111482,
        372.84176902091042,
        381.63354043211979,
        389.80612818648166,
        397.39372402714747,
        404.42769023203329,
        410.93665939536214,
        416.94667035624423,
        422.48132312180229,
        427.56194030699652,
        432.20772651839792,
        436.43592014494504,
        440.26193425442506,
        443.69948486718675,
        446.76070593700769,
        449.45625104639294,
        451.79538222850658,
        453.7860465437476,
        455.43494112700586,
        456.74756742555644,
        457.72827529765902,
        458.38029755867223,
        458.70577545833441,
        280,
        280,
        302.73415029582333,
        302.80262720255581,
        303.07435766466818,
        315.86655250304432,
        327.78225128124876,
        338.87434481434485,
        349.19398076230749,
        358.78982230744549,
        367.70753458833912,
        375.9894894322498,
        383.67465328355803,
        390.79861349438812,
        397.39369936465619,
        403.48916117479155,
        409.11137906557383,
        414.28408177119292,
        419.02856194338358,
        423.36387990788711,
        427.30705132877068,
        430.87321671145452,
        434.07579225036545,
        436.92660247742606,
        439.43599569166145,
        441.61294339348683,
        443.46512501173248,
        444.99899916645273,
        446.21986260165039,
        447.13189777835032,
        447.73820995788998,
        448.04085443824829,
        280,
        280,
        300.92677160384488,
        300.99206259788184,
        301.25123687808298,
        313.40630787932565,
        324.68784674105422,
        335.15440405451994,
        344.86216192521204,
        353.86402612323144,
        362.20906158745038,
        369.94223387790356,
        377.10438691969631,
        383.73238263262039,
        389.85933830465257,
        395.5149130214171,
        400.72560949116962,
        405.51506984995564,
        409.90435302879138,
        413.91218740729533,
        417.55519642072045,
        420.848097163609,
        423.80387335518395,
        426.43392467517589,
        428.74819470418134,
        430.75527968224856,
        432.46252014330565,
        433.87607725939841,
        435.00099547780462,
        435.84125277872187,
        436.39979963290716,
        436.67858750222513,
        280,
        280,
        299.0541019273993,
        299.11590598533621,
        299.36134720894586,
        310.81959682724556,
        321.40842564079662,
        331.19326508335598,
        340.23623955427581,
        348.59506190602264,
        356.32245542413204,
        363.4660016562878,
        370.06828237751995,
        376.16719670063577,
        381.79636347109295,
        386.98554875364874,
        391.76108201784092,
        396.14624143339609,
        400.16159955603865,
        403.82532721133765,
        407.15345700123692,
        410.160109643223,
        412.85768704476135,
        415.25703609236155,
        417.36758688900233,
        419.19746877938803,
        420.7536070578322,
        422.0418028108399,
        423.06679793093133,
        423.83232695956013,
        424.34115707614325,
        424.59511724371532,
        280,
        280,
        297.11962820169839,
        297.17760562150744,
        297.40799059244847,
        308.10230787913525,
        317.93299562765498,
        326.9739932120184,
        335.29447183676007,
        342.95752945561054,
        350.01969987056617,
        356.53103865065867,
        362.53554830933717,
        368.07176275650539,
        373.17337490744012,
        377.86984161956974,
        382.18693382736797,
        386.14722004753389,
        389.77048257579679,
        393.07407117633534,
        396.07320127785528,
        398.78120414563489,
        401.20973608566021,
        403.36895295318749,
        405.26765534776814,
        406.91340901038814,
        408.31264415272813,
        409.47073676155281,
        410.39207433139347,
        411.0801079756763,
        411.53739243663608,
        411.76561514356086,
        280,
        280,
        295.1271770639604,
        295.18093827495335,
        295.39475850837306,
        305.2486496322789,
        314.24725731497716,
        322.47524220549673,
        330.01012739384208,
        336.92088683491414,
        343.26777892739165,
        349.10290808328705,
        354.47111310146852,
        359.41093507908056,
        363.95553606353366,
        368.13351356915939,
        371.96959648934421,
        375.4852274901495,
        378.69904491712578,
        381.62727927625372,
        384.28407865027657,
        386.68177564749124,
        388.83110646324855,
        390.74139072598501,
        392.42067913419186,
        393.8758744954344,
        395.11283063398167,
        396.13643270383977,
        396.95066169014115,
        397.5586452681452,
        397.96269668449196,
        398.16434290383472,
        280,
        280,
        293.08096472207251,
        293.13005292837488,
        293.32555778065381,
        302.25027585069222,
        310.33214314714803,
        317.66959257861225,
        324.34991069572197,
        330.44811387476415,
        336.02762354065851,
        341.1417045740057,
        345.8350641468262,
        350.14533189069039,
        354.10432446768942,
        357.73908413615231,
        361.07271535401816,
        364.12505302918714,
        366.9131950754122,
        369.45192732351262,
        371.75406356086842,
        373.83071867673601,
        375.69152890928245,
        377.34483002267103,
        378.79780177488391,
        380.05658513204384,
        381.12637721763645,
        382.01150785122036,
        382.71550064928158,
        383.24112096746052,
        383.5904124100208,
        383.76472318188445,
        280,
        280,
        290.98566417977906,
        291.02952584360929,
        291.20463461980086,
        299.09475372517835,
        306.16144800997876,
        312.52107571164385,
        318.27181065704502,
        323.49403642009759,
        328.25289904414768,
        332.60122504189872,
        336.58210978966122,
        340.23100688175214,
        343.57734902853656,
        346.64578662830371,
        349.45713145255775,
        352.02907809453302,
        354.37675912439192,
        356.51317565774383,
        358.44953404077557,
        360.19551119618387,
        361.75946522741032,
        363.14860356612951,
        364.36911781509934,
        365.42629214945032,
        366.32459045287834,
        367.06772610938555,
        367.6587174243702,
        368.09993092421502,
        368.39311421864107,
        368.53941965950162,
        280,
        280,
        288.84650197859173,
        288.88443049564552,
        289.03658605083808,
        295.76265715896511,
        301.69782043171534,
        306.98145281100597,
        311.72222124166353,
        316.00332377219809,
        319.88871306381367,
        323.42818886644113,
        326.66114540566963,
        329.61923933992722,
        332.32828085539904,
        334.80958480448987,
        337.080948823217,
        339.15737203444064,
        341.05159109572651,
        342.77448572188507,
        344.33538947153824,
        345.74233069618822,
        347.00222121668696,
        348.12100528902124,
        349.10377795912672,
        349.9548794746658,
        350.67797068755704,
        351.27609312494315,
        351.75171648136757,
        352.10677559145688,
        352.34269841148944,
        352.46042612145698,
        280,
        280,
        286.66941126431402,
        286.70042511333622,
        286.82632231209652,
        292.22146902125195,
        296.88572453035647,
        300.98469483056141,
        304.6322988833391,
        307.90832570226394,
        310.87044208984884,
        313.56167364581853,
        316.01504548736915,
        318.25652624164439,
        320.30695243172369,
        322.18332524628858,
        323.89970948221031,
        325.46787274370064,
        326.89775025279852,
        328.19778953657936,
        329.37521041668407,
        330.43620399662348,
        331.38608684857439,
        332.22942169774507,
        332.97011262121885,
        333.61148053643672,
        334.15632319237159,
        334.60696276695819,
        334.9652833710656,
        335.2327601658713,
        335.41048135206745,
        335.49916394133351,
        280,
        280,
        284.46130670399464,
        284.48385549766294,
        284.57872070325038,
        288.41117252963301,
        291.63891442706858,
        294.43930382818331,
        296.91395462891359,
        299.12723907066857,
        301.12305708860976,
        302.93301941162053,
        304.58083805565565,
        306.08486923022838,
        307.45967902529839,
        308.71705602290001,
        309.86669298128578,
        310.91666117522266,
        311.87374970372014,
        312.74371391631081,
        313.53146091008841,
        314.24119035353226,
        314.87650288555182,
        315.44048450152383,
        315.93577282100847,
        316.36460943988101,
        316.72888140659069,
        317.03015404576195,
        317.26969676710735,
        317.44850306906255,
        317.56730562519692,
        317.62658709383504,
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

    SystemLayout_t fen = SystemLayout::createLeftConductorStripSystemLayout(boardHeight, boardWidth, 2);
    std::vector<int> systemLayouts(fen, fen + boardHeight * boardWidth);
    delete[] fen;

    simulation_value_t *startTemperatures = new simulation_value_t[boardHeight * boardWidth];

    std::fill_n(startTemperatures, boardHeight * boardWidth, (simulation_value_t)280);

    std::vector<simulation_value_t> kValues, invCValues, qGenValues;
    createPerCellArrays(systemLayouts.data(), boardHeight * boardWidth,
                        1, 1.0, 1,
                        100, 1.0, 0,
                        0, 1.0, 0,
                        100, 0, 0,
                        kValues, invCValues, qGenValues);

    SimulationParams testParams = {
        .simulationSteps = expectedEquilibriumMoment + (simulation_steps_index_t)(TEST_EQUILIBRIUM_MOMENT_STEPS_REL_ERROR * expectedEquilibriumMoment) + 2,
        .startTemperatures = startTemperatures,
        .drainTemperature = 280,
        .delta_time = 0.003,
        .ETA = usedEta,
        .k_values = kValues.data(),
        .invC_values = invCValues.data(),
        .qGen_values = qGenValues.data()};

    SequentialHeatSimulation simulator(boardHeight, boardWidth, testParams);
    simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
    simulation_steps_index_t equilibriumStep;

    auto result = simulator.evaluateGeneration(kValues.data(), invCValues.data(), qGenValues.data(), 1, minTemperatures, &equilibriumStep);

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
