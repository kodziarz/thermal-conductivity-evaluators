#include <gtest/gtest.h>
#include <cmath>
#include "SequentialHeatSimulation.h"
#include "Fenotype.h"

TEST(SequentialHeatSimulationTest, EvaluateGenerationWithGeneratorFenotype)
{
    int height = 6;
    int width = 6;
    Fenotype_t fen = Fenotype::createGeneratorFenotype(height, width);
    std::vector<cell_type_t> fen_vec(fen, fen + height * width);
    delete[] fen;

    SequentialHeatSimulation sim(height, width, 100'000);

        auto res = sim.evaluateGeneration(fen_vec);
    ASSERT_EQ(res.size(), 1u);
    EXPECT_TRUE(std::isfinite(res[0]));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
