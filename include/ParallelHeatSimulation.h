#include "GraphicsUtils.h"
#include <vector>
#include <chrono>
#include "Simulation.h"

namespace conductivity_evaluators
{

    /*
    Possible optimizations:
    - strips
    */

    class ParallelHeatSimulation : public Simulation
    {
    public:
        ParallelHeatSimulation(int boardHeight, int boardWidth, const SimulationParams &params = SimulationParams{});
        ~ParallelHeatSimulation();

        std::vector<simulation_value_t> evaluateGeneration(const std::vector<cell_type_t> &systemLayouts, simulation_value_t *minFinalTemperatures = NULL, simulation_steps_index_t *lastEquilibriumMoment = NULL) override;

        void setSimulationParams(const SimulationParams &params) override;

    protected:
        cl_platform_id platform;
        cl_device_id device;
        cl_context context;
        int stripLength, stripsPerColumn;
        size_t kernelMaxWorkGroupSize;
        cl_program program;

        inline void buildProgram();
        inline int calculateStripLength(int maxWorkGroupSize);

        inline std::vector<simulation_value_t> runSimulationKernel(
            const std::vector<cell_type_t> &boards, int individualsNumber, simulation_value_t **returnedMinFinalTemperatures, simulation_steps_index_t *returnedLastEquilibriumMoment);
    };

}

using conductivity_evaluators::ParallelHeatSimulation;
