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
        ParallelHeatSimulation(int boardHeight, int boardWidth, int boardThickness, const SimulationParams &params = SimulationParams{});
        ~ParallelHeatSimulation();

        std::vector<simulation_value_t> evaluateGeneration(const std::vector<cell_type_t> &fenotypes, simulation_value_t *finalTemperatureDistributions = NULL, simulation_steps_index_t *equilibriumMoments = NULL) override;

        void setSimulationParams(const SimulationParams &params) override;

    protected:
        cl_platform_id platform;
        cl_device_id device;
        cl_context context;
        int stripLength, stripsPerColumn;
        size_t kernelMaxWorkGroupSize;
        cl_program program = nullptr;

        inline void buildProgram();
        inline int calculateStripLength(int maxWorkGroupSize);

        inline std::vector<simulation_value_t> runSimulationKernel(
            const std::vector<cell_type_t> &boards, int individualsNumber, simulation_value_t *resultFinalTemperatureDistributions, simulation_steps_index_t *resultEquilibriumMoment);
    };

}

using conductivity_evaluators::ParallelHeatSimulation;
