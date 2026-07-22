#include "GraphicsUtils.h"
#include <cstring>
#include <iostream>
#include <CL/cl.h>
#include <vector>
#include <sstream>
#include <chrono>
#include <exception>
#include <utility>
#include "ParallelHeatSimulation.h"
#include <ctime>
#include "env.h"
#include "simulation_kernel_global_stripped_column_wise.hpp" // makes the std::string_view simulation_kernel_global_stripped_column_wise available

// #define KERNEL_DEBUG

namespace conductivity_evaluators
{

    ParallelHeatSimulation::ParallelHeatSimulation(int boardHeight, int boardWidth, const SimulationParams &params) : Simulation(boardHeight, boardWidth, params)
    {
        // ----------------------------------------------------
        // 2. Platform
        // ----------------------------------------------------
        cl_int status;
        cl_uint numPlatforms = 0;
        status = clGetPlatformIDs(0, nullptr, &numPlatforms);
        if (numPlatforms == 0)
        {
            std::cout << "No platforms found..." << std::endl;
            throw std::runtime_error("Platform not found error");
        }

        std::vector<cl_platform_id> platforms(numPlatforms);
        status = clGetPlatformIDs(numPlatforms, platforms.data(), nullptr);

        platform = platforms[0];

        char platform_name[512];
        size_t name_len;
        clGetPlatformInfo(platforms[0], CL_PLATFORM_NAME, sizeof(platform_name), platform_name, &name_len);

        // ----------------------------------------------------
        // 3. Device (try GPU, fallback to CPU)
        // ----------------------------------------------------
#ifdef PREFERRED_DEVICE_GPU
        device = getDevice(platform, CL_DEVICE_TYPE_GPU, CL_DEVICE_TYPE_CPU);
#else
        device = getDevice(platform, CL_DEVICE_TYPE_CPU, CL_DEVICE_TYPE_GPU);
#endif

        // ----------------------------------------------------
        // 4. Context + Command Queue
        // ----------------------------------------------------
        context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
        checkContextCreationErrors(status);

        stripLength = individualHeight;
        stripsPerColumn = 1;

        // ----------------------------------------------------
        // 6. Build Kernel Program
        // ----------------------------------------------------

#ifdef PARALLEL_SIMULATION_BENCHMARK
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
#endif

        buildProgram();

        cl_kernel kernel = clCreateKernel(program, "simulate_heat", &status);
        checkKernelCreationErrors(status);

#ifdef PARALLEL_SIMULATION_BENCHMARK
        std::chrono::time_point<std::chrono::high_resolution_clock> stop = std::chrono::high_resolution_clock::now();
        std::cout << "Kernel compilation took " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;
#endif

        status = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &kernelMaxWorkGroupSize,
                                          nullptr);
        if (status != CL_SUCCESS)
        {
            std::cout << "Failed to get kernel work group size" << std::endl;
            throw std::runtime_error("clGetKernelWorkGroupInfo error");
        }

        stripLength = calculateStripLength(kernelMaxWorkGroupSize);
        stripsPerColumn = individualHeight / stripLength;

        size_t preferedWorkGroupSizeMultiple;
        status = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &preferedWorkGroupSizeMultiple,
                                          nullptr);
        if (status != CL_SUCCESS)
        {
#ifdef PRINT_DEVICE_INFO
            printf("Failed to get CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE\n");
#endif
        }

        if (individualWidth % preferedWorkGroupSizeMultiple != 0)
        {
            std::cout << "Warning: Width of an individual (which is " << individualWidth << ") is adviced to be a mulltiple of prefered work group size multiple (which is for the device and the kernel " << preferedWorkGroupSizeMultiple << ")." << std::endl;
        }

#ifdef PRINT_DEVICE_INFO
        printf("Platform: %s\n", platform_name);
        printDeviceInfo(device);
        std::cout << "CL_KERNEL_WORK_GROUP_SIZE: " << kernelMaxWorkGroupSize << std::endl;
        std::cout << "CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE: " << preferedWorkGroupSizeMultiple << std::endl;
#endif

        clReleaseKernel(kernel);
        buildProgram();
    }

    ParallelHeatSimulation::~ParallelHeatSimulation()
    {
        // clReleaseMemObject(bufInBoards);
        // clReleaseMemObject(bufOutStripMaxTs);
        // clReleaseMemObject(bufInnerForegoingTemperatures);
        // clReleaseMemObject(bufInnerNewTemperatures);
        // clReleaseKernel(kernel);
        clReleaseProgram(program);
        // clReleaseCommandQueue(queue);
        clReleaseContext(context);
    }

    std::vector<simulation_value_t> ParallelHeatSimulation::evaluateGeneration(const std::vector<cell_type_t> &fenotypes, simulation_value_t *minFinalTemperatures, simulation_steps_index_t *lastEquilibriumMoment)
    {
#ifdef BENCHMARK
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
#endif

        int board_size = boardHeight * boardWidth;
        int solutionsNumber = fenotypes.size() / board_size;

        if (solutionsNumber <= 0)
        {
            std::cout << "evaluateGeneration() received an empty list of solutions..." << std::endl;
            return std::vector<simulation_value_t>(0);
        }

        simulation_value_t *minTs;
        simulation_steps_index_t eqMoment;
        std::vector<simulation_value_t> maxTs = runSimulationKernel(fenotypes, solutionsNumber, &minTs, &eqMoment);

        if (minFinalTemperatures != NULL)
        {
            for (int cellIndex = 0; cellIndex < boardHeight * boardWidth; cellIndex++)
            {
                minFinalTemperatures[cellIndex] = minTs[cellIndex];
            }
        }

        if (lastEquilibriumMoment != NULL)
        {
            *lastEquilibriumMoment = eqMoment;
        }

        for (simulation_value_t &temperature : maxTs)
        {
            temperature *= -1;
        }
#ifdef BENCHMARK
        std::chrono::time_point<std::chrono::high_resolution_clock> stop = std::chrono::high_resolution_clock::now();
        std::cout << "BENCHMARK: generation evaluation took "
                  << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
                  << " microseconds" << std::endl;
#endif
        return maxTs;
    }

    std::vector<simulation_value_t> ParallelHeatSimulation::runSimulationKernel(
        const std::vector<cell_type_t> &boards, int individualsNumber, simulation_value_t **returnedMinFinalTemperatures, simulation_steps_index_t *returnedLastEquilibriumMoment)
    {
        timestamp start, stop;
        std::vector<simulation_value_t> maxTemperatures(individualsNumber);
        const int globalStripsNumber = individualsNumber * individualWidth * stripsPerColumn;
        const int stripsPerIndividual = individualWidth * stripsPerColumn;

        // TODO clCreateCommandQueue is deprecated???
        cl_int status;
        cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, NULL, &status);
        if (status != CL_SUCCESS)
        {
            std::cout << "clCreateCommandQueue failed. Simulation cancelled." << std::endl;
            throw std::runtime_error("Command queue error");
        }

        // ----------------------------------------------------
        // 5. Create Buffers
        // ----------------------------------------------------

#ifdef PARALLEL_SIMULATION_BENCHMARK
        start = std::chrono::high_resolution_clock::now();
#endif

        cl_mem bufInBoards = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                            boards.size() * sizeof(cell_type_t), (void *)boards.data(), &status);

        if (status != CL_SUCCESS)
        {
            printf("Boards input buffer allocation failed with status: %i\n", status);
            throw std::runtime_error("Buffer allocation error");
        }

        std::vector<simulation_value_t> startTs(boards.size());
        {
            int boardSize = boardHeight * boardWidth;
            for (int i = 0; i < individualsNumber; i++)
            {
                std::memcpy(startTs.data() + i * boardSize, startTemperatures, boardSize * sizeof(simulation_value_t));
            }
        }
        cl_mem bufInStartTemperatures = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                                       startTs.size() * sizeof(simulation_value_t), (void *)startTs.data(), &status);

        if (status != CL_SUCCESS)
        {
            printf("Start temperatures input buffer allocation failed with status: %i\n", status);
            throw std::runtime_error("Buffer allocation error");
        }

        std::vector<simulation_value_t> stripsMaxTs(globalStripsNumber);

        cl_mem bufOutStripMaxTs = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                                 stripsMaxTs.size() * sizeof(simulation_value_t), nullptr, &status);

        if (status != CL_SUCCESS)
        {
            printf("Max temperatures output buffer allocation failed with status: %i\n", status);
            throw std::runtime_error("Buffer allocation error");
        }

        std::vector<simulation_value_t> finalTs(boards.size());

        cl_mem bufOutFinalTs = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                              finalTs.size() * sizeof(simulation_value_t), nullptr, &status);

        if (status != CL_SUCCESS)
        {
            printf("Min temperatures output buffer allocation failed with status: %i\n", status);
            throw std::runtime_error("Buffer allocation error");
        }

        std::vector<simulation_steps_index_t> stripsEquilibriumMoments(globalStripsNumber);

        cl_mem bufOutStripEquilibriumMoments = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                                              stripsEquilibriumMoments.size() * sizeof(simulation_steps_index_t), nullptr, &status);

        if (status != CL_SUCCESS)
        {
            printf("Equilibrium moments output buffer allocation failed with status: %i\n", status);
            throw std::runtime_error("Buffer allocation error");
        }

        cl_mem bufInnerForegoingTemperatures = clCreateBuffer(context, CL_MEM_READ_WRITE,
                                                              boards.size() * sizeof(simulation_value_t), nullptr, &status);

        if (status != CL_SUCCESS)
        {
            printf("Foregoing temperatures buffer allocation failed with status: %i\n", status);
            throw std::runtime_error("Buffer allocation error");
        }

        cl_mem bufInnerNewTemperatures = clCreateBuffer(context, CL_MEM_READ_WRITE,
                                                        boards.size() * sizeof(simulation_value_t), nullptr, &status);

        if (status != CL_SUCCESS)
        {
            printf("New temperatures buffer allocation failed with status: %i\n", status);
            throw std::runtime_error("Buffer allocation error");
        }

#ifdef KERNEL_DEBUG
        std::vector<simulation_value_t> debug(globalStripsNumber);

        cl_mem bufOutDebug = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                            globalStripsNumber * sizeof(simulation_value_t), nullptr, &status);

        if (status != CL_SUCCESS)
        {
            std::cout << "Debug output buffer allocation failed with status: " << status << std::endl;
        }
#endif

#ifdef PARALLEL_SIMULATION_BENCHMARK
        stop = std::chrono::high_resolution_clock::now();
        std::cout << "Data copying took " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;
#endif

        cl_kernel kernel = clCreateKernel(program, "simulate_heat", &status);
        checkKernelCreationErrors(status);

        // ----------------------------------------------------
        // 7. Set Kernel Arguments
        // ----------------------------------------------------
        int arg_index = 0;
        clSetKernelArg(kernel, arg_index++, sizeof(cl_mem), &bufInBoards);
        clSetKernelArg(kernel, arg_index++, sizeof(cl_mem), &bufInStartTemperatures);
        clSetKernelArg(kernel, arg_index++, sizeof(cl_mem), &bufOutStripMaxTs);
        clSetKernelArg(kernel, arg_index++, sizeof(cl_mem), &bufOutFinalTs);
        clSetKernelArg(kernel, arg_index++, sizeof(cl_mem), &bufOutStripEquilibriumMoments);
        clSetKernelArg(kernel, arg_index++, sizeof(cl_mem), &bufInnerForegoingTemperatures);
        clSetKernelArg(kernel, arg_index++, sizeof(cl_mem), &bufInnerNewTemperatures);

#ifdef KERNEL_DEBUG
        clSetKernelArg(kernel, arg_index++, sizeof(cl_mem), &bufOutDebug);
        std::cout << "Debug kernel arg set" << std::endl;
#endif

        // ----------------------------------------------------
        // 8. Launch Kernel
        // ----------------------------------------------------

#ifdef PARALLEL_SIMULATION_BENCHMARK
        start = std::chrono::high_resolution_clock::now();
#endif

        size_t globalSize[2] = {(size_t)stripsPerColumn, (size_t)individualsNumber * individualWidth}; // in multiple dimensions (here a single one only)
        size_t work_group_size[2] = {(size_t)stripsPerColumn, (size_t)individualWidth};

        status = clEnqueueNDRangeKernel(queue, kernel, 2, nullptr,
                                        globalSize, work_group_size, 0, nullptr, nullptr);
        checkEnqueueNDRangeErrors(status);
        if (status != CL_SUCCESS)
        {
            std::cout << "EnqueueNDRange failed with status: " << status << std::endl;
        }

        clFinish(queue); // sychronizes

#ifdef PARALLEL_SIMULATION_BENCHMARK
        stop = std::chrono::high_resolution_clock::now();
        std::cout << "Kernel execution took " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
                  << std::endl;
#endif

        // ----------------------------------------------------
        // 9. Read Back Results
        // ----------------------------------------------------

#ifdef PARALLEL_SIMULATION_BENCHMARK
        start = std::chrono::high_resolution_clock::now();
#endif
        status = clEnqueueReadBuffer(queue, bufOutStripMaxTs, CL_TRUE, 0, stripsMaxTs.size() * sizeof(simulation_value_t),
                                     stripsMaxTs.data(),
                                     0, nullptr, nullptr);
        checkBufferReadErrors(status, "max temperatures");

        status = clEnqueueReadBuffer(queue, bufOutFinalTs, CL_TRUE, 0, finalTs.size() * sizeof(simulation_value_t),
                                     finalTs.data(),
                                     0, nullptr, nullptr);
        checkBufferReadErrors(status, "final temperatures");

        status = clEnqueueReadBuffer(queue, bufOutStripEquilibriumMoments, CL_TRUE, 0, stripsEquilibriumMoments.size() * sizeof(simulation_steps_index_t),
                                     stripsEquilibriumMoments.data(),
                                     0, nullptr, nullptr);
        checkBufferReadErrors(status, "equilibrium moments temperatures");

#ifdef PARALLEL_SIMULATION_BENCHMARK
        stop = std::chrono::high_resolution_clock::now();
        std::cout << "Reading kernel results took " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;
#endif

#ifdef KERNEL_DEBUG
        clEnqueueReadBuffer(queue, bufOutDebug, CL_TRUE, 0, debug.size() * sizeof(simulation_value_t), debug.data(),
                            0, nullptr, nullptr);

        for (int i = 0; i < individualsNumber; ++i)
        {
            std::cout << "DEBUG Individual: " << i << std::endl;
            for (int stripRow = 0; stripRow < stripsPerColumn; stripRow++)
            {
                for (int column = 0; column < individualWidth; column++)
                {
                    std::cout << debug[i * stripsPerIndividual + stripRow * individualWidth + column] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
#endif

        // ----------------------------------------------------
        // 10. Print output
        // ----------------------------------------------------

#ifdef DEBUG
        for (int i = 0; i < individualsNumber; ++i)
        {
            std::cout << "Max temperatures inside strips for individual: " << i << std::endl;

            for (int stripRow = 0; stripRow < stripsPerColumn; stripRow++)
            {
                for (int column = 0; column < individualWidth; column++)
                {
                    std::cout << debug[i * stripsPerIndividual + stripRow * individualWidth + column] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
#endif
#ifdef PARALLEL_SIMULATION_BENCHMARK
        start = std::chrono::high_resolution_clock::now();
#endif
        // aggregate obtained parameters
        simulation_value_t *minTemperatures = new simulation_value_t[boardHeight * boardWidth];
        std::memcpy(minTemperatures, finalTs.data(), boardHeight * boardWidth * sizeof(simulation_value_t));

        simulation_steps_index_t lastEquilibriumMoment = stripsEquilibriumMoments[0];
        for (int i = 0; i < individualsNumber; ++i)
        {
            int passedIndividualsStrips = i * stripsPerIndividual;
            simulation_value_t maxT = stripsMaxTs[passedIndividualsStrips];
            for (int stripIndex = 0; stripIndex < stripsPerIndividual; ++stripIndex)
            {
                simulation_value_t checkedT = stripsMaxTs[passedIndividualsStrips + stripIndex];
                maxT = std::max(maxT, checkedT);
                lastEquilibriumMoment = std::min(lastEquilibriumMoment, stripsEquilibriumMoments[passedIndividualsStrips + stripIndex]);
            }
            maxTemperatures[i] = maxT;

            int passedIndividuallsCells = i * boardHeight * boardWidth;
            for (int cellIndex = 0; cellIndex < boardHeight * boardWidth; cellIndex++)
            {
                minTemperatures[cellIndex] = std::min(minTemperatures[cellIndex], finalTs[passedIndividuallsCells + cellIndex]);
            }
        }
#ifdef PARALLEL_SIMULATION_BENCHMARK
        std::cout << "Minimal temperature in equilibrium state was: " << minTemperature << " K" << std::endl;
        std::cout << "Finding equilibrium took: " << lastEquilibriumMoment << " steps" << std::endl;
        stop = std::chrono::high_resolution_clock::now();
        std::cout << "Finding individuals maxima took " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;
#endif

#ifdef DEBUG
        std::cout << "Found maximum temperatures: ";
        for (double f : maxTemperatures)
            std::cout << f << " ";
        std::cout << std::endl;
#endif

        // ----------------------------------------------------
        // 11. Cleanup
        // ----------------------------------------------------

#ifdef PARALLEL_SIMULATION_BENCHMARK
        start = std::chrono::high_resolution_clock::now();
#endif
        clReleaseMemObject(bufInBoards);
        clReleaseMemObject(bufInStartTemperatures);
        clReleaseMemObject(bufOutStripMaxTs);
        clReleaseMemObject(bufOutFinalTs);
        clReleaseMemObject(bufOutStripEquilibriumMoments);
        clReleaseMemObject(bufInnerForegoingTemperatures);
        clReleaseMemObject(bufInnerNewTemperatures);
#ifdef KERNEL_DEBUG
        clReleaseMemObject(bufOutDebug);
#endif
        clReleaseKernel(kernel);
        clReleaseCommandQueue(queue);

#ifdef PARALLEL_SIMULATION_BENCHMARK
        stop = std::chrono::high_resolution_clock::now();
        std::cout << "Releasing kernel took " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count()
                  << std::endl;
#endif

        *returnedMinFinalTemperatures = minTemperatures;
        *returnedLastEquilibriumMoment = lastEquilibriumMoment;
        return maxTemperatures;
    }

    void ParallelHeatSimulation::setSimulationParams(const SimulationParams &params)
    {
        Simulation::setSimulationParams(params);
        buildProgram();
    }

    inline void ParallelHeatSimulation::buildProgram()
    {
        // std::string source = read_file("kernels/simulation_kernel_global_column_wise.cl");
        std::string source = loadKernel(
            kernels::simulation_kernel_global_stripped_column_wise,
            simulationSteps,
            ETA,
            GENERATOR_ALPHA,
            GENERATOR_BETA,
            CONDUCTOR_ALPHA,
            CONDUCTOR_BETA,
            DRAIN_ALPHA,
            delta_time,
            boardHeight, boardWidth, stripLength);
        const char *src = source.c_str();

        if (source.length() == 0)
        {
            std::cout << "Simulation kernel not loaded " << std::endl;
            throw std::runtime_error("Load error");
        }

        cl_int status;

        program = clCreateProgramWithSource(context, 1, &src, nullptr, &status);
        if (status != CL_SUCCESS)
        {
            std::cout << "Program creation from source failed with status: " << status << std::endl;
        }

        status = clBuildProgram(program, 1, &device, nullptr, nullptr, nullptr);

        if (status != CL_SUCCESS)
        {
            std::cerr << "Program build failed with status: " << status << std::endl;
            size_t logSize = 0;
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &logSize);
            std::vector<char> log(logSize);
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, logSize, log.data(), nullptr);
            std::cout << "Build error:\n"
                      << log.data() << "\n";
            throw std::runtime_error("Build error");
        }
    }

    inline int ParallelHeatSimulation::calculateStripLength(int maxWorkGroupSize)
    {
        std::cout << "DEBUG ParallelSim: H=" << individualHeight
                  << ", W=" << individualWidth
                  << ", MaxWG=" << maxWorkGroupSize << std::endl;

        for (int i = 1; i <= individualHeight; i++)
        {
            if (individualHeight % i == 0)
            {
                if ((individualHeight / i) * individualWidth <= maxWorkGroupSize)
                    return i;
            }
        }
        throw std::runtime_error("This line of code should be unreachable");
    }

}
