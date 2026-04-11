#ifndef GRAPHICSUTILS_H
#define GRAPHICSUTILS_H

#include <CL/cl.h>
#include <string>
#include <string_view>
#include <vector>

#include "Types.h"

namespace conductivity_evaluators
{

    // Function declarations
    void setKernelParam(std::string &kernel_code, std::string_view param_string, std::string_view replace_string);
    std::string read_file(std::string_view path);
    std::string loadKernel(std::string_view path,
                           int stepsNumber,
                           simulation_value_t RESULTANT_POWER_TOL,
                           simulation_value_t generatorAlpha,
                           simulation_value_t generatorBeta,
                           simulation_value_t conductorAlpha,
                           simulation_value_t conductorBeta,
                           simulation_value_t drainAlpha,
                           simulation_value_t deltaTime,
                           int height,
                           int width,
                           int thickness,
                           int stripLength);

    cl_device_id getDevice(cl_platform_id platform,
                           cl_device_type preferredDevice,
                           cl_device_type fallbackDeviceType);

    void printDeviceInfo(cl_device_id device);

    void checkContextCreationErrors(cl_int status);
    void checkKernelCreationErrors(cl_int status);
    void checkEnqueueNDRangeErrors(cl_int status);
    void checkBufferReadErrors(cl_int status, std::string bufferName);

}

#endif // GRAPHICSUTILS_H
