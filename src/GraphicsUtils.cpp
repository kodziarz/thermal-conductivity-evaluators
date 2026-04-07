#include "GraphicsUtils.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include <sstream>

namespace conductivity_evaluators
{

    void setKernelParam(std::string &kernel_code, std::string_view param_string, std::string_view replace_string)
    {
        size_t index = kernel_code.find(param_string);
        if (index == std::string::npos)
            throw std::runtime_error("Parameter '" + std::string(param_string) +
                                     "' not found in kernel source code.");

        kernel_code.replace(index, param_string.length(), replace_string);
    }

    std::string read_file(std::string_view path)
    {
        constexpr size_t read_size = 4096;
        std::ifstream stream(path.data());
        stream.exceptions(std::ios_base::badbit);

        if (!stream)
            throw std::ios_base::failure("file does not exist");

        std::string out, buf(read_size, '\0');

        while (stream.read(&buf[0], read_size))
            out.append(buf, 0, stream.gcount());

        out.append(buf, 0, stream.gcount());
        return out;
    }

    std::string loadKernel(std::string_view path,
                           int stepsNumber,
                           simulation_value_t ETA,
                           simulation_value_t generatorAlpha,
                           simulation_value_t generatorBeta,
                           simulation_value_t conductorAlpha,
                           simulation_value_t conductorBeta,
                           simulation_value_t drainAlpha,
                           simulation_value_t deltaTime,
                           int height,
                           int width,
                           int stripLength)
    {
        std::string kernel_code = read_file(path);

        setKernelParam(kernel_code, "#define STEPS_NUMBER 0", "#define STEPS_NUMBER " + std::to_string(stepsNumber));
        std::ostringstream sstream;
        sstream << ETA;
        setKernelParam(kernel_code, "#define ETA 0", "#define ETA " + sstream.str());
        setKernelParam(kernel_code, "#define GENERATOR_ALPHA 0", "#define GENERATOR_ALPHA " + std::to_string(generatorAlpha));
        setKernelParam(kernel_code, "#define GENERATOR_BETA 0", "#define GENERATOR_BETA " + std::to_string(generatorBeta));
        setKernelParam(kernel_code, "#define CONDUCTOR_ALPHA 0", "#define CONDUCTOR_ALPHA " + std::to_string(conductorAlpha));
        setKernelParam(kernel_code, "#define CONDUCTOR_BETA 0", "#define CONDUCTOR_BETA " + std::to_string(conductorBeta));
        setKernelParam(kernel_code, "#define DRAIN_ALPHA 0", "#define DRAIN_ALPHA " + std::to_string(drainAlpha));
        setKernelParam(kernel_code, "#define DELTA_TIME 0", "#define DELTA_TIME " + std::to_string(deltaTime));
        setKernelParam(kernel_code, "#define WIDTH 0", "#define WIDTH " + std::to_string(width));
        setKernelParam(kernel_code, "#define HEIGHT 0", "#define HEIGHT " + std::to_string(height));
        setKernelParam(kernel_code, "#define STRIP_LENGTH 1", "#define STRIP_LENGTH " + std::to_string(stripLength));

        return kernel_code;
    }

    cl_device_id getDevice(cl_platform_id platform, cl_device_type preferredDevice, cl_device_type fallbackDeviceType)
    {
        cl_device_id device;
        cl_uint status = clGetDeviceIDs(platform, preferredDevice, 1, &device, nullptr);

        if (status != CL_SUCCESS)
        {
            std::cout << "Preferred device not found. Trying the fallback one...\n";

            status = clGetDeviceIDs(platform, fallbackDeviceType, 1, &device, nullptr);

            if (status != CL_SUCCESS)
            {
                std::cout << "Fallback device also not found\n";
                throw std::runtime_error("Device not found");
            }
        }

        return device;
    }

    void printDeviceInfo(cl_device_id device)
    {
        char device_name[512];
        clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_name), device_name, nullptr);
        std::cout << "Device: " << device_name << "\n";

        cl_uint compute_units;
        clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(compute_units), &compute_units, nullptr);
        std::cout << "CL_DEVICE_MAX_COMPUTE_UNITS: " << compute_units << "\n";

        cl_ulong buffer_size;
        clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(buffer_size), &buffer_size, nullptr);
        std::cout << "CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE: " << buffer_size << "\n";

        clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, sizeof(buffer_size), &buffer_size, nullptr);
        std::cout << "CL_DEVICE_GLOBAL_MEM_CACHE_SIZE: " << buffer_size << "\n";

        size_t max_workgroup;
        clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_workgroup), &max_workgroup, nullptr);
        std::cout << "CL_DEVICE_MAX_WORK_GROUP_SIZE: " << max_workgroup << "\n";

        cl_uint dims;
        clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(dims), &dims, nullptr);
        std::vector<size_t> sizes(dims);
        clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t) * dims, sizes.data(), nullptr);

        std::cout << "Device work-group sizes: (";
        for (unsigned i = 0; i < dims; i++)
        {
            std::cout << sizes[i];
            if (i + 1 < dims)
                std::cout << ", ";
        }
        std::cout << ")\n";
    }

    //
    // Error checking helpers
    //

    void checkContextCreationErrors(cl_int status)
    {
        if (status != CL_SUCCESS)
        {
            std::cout << "clCreateContext failed. Error: " << status << "\n";
            throw std::runtime_error("Context creation error");
        }
    }

    void checkKernelCreationErrors(cl_int status)
    {
        if (status != CL_SUCCESS)
        {
            std::cout << "Kernel creation failed with status: " << status << "\n";
            throw std::runtime_error("Kernel creation error");
        }
    }

    void checkEnqueueNDRangeErrors(cl_int status)
    {
        if (status == CL_SUCCESS)
        {
            return;
        }

        std::cout << "clEnqueueNDRangeKernel failed. Kernel: "
                  << " Error: " << status << std::endl;

        switch (status)
        {
        case CL_INVALID_PROGRAM_EXECUTABLE:
            std::cout << "Error: CL_INVALID_PROGRAM_EXECUTABLE - No successfully built program executable is available for the device associated with the command queue." << std::endl;
            break;

        case CL_INVALID_COMMAND_QUEUE:
            std::cout << "Error: CL_INVALID_COMMAND_QUEUE - Command queue is not a valid host command queue." << std::endl;
            break;

        case CL_INVALID_KERNEL:
            std::cout << "Error: CL_INVALID_KERNEL - Kernel is not a valid kernel object." << std::endl;
            break;

        case CL_INVALID_CONTEXT:
            std::cout << "Error: CL_INVALID_CONTEXT - Context associated with the command queue, kernel, or event wait list does not match." << std::endl;
            break;

        case CL_INVALID_KERNEL_ARGS:
            std::cout << "Error: CL_INVALID_KERNEL_ARGS - Kernel argument values have not been specified correctly." << std::endl;
            break;

        case CL_INVALID_WORK_DIMENSION:
            std::cout << "Error: CL_INVALID_WORK_DIMENSION - work_dim is not between 1 and CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS." << std::endl;
            break;

        case CL_INVALID_GLOBAL_WORK_SIZE:
            std::cout << "Error: CL_INVALID_GLOBAL_WORK_SIZE - global_work_size is NULL, contains zero values, or exceeds device limits." << std::endl;
            break;

        case CL_INVALID_GLOBAL_OFFSET:
            std::cout << "Error: CL_INVALID_GLOBAL_OFFSET - global_work_offset is invalid, exceeds device limits, or is non-NULL before OpenCL 1.1." << std::endl;
            break;

        case CL_INVALID_WORK_GROUP_SIZE:
            std::cout << "Error: CL_INVALID_WORK_GROUP_SIZE - local_work_size is invalid, inconsistent with kernel requirements, or exceeds CL_KERNEL_WORK_GROUP_SIZE." << std::endl;
            break;

        case CL_INVALID_WORK_ITEM_SIZE:
            std::cout << "Error: CL_INVALID_WORK_ITEM_SIZE - One or more dimensions of local_work_size exceed CL_DEVICE_MAX_WORK_ITEM_SIZES." << std::endl;
            break;

        case CL_MISALIGNED_SUB_BUFFER_OFFSET:
            std::cout << "Error: CL_MISALIGNED_SUB_BUFFER_OFFSET - Sub-buffer offset is not aligned to CL_DEVICE_MEM_BASE_ADDR_ALIGN." << std::endl;
            break;

        case CL_INVALID_IMAGE_SIZE:
            std::cout << "Error: CL_INVALID_IMAGE_SIZE - Image dimensions or pitches are not supported by the device." << std::endl;
            break;

        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
            std::cout << "Error: CL_IMAGE_FORMAT_NOT_SUPPORTED - Image format is not supported by the device." << std::endl;
            break;

        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            std::cout << "Error: CL_MEM_OBJECT_ALLOCATION_FAILURE - Failed to allocate memory for buffer or image objects used by the kernel." << std::endl;
            break;

        case CL_INVALID_EVENT_WAIT_LIST:
            std::cout << "Error: CL_INVALID_EVENT_WAIT_LIST - event_wait_list or num_events_in_wait_list is invalid, or contains invalid events." << std::endl;
            break;

        case CL_INVALID_OPERATION:
            std::cout << "Error: CL_INVALID_OPERATION - Unsupported SVM or system pointers were passed as kernel arguments." << std::endl;
            break;

        case CL_OUT_OF_RESOURCES:
            std::cout << "Error: CL_OUT_OF_RESOURCES - Insufficient device resources to enqueue or execute the kernel." << std::endl;
            break;

        case CL_OUT_OF_HOST_MEMORY:
            std::cout << "Error: CL_OUT_OF_HOST_MEMORY - Not enough host memory to enqueue the kernel." << std::endl;
            break;

        default:
            std::cout << "Error: Unknown status code " << status << std::endl;
            break;
        }

        throw std::runtime_error("Kernel enqueue error occurred.");
    }

    void checkBufferReadErrors(cl_int status, std::string bufferName)
    {
        if (status == CL_SUCCESS)
        {
            return;
        }

        std::cout << "clEnqueueReadBuffer failed. Buffer: " << bufferName << "Error: " << status << std::endl;

        switch (status)
        {
        case CL_INVALID_COMMAND_QUEUE:
            std::cout << "Error: CL_INVALID_COMMAND_QUEUE - Command queue is not valid." << std::endl;
            break;

        case CL_INVALID_CONTEXT:
            std::cout << "Error: CL_INVALID_CONTEXT - The context associated with the command queue or buffer is invalid." << std::endl;
            break;

        case CL_INVALID_MEM_OBJECT:
            std::cout << "Error: CL_INVALID_MEM_OBJECT - The buffer is not a valid buffer object." << std::endl;
            break;

        case CL_INVALID_VALUE:
            std::cout << "Error: CL_INVALID_VALUE - The region specified by (offset, size) is out of bounds, or ptr is NULL." << std::endl;
            break;

        case CL_INVALID_EVENT_WAIT_LIST:
            std::cout << "Error: CL_INVALID_EVENT_WAIT_LIST - Invalid event wait list or mismatched event count." << std::endl;
            break;

        case CL_MISALIGNED_SUB_BUFFER_OFFSET:
            std::cout << "Error: CL_MISALIGNED_SUB_BUFFER_OFFSET - The sub-buffer offset is misaligned." << std::endl;
            break;

        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
            std::cout << "Error: CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST - Negative status in event wait list." << std::endl;
            break;

        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            std::cout << "Error: CL_MEM_OBJECT_ALLOCATION_FAILURE - Failure to allocate memory for buffer." << std::endl;
            break;

        case CL_INVALID_OPERATION:
            std::cout << "Error: CL_INVALID_OPERATION - Invalid operation for the buffer. It may be created with CL_MEM_HOST_WRITE_ONLY, CL_MEM_HOST_NO_ACCESS, CL_MEM_HOST_READ_ONLY, or CL_MEM_IMMUTABLE_EXT." << std::endl;
            break;

        case CL_OUT_OF_RESOURCES:
            std::cout << "Error: CL_OUT_OF_RESOURCES - Not enough resources to complete the operation on the device." << std::endl;
            break;

        case CL_OUT_OF_HOST_MEMORY:
            std::cout << "Error: CL_OUT_OF_HOST_MEMORY - Not enough memory on the host to perform the operation." << std::endl;
            break;

        default:
            std::cout << "Error: Unknown status code " << status << std::endl;
            break;
        }

        throw std::runtime_error("Buffer read error occurred.");
    }

}
