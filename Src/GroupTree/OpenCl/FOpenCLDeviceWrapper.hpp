// @SCALFMM_PRIVATE
#ifndef FOPENCLDEVICEWRAPPER_HPP
#define FOPENCLDEVICEWRAPPER_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Core/FCoreCommon.hpp"
#include "../../Utils/FQuickSort.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../../Utils/FLog.hpp"
#include "../../Utils/FTic.hpp"
#include "../../Utils/FAssert.hpp"
#include "../../Utils/FAlignedMemory.hpp"
#include "../../Utils/FAssert.hpp"

#include "../FOutOfBlockInteraction.hpp"

#include <starpu.h>

struct FEmptyOpenCLFilename{
    operator const char*(){
        return nullptr;
    }
};

template <class OriginalKernelClass, class KernelFilenameClass = FEmptyOpenCLFilename>
class FOpenCLDeviceWrapper {
protected:
    static void SetKernelArgs(cl_kernel& kernel, const int pos){
    }
    template <class ParamClass, class... Args>
    static void SetKernelArgs(cl_kernel& kernel, const int pos, ParamClass* param, Args... args){
        FAssertLF(clSetKernelArg(kernel, pos, sizeof(*param), param) == 0);
        SetKernelArgs(kernel, pos+1, args...);
    }

    int workerId;
    int workerDevid;

    struct starpu_opencl_program opencl_code;

    cl_context context;

    cl_kernel kernel_bottomPassPerform;
    cl_command_queue queue_bottomPassPerform;

    cl_kernel kernel_upwardPassPerform;
    cl_command_queue queue_upwardPassPerform;

    cl_kernel kernel_transferInoutPassPerformMpi;
    cl_command_queue queue_transferInoutPassPerformMpi;

    cl_kernel kernel_transferInPassPerform;
    cl_command_queue queue_transferInPassPerform;

    cl_kernel kernel_transferInoutPassPerform;
    cl_command_queue queue_transferInoutPassPerform;

    cl_kernel kernel_downardPassPerform;
    cl_command_queue queue_downardPassPerform;

    cl_kernel kernel_directInoutPassPerformMpi;
    cl_command_queue queue_directInoutPassPerformMpi;

    cl_kernel kernel_directInoutPassPerform;
    cl_command_queue queue_directInoutPassPerform;

    cl_kernel kernel_directInPassPerform;
    cl_command_queue queue_directInPassPerform;

    cl_kernel kernel_mergePassPerform;
    cl_command_queue queue_mergePassPerform;

public:
    FOpenCLDeviceWrapper() : workerId(0) , workerDevid(0){
        workerId = starpu_worker_get_id();
        workerDevid = starpu_worker_get_devid(workerId);

        KernelFilenameClass kernelFilename;
        const char* filename = kernelFilename;
        if(filename){
            starpu_opencl_get_context (workerDevid, &context);

            const int err = starpu_opencl_load_opencl_from_file(filename, &opencl_code, NULL);
            if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

            FAssertLF( starpu_opencl_load_kernel(&kernel_bottomPassPerform, &queue_bottomPassPerform, &opencl_code, "bottomPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_upwardPassPerform, &queue_upwardPassPerform, &opencl_code, "upwardPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_transferInoutPassPerformMpi, &queue_transferInoutPassPerformMpi, &opencl_code, "transferInoutPassPerformMpi", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_transferInPassPerform, &queue_transferInPassPerform, &opencl_code, "transferInPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_transferInoutPassPerform, &queue_transferInoutPassPerform, &opencl_code, "transferInoutPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_downardPassPerform, &queue_downardPassPerform, &opencl_code, "downardPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_directInoutPassPerformMpi, &queue_directInoutPassPerformMpi, &opencl_code, "directInoutPassPerformMpi", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_directInoutPassPerform, &queue_directInoutPassPerform, &opencl_code, "directInoutPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_directInPassPerform, &queue_directInPassPerform, &opencl_code, "directInPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_mergePassPerform, &queue_mergePassPerform, &opencl_code, "mergePassPerform", workerDevid) == CL_SUCCESS);
        }
    }

    virtual void initDeviceFromKernel(const OriginalKernelClass& /*originalKernel*/){
    }

    virtual void releaseKernel(){
    }

    virtual ~FOpenCLDeviceWrapper(){
        // Release
        releaseKernel();
        KernelFilenameClass kernelFilename;
        const char* filename = kernelFilename;
        if(filename){
            const int err = starpu_opencl_unload_opencl(&opencl_code);
            if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
        }
    }

    cl_context& getOpenCLContext(){
        return context;
    }

    void bottomPassPerform(cl_mem leafCellsPtr,  size_t leafCellsSize, cl_mem containersPtr,  size_t containersSize){
        SetKernelArgs(kernel_bottomPassPerform, 0, &leafCellsPtr,  &leafCellsSize, &containersPtr,  &containersSize);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_bottomPassPerform, kernel_bottomPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }


    void upwardPassPerform(cl_mem currentCellsPtr,  size_t currentCellsSize, cl_mem subCellGroupsPtr[9],  size_t subCellGroupsSize[9], int nbSubCellGroups, int idxLevel){
        SetKernelArgs(kernel_upwardPassPerform, 0, &currentCellsPtr, &currentCellsSize, &subCellGroupsPtr,  &subCellGroupsSize, &nbSubCellGroups, &idxLevel);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_upwardPassPerform, kernel_upwardPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void transferInoutPassPerformMpi(cl_mem currentCellsPtr,
                         size_t currentCellsSize, cl_mem externalCellsPtr,  size_t externalCellsSize, int idxLevel, cl_mem outsideInteractionsCl,
                                                                                     size_t  outsideInteractionsSize){
        SetKernelArgs(kernel_transferInoutPassPerformMpi, 0, &currentCellsPtr,&currentCellsSize, &externalCellsPtr,  &externalCellsSize, &idxLevel, &outsideInteractionsCl,
                                                                                  &outsideInteractionsSize);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_transferInoutPassPerformMpi, kernel_transferInoutPassPerformMpi, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void transferInPassPerform(cl_mem currentCellsPtr, size_t currentCellsSize, int idxLevel){
        SetKernelArgs(kernel_transferInPassPerform, 0, &currentCellsPtr, &currentCellsSize, &idxLevel);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_transferInPassPerform, kernel_transferInPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void transferInoutPassPerform(cl_mem currentCellsPtr,
                          size_t currentCellsSize, cl_mem externalCellsPtr,  size_t externalCellsSize, int idxLevel, cl_mem outsideInteractionsCl,
                                                                                  size_t outsideInteractionsSize){
        SetKernelArgs(kernel_transferInoutPassPerform, 0, &currentCellsPtr,&currentCellsSize, &externalCellsPtr, &externalCellsSize, &idxLevel, &outsideInteractionsCl,&outsideInteractionsSize);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_transferInoutPassPerform, kernel_transferInoutPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void downardPassPerform(cl_mem currentCellsPtr,
                                     size_t currentCellsSize, cl_mem subCellGroupsPtr[9],  size_t subCellGroupsSize[9], int nbSubCellGroups, int idxLevel){
        SetKernelArgs(kernel_downardPassPerform, 0, &currentCellsPtr,
                      &currentCellsSize, &subCellGroupsPtr,  &subCellGroupsSize, &nbSubCellGroups, &idxLevel);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_downardPassPerform, kernel_downardPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void directInoutPassPerformMpi(cl_mem containersPtr,
                   size_t containersSize, cl_mem externalContainersPtr,  size_t externalContainersSize, cl_mem outsideInteractionsCl,
                                                                                size_t outsideInteractionsSize){
        SetKernelArgs(kernel_directInoutPassPerformMpi, 0, &containersPtr,
                      &containersSize, &externalContainersPtr, &externalContainersSize, &outsideInteractionsCl,&outsideInteractionsSize);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_directInoutPassPerformMpi, kernel_directInoutPassPerformMpi, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void directInPassPerform(cl_mem containersPtr,  size_t containerSize){
        SetKernelArgs(kernel_directInPassPerform, 0, &containersPtr, &containerSize);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_directInPassPerform, kernel_directInPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void directInoutPassPerform(cl_mem containersPtr,
                          size_t containerSize, cl_mem externalContainersPtr,  size_t externalContainersSize, cl_mem outsideInteractionsCl,
                          size_t  outsideInteractionsSize){
        SetKernelArgs(kernel_directInoutPassPerform, 0, &containersPtr,
                      &containerSize, &externalContainersPtr, &externalContainersSize, &outsideInteractionsCl, &outsideInteractionsSize);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_directInoutPassPerform, kernel_directInoutPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void mergePassPerform(cl_mem leafCellsPtr,
                            size_t leafCellsSize, cl_mem containersPtr, size_t containersSize){
        SetKernelArgs(kernel_mergePassPerform, 0, &leafCellsPtr, &leafCellsSize, &containersPtr, &containersSize);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_mergePassPerform, kernel_mergePassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }
};

#endif // FOPENCLDEVICEWRAPPER_HPP

