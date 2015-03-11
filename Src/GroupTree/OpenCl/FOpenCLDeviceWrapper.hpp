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
    struct Uptr9{
        cl_mem ptrs[9];
    };

    struct size_t9{
        size_t v[9];
    };

    static void SetKernelArgs(cl_kernel& kernel, const int pos){
    }
    template <class ParamClass, class... Args>
    static void SetKernelArgs(cl_kernel& kernel, const int pos, ParamClass* param, Args... args){
        FAssertLF(clSetKernelArg(kernel, pos, sizeof(*param), param) == 0,
                  "Error when assigning opencl argument ", pos);
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

    cl_mem user_data;

    int treeHeight;
public:
    FOpenCLDeviceWrapper(const int inTreeHeight) : workerId(0) , workerDevid(0), user_data(0), treeHeight(inTreeHeight){
        workerId = starpu_worker_get_id();
        workerDevid = starpu_worker_get_devid(workerId);

        KernelFilenameClass kernelFilename;
        const char* filename = kernelFilename;
        if(filename){
            starpu_opencl_get_context (workerDevid, &context);

            const int err = starpu_opencl_load_opencl_from_string(filename, &opencl_code, "-cl-std=CL2.0 -cl-mad-enable -Werror");
            if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

            FAssertLF( starpu_opencl_load_kernel(&kernel_bottomPassPerform, &queue_bottomPassPerform, &opencl_code, "FOpenCL__bottomPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_upwardPassPerform, &queue_upwardPassPerform, &opencl_code, "FOpenCL__upwardPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_transferInoutPassPerformMpi, &queue_transferInoutPassPerformMpi, &opencl_code, "FOpenCL__transferInoutPassPerformMpi", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_transferInPassPerform, &queue_transferInPassPerform, &opencl_code, "FOpenCL__transferInPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_transferInoutPassPerform, &queue_transferInoutPassPerform, &opencl_code, "FOpenCL__transferInoutPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_downardPassPerform, &queue_downardPassPerform, &opencl_code, "FOpenCL__downardPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_directInoutPassPerformMpi, &queue_directInoutPassPerformMpi, &opencl_code, "FOpenCL__directInoutPassPerformMpi", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_directInoutPassPerform, &queue_directInoutPassPerform, &opencl_code, "FOpenCL__directInoutPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_directInPassPerform, &queue_directInPassPerform, &opencl_code, "FOpenCL__directInPassPerform", workerDevid) == CL_SUCCESS);
            FAssertLF( starpu_opencl_load_kernel(&kernel_mergePassPerform, &queue_mergePassPerform, &opencl_code, "FOpenCL__mergePassPerform", workerDevid) == CL_SUCCESS);
        }
    }

    virtual void initDeviceFromKernel(const OriginalKernelClass& /*originalKernel*/){
    }

    virtual void releaseKernel(){
        int err;
        err = starpu_opencl_release_kernel(kernel_bottomPassPerform);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

        err = starpu_opencl_release_kernel(kernel_upwardPassPerform);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

        err = starpu_opencl_release_kernel(kernel_transferInoutPassPerformMpi);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

        err = starpu_opencl_release_kernel(kernel_transferInPassPerform);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

        err = starpu_opencl_release_kernel(kernel_transferInoutPassPerform);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

        err = starpu_opencl_release_kernel(kernel_downardPassPerform);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

        err = starpu_opencl_release_kernel(kernel_directInoutPassPerformMpi);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

        err = starpu_opencl_release_kernel(kernel_directInoutPassPerform);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

        err = starpu_opencl_release_kernel(kernel_directInPassPerform);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

        err = starpu_opencl_release_kernel(kernel_mergePassPerform);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);

        err = starpu_opencl_unload_opencl(&opencl_code);
        if(err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    virtual ~FOpenCLDeviceWrapper(){
    }

    cl_context& getOpenCLContext(){
        return context;
    }

    void bottomPassPerform(cl_mem leafCellsPtr,  size_t leafCellsSize, cl_mem containersPtr,  size_t containersSize){
        /*cl_int errcode_ret;
        const int size = sizeof(FTestCell);
        int* output = new int[size];
        cl_mem outputcl = clCreateBuffer(getOpenCLContext(),
           CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR,
           size*sizeof(int),
           (void*)output, &errcode_ret);
        FAssertLF(outputcl && errcode_ret == CL_SUCCESS, "OpenCL error code " , errcode_ret);*/

        SetKernelArgs(kernel_bottomPassPerform, 0, &leafCellsPtr,  &leafCellsSize, &containersPtr,  &containersSize, &user_data/*, &outputcl*/);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_bottomPassPerform, kernel_bottomPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
        /*errcode_ret = clEnqueueReadBuffer(queue_bottomPassPerform, outputcl,
                CL_TRUE, // blocking read
                0, // write from the start
                sizeof(int) * size,
                output, 0, NULL, NULL);
       FAssertLF(errcode_ret == CL_SUCCESS, "OpenCL error code " , errcode_ret);

        for(int idx = 0 ; idx < 10 ; ++idx){
            std::cout << "value " << idx << " = " << output[idx] << "\n";
        }
//       FTestCell* cell = (FTestCell*)output;
//       std::cout << " cell->getDataUp() " << cell->getDataUp() << "\n";

        clReleaseMemObject(outputcl);
        delete output;*/
    }


    void upwardPassPerform(cl_mem currentCellsPtr,  size_t currentCellsSize, cl_mem subCellGroupsPtr[9],  size_t subCellGroupsSize[9], int nbSubCellGroups, int idxLevel){
        return; // TODO
        Uptr9 ptrs;
        memcpy(ptrs.ptrs, subCellGroupsPtr, sizeof(cl_mem)*9);
        size_t9 sizes;
        memcpy(sizes.v, subCellGroupsSize, sizeof(size_t)*9);

        SetKernelArgs(kernel_upwardPassPerform, 0, &currentCellsPtr, &currentCellsSize, &ptrs,  &sizes, &nbSubCellGroups, &idxLevel, &user_data);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_upwardPassPerform, kernel_upwardPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void transferInoutPassPerformMpi(cl_mem currentCellsPtr,
                         size_t currentCellsSize, cl_mem externalCellsPtr,  size_t externalCellsSize, int idxLevel, cl_mem outsideInteractionsCl,
                                                                                     size_t  outsideInteractionsSize){
        return; // TODO
        SetKernelArgs(kernel_transferInoutPassPerformMpi, 0, &currentCellsPtr,&currentCellsSize, &externalCellsPtr,  &externalCellsSize, &idxLevel, &outsideInteractionsCl,
                                                                                  &outsideInteractionsSize, &user_data);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_transferInoutPassPerformMpi, kernel_transferInoutPassPerformMpi, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void transferInPassPerform(cl_mem currentCellsPtr, size_t currentCellsSize, int idxLevel){
        SetKernelArgs(kernel_transferInPassPerform, 0, &currentCellsPtr, &currentCellsSize, &idxLevel, &user_data);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_transferInPassPerform, kernel_transferInPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void transferInoutPassPerform(cl_mem currentCellsPtr,
                          size_t currentCellsSize, cl_mem externalCellsPtr,  size_t externalCellsSize, int idxLevel, cl_mem outsideInteractionsCl,
                                                                                  size_t outsideInteractionsSize){
        return; // TODO
        SetKernelArgs(kernel_transferInoutPassPerform, 0, &currentCellsPtr,&currentCellsSize, &externalCellsPtr, &externalCellsSize, &idxLevel,
                      &outsideInteractionsCl,&outsideInteractionsSize, &user_data);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_transferInoutPassPerform, kernel_transferInoutPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void downardPassPerform(cl_mem currentCellsPtr,
                                     size_t currentCellsSize, cl_mem subCellGroupsPtr[9],  size_t subCellGroupsSize[9], int nbSubCellGroups, int idxLevel){
        return; // TODO
        Uptr9 ptrs;
        memcpy(ptrs.ptrs, subCellGroupsPtr, sizeof(cl_mem)*9);
        size_t9 sizes;
        memcpy(sizes.v, subCellGroupsSize, sizeof(size_t)*9);

        SetKernelArgs(kernel_downardPassPerform, 0, &currentCellsPtr,
                      &currentCellsSize,  &ptrs, &sizes, &nbSubCellGroups, &idxLevel, &user_data);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_downardPassPerform, kernel_downardPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void directInoutPassPerformMpi(cl_mem containersPtr,
                   size_t containersSize, cl_mem externalContainersPtr,  size_t externalContainersSize, cl_mem outsideInteractionsCl,
                                                                                size_t outsideInteractionsSize){
        return; // TODO
        SetKernelArgs(kernel_directInoutPassPerformMpi, 0, &containersPtr,
                      &containersSize, &externalContainersPtr, &externalContainersSize, &outsideInteractionsCl,&outsideInteractionsSize, &treeHeight, &user_data);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_directInoutPassPerformMpi, kernel_directInoutPassPerformMpi, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void directInPassPerform(cl_mem containersPtr,  size_t containerSize){
        return; // TODO
        SetKernelArgs(kernel_directInPassPerform, 0, &containersPtr, &containerSize, &treeHeight, &user_data);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_directInPassPerform, kernel_directInPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void directInoutPassPerform(cl_mem containersPtr,
                          size_t containerSize, cl_mem externalContainersPtr,  size_t externalContainersSize, cl_mem outsideInteractionsCl,
                          size_t  outsideInteractionsSize){
        return; // TODO
        SetKernelArgs(kernel_directInoutPassPerform, 0, &containersPtr,
                      &containerSize, &externalContainersPtr, &externalContainersSize, &outsideInteractionsCl, &outsideInteractionsSize, &treeHeight, &user_data);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_directInoutPassPerform, kernel_directInoutPassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }

    void mergePassPerform(cl_mem leafCellsPtr,
                            size_t leafCellsSize, cl_mem containersPtr, size_t containersSize){
        return; // TODO
        SetKernelArgs(kernel_mergePassPerform, 0, &leafCellsPtr, &leafCellsSize, &containersPtr, &containersSize, &user_data);
        size_t dim = 1;
        const int err = clEnqueueNDRangeKernel(queue_mergePassPerform, kernel_mergePassPerform, 1, NULL, &dim, NULL, 0, NULL, NULL);
        if (err != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(err);
    }
};

#endif // FOPENCLDEVICEWRAPPER_HPP

