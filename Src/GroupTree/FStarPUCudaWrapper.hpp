// @SCALFMM_PRIVATE
#ifndef FSTARPUCUDAWRAPPER_HPP
#define FSTARPUCUDAWRAPPER_HPP

#include "../Utils/FGlobal.hpp"
#include "../Core/FCoreCommon.hpp"
#include "../Utils/FQuickSort.hpp"
#include "../Containers/FTreeCoordinate.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FAssert.hpp"
#include "../Utils/FAlignedMemory.hpp"
#include "../Utils/FAssert.hpp"

#include "FOutOfBlockInteraction.hpp"

#ifdef ScalFMM_USE_MPI
#include "../Utils/FMpi.hpp"
#endif

#include <vector>
#include <memory>

#include <omp.h>

#include <starpu.h>

#ifdef STARPU_USE_MPI
#include <starpu_mpi.h>
#endif

#include "Cuda/FCudaDeviceWrapper.hpp"

#include "FStarPUUtils.hpp"

template <class KernelClass, class CudaCellGroupClass,
          class CudaParticleGroupClass, class CudaParticleContainerClass,
          class CudaKernelClass>
class FStarPUCudaWrapper {
protected:
    typedef FStarPUCudaWrapper<KernelClass, CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass> ThisClass;

    template <class OtherBlockClass>
    struct BlockInteractions{
        OtherBlockClass* otherBlock;
        int otherBlockId;
        std::vector<OutOfBlockInteraction> interactions;
    };

    const int treeHeight;
    CudaKernelClass* kernels[STARPU_MAXCUDADEVS];        //< The kernels

public:
    FStarPUCudaWrapper(const int inTreeHeight): treeHeight(inTreeHeight){
        memset(kernels, 0, sizeof(CudaKernelClass*)*STARPU_MAXCUDADEVS);
    }

    void initKernel(const int workerId, KernelClass* originalKernel){
        FAssertLF(kernels[workerId] == nullptr);
        kernels[workerId] = FCuda__BuildCudaKernel<CudaKernelClass>(originalKernel);
    }

    ~FStarPUCudaWrapper(){
        for(int idxKernel = 0 ; idxKernel < STARPU_MAXCUDADEVS ; ++idxKernel ){
            if(kernels[idxKernel]){
                FCuda__ReleaseCudaKernel(kernels[idxKernel]);
            }
        }
    }

    static void bottomPassCallback(void *buffers[], void *cl_arg){
        //CudaCellGroupClass leafCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
        //                    STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        //CudaParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
        //                    STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));
        FStarPUPtrInterface* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__bottomPassCallback<CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]),
                kernel);
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Upward Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void upwardPassCallback(void *buffers[], void *cl_arg){
        //CudaCellGroupClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
        //                                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));

        FStarPUPtrInterface* worker = nullptr;
        int nbSubCellGroups = 0;
        int idxLevel = 0;
        starpu_codelet_unpack_args(cl_arg, &worker, &nbSubCellGroups, &idxLevel);

        unsigned char* subCellGroupsPtr[9] ;
        memset(subCellGroupsPtr, 0, 9*sizeof(unsigned char*));
        size_t subCellGroupsSize[9] ;
        memset(subCellGroupsPtr, 0, 9*sizeof(unsigned char*));
        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            subCellGroupsPtr[idxSubGroup] = ((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[idxSubGroup+1]));
            subCellGroupsSize[idxSubGroup] = STARPU_VARIABLE_GET_ELEMSIZE(buffers[idxSubGroup+1]);
        }

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__upwardPassCallback<CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                subCellGroupsPtr,subCellGroupsSize,
                kernel, nbSubCellGroups, idxLevel);
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Transfer Pass Mpi
    /////////////////////////////////////////////////////////////////////////////////////
#ifdef STARPU_USE_MPI
    static void transferInoutPassCallbackMpi(void *buffers[], void *cl_arg){
        // CudaCellGroupClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
        //                                 STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        // CudaCellGroupClass externalCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
        //                                STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        FStarPUPtrInterface* worker = nullptr;
        int idxLevel = 0;
        const std::vector<OutOfBlockInteraction>* outsideInteractions;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel, &outsideInteractions);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__transferInoutPassCallbackMpi<CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]),
                kernel, idxLevel, outsideInteractions->data(), outsideInteractions->size());
    }
#endif
    /////////////////////////////////////////////////////////////////////////////////////
    /// Transfer Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void transferInPassCallback(void *buffers[], void *cl_arg){
        FAssertLF(STARPU_VARIABLE_GET_PTR(buffers[0]) == STARPU_VARIABLE_GET_PTR(buffers[1]));
        //CudaCellGroupClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
        //                                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));

        FStarPUPtrInterface* worker = nullptr;
        int idxLevel = 0;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__transferInPassCallback<CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                kernel, idxLevel);
    }

    static void transferInoutPassCallback(void *buffers[], void *cl_arg){
        FAssertLF(STARPU_VARIABLE_GET_PTR(buffers[0]) == STARPU_VARIABLE_GET_PTR(buffers[2]));
        FAssertLF(STARPU_VARIABLE_GET_PTR(buffers[1]) == STARPU_VARIABLE_GET_PTR(buffers[3]));

        // CudaCellGroupClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
        //                                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        // CudaCellGroupClass externalCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
        //                                STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        FStarPUPtrInterface* worker = nullptr;
        int idxLevel = 0;
        const std::vector<OutOfBlockInteraction>* outsideInteractions;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel, &outsideInteractions);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__transferInoutPassCallback<CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]),
                kernel, idxLevel, outsideInteractions->data(), outsideInteractions->size());
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Downard Pass
    /////////////////////////////////////////////////////////////////////////////////////
    static void downardPassCallback(void *buffers[], void *cl_arg){
        //CudaCellGroupClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
        //                                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));

        FStarPUPtrInterface* worker = nullptr;
        int nbSubCellGroups = 0;
        int idxLevel = 0;
        starpu_codelet_unpack_args(cl_arg, &worker, &nbSubCellGroups, &idxLevel);

        unsigned char* subCellGroupsPtr[9];
        memset(subCellGroupsPtr, 0, 9*sizeof(unsigned char*));
        size_t subCellGroupsSize[9];
        memset(subCellGroupsPtr, 0, 9*sizeof(size_t));
        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            subCellGroupsPtr[idxSubGroup] = ((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[idxSubGroup+1]));
            subCellGroupsSize[idxSubGroup] = (STARPU_VARIABLE_GET_ELEMSIZE(buffers[idxSubGroup+1]));
        }

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__downardPassCallback<CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                subCellGroupsPtr,subCellGroupsSize,
                kernel, nbSubCellGroups, idxLevel);
    }
    /////////////////////////////////////////////////////////////////////////////////////
    /// Direct Pass MPI
    /////////////////////////////////////////////////////////////////////////////////////

#ifdef STARPU_USE_MPI
    static void directInoutPassCallbackMpi(void *buffers[], void *cl_arg){
        //CudaParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
        //                              STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        //CudaParticleGroupClass externalContainers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
        //                              STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        FStarPUPtrInterface* worker = nullptr;
        const std::vector<OutOfBlockInteraction>* outsideInteractions = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker, &outsideInteractions);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__directInoutPassCallbackMpi<CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]),
                kernel, outsideInteractions->data(), outsideInteractions->size(), worker->get<ThisClass>(FSTARPU_CPU_IDX)->treeHeight);
    }
#endif
    /////////////////////////////////////////////////////////////////////////////////////
    /// Direct Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void directInPassCallback(void *buffers[], void *cl_arg){
        // CudaParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
        //                              STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));

        FStarPUPtrInterface* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);
        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__directInPassCallback<CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                kernel, worker->get<ThisClass>(FSTARPU_CPU_IDX)->treeHeight);
    }

    static void directInoutPassCallback(void *buffers[], void *cl_arg){
        // CudaParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
        //                              STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        // CudaParticleGroupClass externalContainers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
        //                              STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        FStarPUPtrInterface* worker = nullptr;
        const std::vector<OutOfBlockInteraction>* outsideInteractions = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker, &outsideInteractions);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__directInoutPassCallback<CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]),
                kernel, outsideInteractions->data(), outsideInteractions->size(), worker->get<ThisClass>(FSTARPU_CPU_IDX)->treeHeight);
    }


    /////////////////////////////////////////////////////////////////////////////////////
    /// Merge Pass
    /////////////////////////////////////////////////////////////////////////////////////

    static void mergePassCallback(void *buffers[], void *cl_arg){
        // CudaCellGroupClass leafCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
        //                             STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        // CudaParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
        //                             STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        FStarPUPtrInterface* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);

        CudaKernelClass* kernel = worker->get<ThisClass>(FSTARPU_CPU_IDX)->kernels[starpu_worker_get_id()];

        FCuda__mergePassCallback<CudaCellGroupClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass>((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]),
                (unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]),
                kernel);
    }
};


#endif // FSTARPUCUDAWRAPPER_HPP

