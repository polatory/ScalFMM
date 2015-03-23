// Keep in private GIT
// @SCALFMM_PRIVATE
#ifndef FGROUPTASKSTARPUMPIALGORITHM_HPP
#define FGROUPTASKSTARPUMPIALGORITHM_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Core/FCoreCommon.hpp"
#include "../../Utils/FQuickSort.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../../Utils/FLog.hpp"
#include "../../Utils/FTic.hpp"
#include "../../Utils/FAssert.hpp"
#include "../../Utils/FAlignedMemory.hpp"
#include "../../Utils/FAssert.hpp"

#include "../../Utils/FMpi.hpp"

#include "FOutOfBlockInteraction.hpp"

#include <vector>
#include <memory>

#include <omp.h>

#include <starpu.h>
#include <starpu_mpi.h>
#include "../StarPUUtils/FStarPUUtils.hpp"

#ifdef STARPU_USE_CPU
#include "../StarPUUtils/FStarPUCpuWrapper.hpp"
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
#include "../StarPUUtils/FStarPUCudaWrapper.hpp"
#include "../Cuda/FCudaEmptyKernel.hpp"
#include "../Cuda/FCudaGroupAttachedLeaf.hpp"
#include "../Cuda/FCudaGroupOfParticles.hpp"
#include "../Cuda/FCudaGroupOfCells.hpp"
#include "../Cuda/FCudaEmptyCellSymb.hpp"
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
#include "../StarPUUtils/FStarPUOpenClWrapper.hpp"
#include "../OpenCl/FOpenCLDeviceWrapper.hpp"
#include "../OpenCl/FEmptyOpenCLCode.hpp"
#endif


template <class OctreeClass, class CellContainerClass, class KernelClass, class ParticleGroupClass, class StarPUCpuWrapperClass
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
    , class StarPUCudaWrapperClass = FStarPUCudaWrapper<KernelClass, FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                                        FCudaGroupOfParticles<0, 0, int>, FCudaGroupAttachedLeaf<0, 0, int>, FCudaEmptyKernel>
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
    , class StarPUOpenClWrapperClass = FStarPUOpenClWrapper<KernelClass, FOpenCLDeviceWrapper<KernelClass>>
#endif
          >
class FGroupTaskStarPUMpiAlgorithm {
protected:
    typedef FGroupTaskStarPUMpiAlgorithm<OctreeClass, CellContainerClass, KernelClass, ParticleGroupClass, StarPUCpuWrapperClass
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        , StarPUCudaWrapperClass
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
    , StarPUOpenClWrapperClass
#endif
    > ThisClass;

    int getTag(const int inLevel, const MortonIndex mindex, const int mode) const{
        int shift = 0;
        int height = tree->getHeight();
        while(height) { shift += 1; height >>= 1; }
        return (((mindex<<shift) + inLevel) << 5) + mode;
    }

    const FMpi::FComm& comm;

    template <class OtherBlockClass>
    struct BlockInteractions{
        OtherBlockClass* otherBlock;
        int otherBlockId;
        std::vector<OutOfBlockInteraction> interactions;
    };

    struct CellHandles{
        starpu_data_handle_t symb;
        starpu_data_handle_t up;
        starpu_data_handle_t down;
    };

    struct ParticleHandles{
        starpu_data_handle_t symb;
        starpu_data_handle_t down;
    };

    std::vector< std::vector< std::vector<BlockInteractions<CellContainerClass>>>> externalInteractionsAllLevel;
    std::vector< std::vector<BlockInteractions<ParticleGroupClass>>> externalInteractionsLeafLevel;

    int MaxThreads;         //< The number of threads
    OctreeClass*const tree;       //< The Tree
    KernelClass*const originalCpuKernel;

    std::vector<CellHandles>* cellHandles;
    std::vector<ParticleHandles> particleHandles;

    starpu_codelet p2m_cl;
    starpu_codelet m2m_cl[9];
    starpu_codelet l2l_cl[9];
    starpu_codelet l2p_cl;

    starpu_codelet m2l_cl_in;
    starpu_codelet m2l_cl_inout;
    starpu_codelet m2l_cl_inout_mpi;

    starpu_codelet p2p_cl_in;
    starpu_codelet p2p_cl_inout;
    starpu_codelet p2p_cl_inout_mpi;

#ifdef STARPU_USE_CPU
    StarPUCpuWrapperClass cpuWrapper;
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
    StarPUCudaWrapperClass cudaWrapper;
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
    StarPUOpenClWrapperClass openclWrapper;
#endif

    FStarPUPtrInterface wrappers;
    FStarPUPtrInterface* wrapperptr;

public:
    FGroupTaskStarPUMpiAlgorithm(const FMpi::FComm& inComm, OctreeClass*const inTree, KernelClass* inKernels, const int inMaxThreads = -1)
        :   comm(inComm), MaxThreads(inMaxThreads), tree(inTree), originalCpuKernel(inKernels),
          cellHandles(nullptr),
#ifdef STARPU_USE_CPU
            cpuWrapper(tree->getHeight()),
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
            cudaWrapper(tree->getHeight()),
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
            openclWrapper(tree->getHeight()),
#endif
            wrapperptr(&wrappers){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");
        FAssertLF(MaxThreads <= STARPU_MAXCPUS, "number of threads to high");

        struct starpu_conf conf;
        FAssertLF(starpu_conf_init(&conf) == 0);
        //conf.ncpus = MaxThreads;
        FAssertLF(starpu_init(&conf) == 0);
        FAssertLF(starpu_mpi_init ( 0, 0, 0 ) == 0);

        starpu_pthread_mutex_t initMutex;
        starpu_pthread_mutex_init(&initMutex, NULL);
#ifdef STARPU_USE_CPU
        FStarPUUtils::ExecOnWorkers(STARPU_CPU, [&](){
            starpu_pthread_mutex_lock(&initMutex);
            cpuWrapper.initKernel(starpu_worker_get_id(), inKernels);
            starpu_pthread_mutex_unlock(&initMutex);
        });
        wrappers.set(FSTARPU_CPU_IDX, &cpuWrapper);
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        FStarPUUtils::ExecOnWorkers(STARPU_CUDA, [&](){
            starpu_pthread_mutex_lock(&initMutex);
            cudaWrapper.initKernel(starpu_worker_get_id(), inKernels);
            starpu_pthread_mutex_unlock(&initMutex);
        });
        wrappers.set(FSTARPU_CUDA_IDX, &cudaWrapper);
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        FStarPUUtils::ExecOnWorkers(STARPU_OPENCL, [&](){
            starpu_pthread_mutex_lock(&initMutex);
            openclWrapper.initKernel(starpu_worker_get_id(), inKernels);
            starpu_pthread_mutex_unlock(&initMutex);
        });
        wrappers.set(FSTARPU_OPENCL_IDX, &openclWrapper);
#endif
        starpu_pthread_mutex_destroy(&initMutex);

        starpu_pause();

        MaxThreads = starpu_worker_get_count();//starpu_cpu_worker_get_count();

        cellHandles   = new std::vector<CellHandles>[tree->getHeight()];

        initCodelet();
        initCodeletMpi();

        FLOG(FLog::Controller << "FGroupTaskStarPUAlgorithm (Max Worker " << starpu_worker_get_count() << ")\n");
#ifdef STARPU_USE_CPU
        FLOG(FLog::Controller << "FGroupTaskStarPUAlgorithm (Max CPU " << starpu_cpu_worker_get_count() << ")\n");
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        FLOG(FLog::Controller << "FGroupTaskStarPUAlgorithm (Max OpenCL " << starpu_opencl_worker_get_count() << ")\n");
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        FLOG(FLog::Controller << "FGroupTaskStarPUAlgorithm (Max CUDA " << starpu_cuda_worker_get_count() << ")\n");
#endif
    }

    ~FGroupTaskStarPUMpiAlgorithm(){
        starpu_resume();

        cleanHandle();
        cleanHandleMpi();
        delete[] cellHandles;

        starpu_pthread_mutex_t releaseMutex;
        starpu_pthread_mutex_init(&releaseMutex, NULL);
#ifdef STARPU_USE_CPU
        FStarPUUtils::ExecOnWorkers(STARPU_CPU, [&](){
            starpu_pthread_mutex_lock(&releaseMutex);
            cpuWrapper.releaseKernel(starpu_worker_get_id());
            starpu_pthread_mutex_unlock(&releaseMutex);
        });
        wrappers.set(FSTARPU_CPU_IDX, &cpuWrapper);
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        FStarPUUtils::ExecOnWorkers(STARPU_CUDA, [&](){
            starpu_pthread_mutex_lock(&releaseMutex);
            cudaWrapper.releaseKernel(starpu_worker_get_id());
            starpu_pthread_mutex_unlock(&releaseMutex);
        });
        wrappers.set(FSTARPU_CUDA_IDX, &cudaWrapper);
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        FStarPUUtils::ExecOnWorkers(STARPU_OPENCL, [&](){
            starpu_pthread_mutex_lock(&releaseMutex);
            openclWrapper.releaseKernel(starpu_worker_get_id());
            starpu_pthread_mutex_unlock(&releaseMutex);
        });
        wrappers.set(FSTARPU_OPENCL_IDX, &openclWrapper);
#endif
        starpu_pthread_mutex_destroy(&releaseMutex);

        starpu_mpi_shutdown();
        starpu_shutdown();
    }

    void execute(const unsigned operationsToProceed = FFmmNearAndFarFields){
        FLOG( FLog::Controller << "\tStart FGroupTaskStarPUMpiAlgorithm\n" );

        #pragma omp parallel
        #pragma omp single
        buildExternalInteractionVecs();
        buildHandles();

        #pragma omp parallel
        #pragma omp single
        buildRemoteInteractionsAndHandles();

        starpu_resume();
        postRecvAllocatedBlocks();

        if( operationsToProceed & FFmmP2P ) insertParticlesSend();

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();
        if(operationsToProceed & FFmmM2L) insertCellsSend();

        if(operationsToProceed & FFmmM2L) transferPass();
        if(operationsToProceed & FFmmM2L) transferPassMpi();

        if(operationsToProceed & FFmmL2L) downardPass();

        if( operationsToProceed & FFmmP2P ) directPass();
        if( operationsToProceed & FFmmP2P ) directPassMpi();

        if( operationsToProceed & FFmmL2P ) mergePass();

        starpu_task_wait_for_all();
        starpu_pause();
    }

protected:
    void initCodelet(){
        memset(&p2m_cl, 0, sizeof(p2m_cl));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportP2M(FSTARPU_CPU_IDX)){
            p2m_cl.cpu_funcs[0] = StarPUCpuWrapperClass::bottomPassCallback;
            p2m_cl.where |= STARPU_CPU;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportP2M(FSTARPU_CUDA_IDX)){
            p2m_cl.cuda_funcs[0] = StarPUCudaWrapperClass::bottomPassCallback;
            p2m_cl.where |= STARPU_CUDA;
        }
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        if(originalCpuKernel->supportP2M(FSTARPU_OPENCL_IDX)){
            p2m_cl.opencl_funcs[0] = StarPUOpenClWrapperClass::bottomPassCallback;
            p2m_cl.where |= STARPU_OPENCL;
        }
#endif
        p2m_cl.nbuffers = 3;
        p2m_cl.modes[0] = STARPU_R;
        p2m_cl.modes[1] = STARPU_RW;
        p2m_cl.modes[2] = STARPU_R;
        p2m_cl.name = "p2m_cl";

        memset(m2m_cl, 0, sizeof(m2m_cl[0])*9);
        memset(l2l_cl, 0, sizeof(l2l_cl[0])*9);
        for(int idx = 0 ; idx < 9 ; ++idx){
#ifdef STARPU_USE_CPU
            if(originalCpuKernel->supportM2M(FSTARPU_CPU_IDX)){
                m2m_cl[idx].cpu_funcs[0] = StarPUCpuWrapperClass::upwardPassCallback;
                m2m_cl[idx].where |= STARPU_CPU;
            }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
            if(originalCpuKernel->supportM2M(FSTARPU_CUDA_IDX)){
                m2m_cl[idx].cuda_funcs[0] = StarPUCudaWrapperClass::upwardPassCallback;
                m2m_cl[idx].where |= STARPU_CUDA;
            }
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
            if(originalCpuKernel->supportM2M(FSTARPU_OPENCL_IDX)){
                m2m_cl[idx].opencl_funcs[0] = StarPUOpenClWrapperClass::upwardPassCallback;
                m2m_cl[idx].where |= STARPU_OPENCL;
            }
#endif
            m2m_cl[idx].nbuffers = (idx+2)*2;
            m2m_cl[idx].dyn_modes = (starpu_data_access_mode*)malloc(m2m_cl[idx].nbuffers*sizeof(starpu_data_access_mode));
            m2m_cl[idx].dyn_modes[0] = STARPU_R;
            m2m_cl[idx].dyn_modes[1] = STARPU_RW;
            m2m_cl[idx].name = "m2m_cl";

#ifdef STARPU_USE_CPU
            if(originalCpuKernel->supportL2L(FSTARPU_CPU_IDX)){
                l2l_cl[idx].cpu_funcs[0] = StarPUCpuWrapperClass::downardPassCallback;
                l2l_cl[idx].where |= STARPU_CPU;
            }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
            if(originalCpuKernel->supportL2L(FSTARPU_CUDA_IDX)){
                l2l_cl[idx].cuda_funcs[0] = StarPUCudaWrapperClass::downardPassCallback;
                l2l_cl[idx].where |= STARPU_CUDA;
            }
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
            if(originalCpuKernel->supportL2L(FSTARPU_OPENCL_IDX)){
                l2l_cl[idx].opencl_funcs[0] = StarPUOpenClWrapperClass::downardPassCallback;
                l2l_cl[idx].where |= STARPU_OPENCL;
            }
#endif
            l2l_cl[idx].nbuffers = (idx+2)*2;
            l2l_cl[idx].dyn_modes = (starpu_data_access_mode*)malloc(l2l_cl[idx].nbuffers*sizeof(starpu_data_access_mode));
            l2l_cl[idx].dyn_modes[0] = STARPU_R;
            l2l_cl[idx].dyn_modes[1] = STARPU_R;
            l2l_cl[idx].name = "l2l_cl";

            for(int idxBuffer = 0 ; idxBuffer <= idx ; ++idxBuffer){
                m2m_cl[idx].dyn_modes[(idxBuffer*2)+2] = STARPU_R;
                m2m_cl[idx].dyn_modes[(idxBuffer*2)+3] = STARPU_R;
                l2l_cl[idx].dyn_modes[(idxBuffer*2)+2] = STARPU_R;
                l2l_cl[idx].dyn_modes[(idxBuffer*2)+3] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
            }
        }

        memset(&l2p_cl, 0, sizeof(l2p_cl));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportL2P(FSTARPU_CPU_IDX)){
            l2p_cl.cpu_funcs[0] = StarPUCpuWrapperClass::mergePassCallback;
            l2p_cl.where |= STARPU_CPU;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportL2P(FSTARPU_CUDA_IDX)){
            l2p_cl.cuda_funcs[0] = StarPUCudaWrapperClass::mergePassCallback;
            l2p_cl.where |= STARPU_CUDA;
        }
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        if(originalCpuKernel->supportL2P(FSTARPU_OPENCL_IDX)){
            l2p_cl.opencl_funcs[0] = StarPUOpenClWrapperClass::mergePassCallback;
            l2p_cl.where |= STARPU_OPENCL;
        }
#endif
        l2p_cl.nbuffers = 4;
        l2p_cl.modes[0] = STARPU_R;
        l2p_cl.modes[1] = STARPU_R;
        l2p_cl.modes[2] = STARPU_R;
        l2p_cl.modes[3] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        l2p_cl.name = "l2p_cl";

        memset(&p2p_cl_in, 0, sizeof(p2p_cl_in));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportP2P(FSTARPU_CPU_IDX)){
            p2p_cl_in.cpu_funcs[0] = StarPUCpuWrapperClass::directInPassCallback;
            p2p_cl_in.where |= STARPU_CPU;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportP2P(FSTARPU_CUDA_IDX)){
            p2p_cl_in.cuda_funcs[0] = StarPUCudaWrapperClass::directInPassCallback;
            p2p_cl_in.where |= STARPU_CUDA;
        }
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        if(originalCpuKernel->supportP2P(FSTARPU_OPENCL_IDX)){
            p2p_cl_in.opencl_funcs[0] = StarPUOpenClWrapperClass::directInPassCallback;
            p2p_cl_in.where |= STARPU_OPENCL;
        }
#endif
        p2p_cl_in.nbuffers = 2;
        p2p_cl_in.modes[0] = STARPU_R;
        p2p_cl_in.modes[1] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        p2p_cl_in.name = "p2p_cl_in";
        memset(&p2p_cl_inout, 0, sizeof(p2p_cl_inout));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportP2PExtern(FSTARPU_CPU_IDX)){
            p2p_cl_inout.cpu_funcs[0] = StarPUCpuWrapperClass::directInoutPassCallback;
            p2p_cl_inout.where |= STARPU_CPU;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportP2PExtern(FSTARPU_CUDA_IDX)){
            p2p_cl_inout.cuda_funcs[0] = StarPUCudaWrapperClass::directInoutPassCallback;
            p2p_cl_inout.where |= STARPU_CUDA;
        }
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        if(originalCpuKernel->supportP2PExtern(FSTARPU_OPENCL_IDX)){
            p2p_cl_inout.opencl_funcs[0] = StarPUOpenClWrapperClass::directInoutPassCallback;
            p2p_cl_inout.where |= STARPU_OPENCL;
        }
#endif
        p2p_cl_inout.nbuffers = 4;
        p2p_cl_inout.modes[0] = STARPU_R;
        p2p_cl_inout.modes[1] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        p2p_cl_inout.modes[2] = STARPU_R;
        p2p_cl_inout.modes[3] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        p2p_cl_inout.name = "p2p_cl_inout";

        memset(&m2l_cl_in, 0, sizeof(m2l_cl_in));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportM2L(FSTARPU_CPU_IDX)){
            m2l_cl_in.cpu_funcs[0] = StarPUCpuWrapperClass::transferInPassCallback;
            m2l_cl_in.where |= STARPU_CPU;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportM2L(FSTARPU_CUDA_IDX)){
            m2l_cl_in.cuda_funcs[0] = StarPUCudaWrapperClass::transferInPassCallback;
            m2l_cl_in.where |= STARPU_CUDA;
        }
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        if(originalCpuKernel->supportM2L(FSTARPU_OPENCL_IDX)){
            m2l_cl_in.opencl_funcs[0] = StarPUOpenClWrapperClass::transferInPassCallback;
            m2l_cl_in.where |= STARPU_OPENCL;
        }
#endif
        m2l_cl_in.nbuffers = 3;
        m2l_cl_in.modes[0] = STARPU_R;
        m2l_cl_in.modes[1] = STARPU_R;
        m2l_cl_in.modes[2] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        m2l_cl_in.name = "m2l_cl_in";
        memset(&m2l_cl_inout, 0, sizeof(m2l_cl_inout));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportM2LExtern(FSTARPU_CPU_IDX)){
            m2l_cl_inout.cpu_funcs[0] = StarPUCpuWrapperClass::transferInoutPassCallback;
            m2l_cl_inout.where |= STARPU_CPU;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportM2LExtern(FSTARPU_CUDA_IDX)){
            m2l_cl_inout.cuda_funcs[0] = StarPUCudaWrapperClass::transferInoutPassCallback;
            m2l_cl_inout.where |= STARPU_CUDA;
        }
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        if(originalCpuKernel->supportM2LExtern(FSTARPU_OPENCL_IDX)){
            m2l_cl_inout.opencl_funcs[0] = StarPUOpenClWrapperClass::transferInoutPassCallback;
            m2l_cl_inout.where |= STARPU_OPENCL;
        }
#endif
        m2l_cl_inout.nbuffers = 6;
        m2l_cl_inout.modes[0] = STARPU_R;
        m2l_cl_inout.modes[1] = STARPU_R;
        m2l_cl_inout.modes[2] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        m2l_cl_inout.modes[3] = STARPU_R;
        m2l_cl_inout.modes[4] = STARPU_R;
        m2l_cl_inout.modes[5] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        m2l_cl_inout.name = "m2l_cl_inout";
    }

    /** dealloc in a starpu way all the defined handles */
    void cleanHandle(){
        for(int idxLevel = 0 ; idxLevel < tree->getHeight() ; ++idxLevel){
            for(int idxHandle = 0 ; idxHandle < int(cellHandles[idxLevel].size()) ; ++idxHandle){
                starpu_data_unregister(cellHandles[idxLevel][idxHandle].symb);
                starpu_data_unregister(cellHandles[idxLevel][idxHandle].up);
                starpu_data_unregister(cellHandles[idxLevel][idxHandle].down);
            }
            cellHandles[idxLevel].clear();
        }
        {
            for(int idxHandle = 0 ; idxHandle < int(particleHandles.size()) ; ++idxHandle){
                starpu_data_unregister(particleHandles[idxHandle].symb);
                starpu_data_unregister(particleHandles[idxHandle].down);
            }
            particleHandles.clear();
        }
    }

    ////////////////////////////////////////////////////////////////////////////

    void initCodeletMpi(){
        memset(&p2p_cl_inout_mpi, 0, sizeof(p2p_cl_inout_mpi));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportP2PMpi(FSTARPU_CPU_IDX)){
            p2p_cl_inout_mpi.where |= STARPU_CPU;
            p2p_cl_inout_mpi.cpu_funcs[0] = StarPUCpuWrapperClass::directInoutPassCallbackMpi;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportP2PMpi(FSTARPU_CUDA_IDX)){
            p2p_cl_inout_mpi.where |= STARPU_CUDA;
            p2p_cl_inout_mpi.cuda_funcs[0] = StarPUCudaWrapperClass::directInoutPassCallbackMpi;
        }
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        if(originalCpuKernel->supportP2PMpi(FSTARPU_OPENCL_IDX)){
            p2p_cl_inout_mpi.where |= STARPU_OPENCL;
            p2p_cl_inout_mpi.opencl_funcs[0] = StarPUOpenClWrapperClass::directInoutPassCallbackMpi;
        }
#endif
        p2p_cl_inout_mpi.nbuffers = 3;
        p2p_cl_inout_mpi.modes[0] = STARPU_R;
        p2p_cl_inout_mpi.modes[1] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        p2p_cl_inout_mpi.modes[2] = STARPU_R;
        p2p_cl_inout_mpi.name = "p2p_cl_inout_mpi";

        memset(&m2l_cl_inout_mpi, 0, sizeof(m2l_cl_inout_mpi));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportM2LMpi(FSTARPU_CPU_IDX)){
            m2l_cl_inout_mpi.where |= STARPU_CPU;
            m2l_cl_inout_mpi.cpu_funcs[0] = StarPUCpuWrapperClass::transferInoutPassCallbackMpi;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportM2LMpi(FSTARPU_CUDA_IDX)){
            m2l_cl_inout_mpi.where |= STARPU_CUDA;
            m2l_cl_inout_mpi.cuda_funcs[0] = StarPUCudaWrapperClass::transferInoutPassCallbackMpi;
        }
#endif
#ifdef ScalFMM_ENABLE_OPENCL_KERNEL
        if(originalCpuKernel->supportM2LMpi(FSTARPU_OPENCL_IDX)){
            m2l_cl_inout_mpi.where |= STARPU_OPENCL;
            m2l_cl_inout_mpi.opencl_funcs[0] = StarPUOpenClWrapperClass::transferInoutPassCallbackMpi;
        }
#endif
        m2l_cl_inout_mpi.nbuffers = 4;
        m2l_cl_inout_mpi.modes[0] = STARPU_R;
        m2l_cl_inout_mpi.modes[1] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        m2l_cl_inout_mpi.modes[2] = STARPU_R;
        m2l_cl_inout_mpi.modes[3] = STARPU_R;
        m2l_cl_inout_mpi.name = "m2l_cl_inout_mpi";
    }

    std::vector<std::pair<MortonIndex,MortonIndex>> processesIntervalPerLevels;
    struct BlockDescriptor{
        MortonIndex firstIndex;
        MortonIndex lastIndex;
        int owner;
        int bufferSize;
        size_t bufferSizeSymb;
        size_t bufferSizeUp;
        size_t bufferSizeDown;
        size_t leavesBufferSize;
    };
    std::vector<std::vector<BlockDescriptor>> processesBlockInfos;
    std::vector<int> nbBlocksPerLevelAll;
    std::vector<int> nbBlocksBeforeMinPerLevel;

    std::vector< std::vector< std::vector<BlockInteractions<CellContainerClass>>>> externalInteractionsAllLevelMpi;
    std::vector< std::vector<BlockInteractions<ParticleGroupClass>>> externalInteractionsLeafLevelMpi;

    struct RemoteHandle{
        RemoteHandle() : ptrSymb(nullptr), ptrUp(nullptr), ptrDown(nullptr){
            memset(&handleSymb, 0, sizeof(handleSymb));
            memset(&handleUp, 0, sizeof(handleUp));
            memset(&handleDown, 0, sizeof(handleDown));
        }

        unsigned char * ptrSymb;
        starpu_data_handle_t handleSymb;
        unsigned char * ptrUp;
        starpu_data_handle_t handleUp;
        unsigned char * ptrDown;
        starpu_data_handle_t handleDown;
    };

    std::vector<std::vector<RemoteHandle>> remoteCellGroups;
    std::vector<RemoteHandle> remoteParticleGroupss;

    void buildRemoteInteractionsAndHandles(){
        cleanHandleMpi();

        // We need to have information about all other blocks
        std::unique_ptr<int[]> nbBlocksPerLevel(new int[tree->getHeight()]);
        nbBlocksPerLevel[0] = 0;
        for(int idxLevel = 1 ; idxLevel < tree->getHeight() ; ++idxLevel){
            nbBlocksPerLevel[idxLevel] = tree->getNbCellGroupAtLevel(idxLevel);
        }
        // Exchange the number of blocks per proc
        nbBlocksPerLevelAll.resize(tree->getHeight() * comm.processCount());
        FMpi::Assert(MPI_Allgather(nbBlocksPerLevel.get(), tree->getHeight(), MPI_INT,
                                   nbBlocksPerLevelAll.data(), tree->getHeight(), MPI_INT,
                                   comm.getComm()), __LINE__);
        // Compute the number of blocks before mine
        nbBlocksBeforeMinPerLevel.resize(tree->getHeight());
        for(int idxLevel = 1 ; idxLevel < tree->getHeight() ; ++idxLevel){
            nbBlocksBeforeMinPerLevel[idxLevel] = 0;
            for(int idxProc = 0 ; idxProc < comm.processId() ; ++idxProc){
                nbBlocksBeforeMinPerLevel[idxLevel] += nbBlocksPerLevelAll[idxProc*tree->getHeight() + idxLevel];
            }
        }
        // Prepare the block infos
        processesBlockInfos.resize(tree->getHeight());
        std::unique_ptr<int[]> recvBlocksCount(new int[comm.processCount()]);
        std::unique_ptr<int[]> recvBlockDispl(new int[comm.processCount()]);
        // Exchange the block info per level
        for(int idxLevel = 1 ; idxLevel < tree->getHeight() ; ++idxLevel){
            // Count the total number of blocks
            int nbBlocksInLevel = 0;
            recvBlockDispl[0] = 0;
            for(int idxProc = 0 ; idxProc < comm.processCount() ; ++idxProc){
                nbBlocksInLevel += nbBlocksPerLevelAll[idxProc*tree->getHeight() + idxLevel];
                // Count and displacement for the MPI all gatherv
                recvBlocksCount[idxProc] = nbBlocksPerLevelAll[idxProc*tree->getHeight() + idxLevel] * int(sizeof(BlockDescriptor));
                if(idxProc) recvBlockDispl[idxProc] = recvBlockDispl[idxProc-1] + recvBlocksCount[idxProc-1];
            }
            processesBlockInfos[idxLevel].resize(nbBlocksInLevel);
            // Fill my blocks
            std::vector<BlockDescriptor> myBlocksAtLevel;
            myBlocksAtLevel.resize(nbBlocksPerLevel[idxLevel]);
            FAssertLF(tree->getNbCellGroupAtLevel(idxLevel) == int(myBlocksAtLevel.size()));
            FAssertLF(nbBlocksPerLevel[idxLevel] == nbBlocksPerLevelAll[comm.processId()*tree->getHeight() + idxLevel]);

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                CellContainerClass*const currentCells = tree->getCellGroup(idxLevel, idxGroup);
                myBlocksAtLevel[idxGroup].firstIndex = currentCells->getStartingIndex();
                myBlocksAtLevel[idxGroup].lastIndex  = currentCells->getEndingIndex();
                myBlocksAtLevel[idxGroup].owner = comm.processId();
                myBlocksAtLevel[idxGroup].bufferSizeSymb = currentCells->getBufferSizeInByte();
                myBlocksAtLevel[idxGroup].bufferSizeUp   = currentCells->getMultipoleBufferSizeInByte();
                myBlocksAtLevel[idxGroup].bufferSizeDown = currentCells->getLocalBufferSizeInByte();

                if(idxLevel == tree->getHeight() - 1){
                    myBlocksAtLevel[idxGroup].leavesBufferSize = tree->getParticleGroup(idxGroup)->getBufferSizeInByte();
                }
                else{
                    myBlocksAtLevel[idxGroup].leavesBufferSize = 0;
                }
            }
            // Exchange with all other
            FMpi::Assert(MPI_Allgatherv(myBlocksAtLevel.data(), int(myBlocksAtLevel.size()*sizeof(BlockDescriptor)), MPI_BYTE,
                                        processesBlockInfos[idxLevel].data(), recvBlocksCount.get(), recvBlockDispl.get(), MPI_BYTE,
                                        comm.getComm()), __LINE__);
        }
        // Prepare remate ptr and handles
        remoteCellGroups.resize( tree->getHeight() );
        for(int idxLevel = 1 ; idxLevel < tree->getHeight() ; ++idxLevel){
            remoteCellGroups[idxLevel].resize( processesBlockInfos[idxLevel].size());
        }
        remoteParticleGroupss.resize(processesBlockInfos[tree->getHeight()-1].size());

        // From now we have the number of blocks for all process
        // we also have the size of the blocks therefor we can
        // create the handles we need
        // We will now detect the relation between our blocks and others
        // During the M2M (which is the same for the L2L)
        // During the M2L and during the P2P
        // I need to insert the task that read my data or that write the data I need.
        // M2L
        externalInteractionsAllLevelMpi.clear();
        externalInteractionsAllLevelMpi.resize(tree->getHeight());
        for(int idxLevel = tree->getHeight()-1 ; idxLevel >= 2 ; --idxLevel){
            // From this level there are no more blocks
            if(tree->getNbCellGroupAtLevel(idxLevel) == 0){
                // We stop here
                break;
            }
            // What are my morton interval at this level
            const MortonIndex myFirstIndex = tree->getCellGroup(idxLevel, 0)->getStartingIndex();
            const MortonIndex myLastIndex = tree->getCellGroup(idxLevel, tree->getNbCellGroupAtLevel(idxLevel)-1)->getEndingIndex();

            externalInteractionsAllLevelMpi[idxLevel].resize(tree->getNbCellGroupAtLevel(idxLevel));

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                CellContainerClass* currentCells = tree->getCellGroup(idxLevel, idxGroup);

                std::vector<BlockInteractions<CellContainerClass>>* externalInteractions = &externalInteractionsAllLevelMpi[idxLevel][idxGroup];

                #pragma omp task default(none) firstprivate(idxGroup, currentCells, idxLevel, externalInteractions)
                {
                    std::vector<OutOfBlockInteraction> outsideInteractions;
                    const MortonIndex blockStartIdx = currentCells->getStartingIndex();
                    const MortonIndex blockEndIdx   = currentCells->getEndingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                        if(currentCells->exists(mindex)){
                            const typename CellContainerClass::CompleteCellClass cell = currentCells->getCompleteCell(mindex);
                            FAssertLF(cell.getMortonIndex() == mindex);
                            MortonIndex interactionsIndexes[189];
                            int interactionsPosition[189];
                            const FTreeCoordinate coord(cell.getCoordinate());
                            int counter = coord.getInteractionNeighbors(idxLevel,interactionsIndexes,interactionsPosition);

                            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                                // This interactions need a block owned by someoneelse
                                if(interactionsIndexes[idxInter] < myFirstIndex || myLastIndex <= interactionsIndexes[idxInter]){
                                    OutOfBlockInteraction property;
                                    property.insideIndex = mindex;
                                    property.outIndex    = interactionsIndexes[idxInter];
                                    property.outPosition = interactionsPosition[idxInter];
                                    outsideInteractions.push_back(property);
                                }
                            }
                        }
                    }

                    // Manage outofblock interaction
                    FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                    int currentOutInteraction = 0;
                    for(int idxOtherGroup = 0 ; idxOtherGroup < int(processesBlockInfos[idxLevel].size())
                                                && currentOutInteraction < int(outsideInteractions.size()) ; ++idxOtherGroup){
                        // Skip my blocks
                        if(idxOtherGroup == nbBlocksBeforeMinPerLevel[idxLevel]){
                            idxOtherGroup += tree->getNbCellGroupAtLevel(idxLevel);
                            if(idxOtherGroup == int(processesBlockInfos[idxLevel].size())){
                                break;
                            }
                            FAssertLF(idxOtherGroup < int(processesBlockInfos[idxLevel].size()));
                        }

                        const MortonIndex blockStartIdx = processesBlockInfos[idxLevel][idxOtherGroup].firstIndex;
                        const MortonIndex blockEndIdx   = processesBlockInfos[idxLevel][idxOtherGroup].lastIndex;

                        while(currentOutInteraction < int(outsideInteractions.size()) && outsideInteractions[currentOutInteraction].outIndex < blockStartIdx){
                            currentOutInteraction += 1;
                        }

                        int lastOutInteraction = currentOutInteraction;
                        while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdx){
                            lastOutInteraction += 1;
                        }

                        // Create interactions
                        const int nbInteractionsBetweenBlocks = (lastOutInteraction-currentOutInteraction);
                        if(nbInteractionsBetweenBlocks){
                            if(remoteCellGroups[idxLevel][idxOtherGroup].ptrSymb == nullptr){
                                #pragma omp critical(CreateM2LRemotes)
                                {
                                    if(remoteCellGroups[idxLevel][idxOtherGroup].ptrSymb == nullptr){
                                        const int nbBytesInBlockSymb = processesBlockInfos[idxLevel][idxOtherGroup].bufferSizeSymb;
                                        unsigned char* memoryBlockSymb = (unsigned char*)FAlignedMemory::AllocateBytes<32>(nbBytesInBlockSymb);
                                        remoteCellGroups[idxLevel][idxOtherGroup].ptrSymb = memoryBlockSymb;
                                        starpu_variable_data_register(&remoteCellGroups[idxLevel][idxOtherGroup].handleSymb, 0,
                                                                      (uintptr_t)remoteCellGroups[idxLevel][idxOtherGroup].ptrSymb, nbBytesInBlockSymb);
                                        const int nbBytesInBlockUp = processesBlockInfos[idxLevel][idxOtherGroup].bufferSizeUp;
                                        unsigned char* memoryBlockUp = (unsigned char*)FAlignedMemory::AllocateBytes<32>(nbBytesInBlockUp);
                                        remoteCellGroups[idxLevel][idxOtherGroup].ptrUp = memoryBlockUp;
                                        starpu_variable_data_register(&remoteCellGroups[idxLevel][idxOtherGroup].handleUp, 0,
                                                                      (uintptr_t)remoteCellGroups[idxLevel][idxOtherGroup].ptrUp, nbBytesInBlockUp);
                                    }
                                }
                            }

                            externalInteractions->emplace_back();
                            BlockInteractions<CellContainerClass>* interactions = &externalInteractions->back();
                            //interactions->otherBlock = remoteCellGroups[idxLevel][idxOtherGroup].ptr;
                            interactions->otherBlockId = idxOtherGroup;
                            interactions->interactions.resize(nbInteractionsBetweenBlocks);
                            std::copy(outsideInteractions.begin() + currentOutInteraction,
                                      outsideInteractions.begin() + lastOutInteraction,
                                      interactions->interactions.begin());
                        }

                        currentOutInteraction = lastOutInteraction;
                    }
                }
            }
        }
        // P2P
        // We create one big vector per block
        {
            const MortonIndex myFirstIndex = tree->getParticleGroup(0)->getStartingIndex();
            const MortonIndex myLastIndex = tree->getParticleGroup(tree->getNbParticleGroup()-1)->getEndingIndex();

            externalInteractionsLeafLevelMpi.clear();
            externalInteractionsLeafLevelMpi.resize(tree->getNbParticleGroup());
            for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
                // Create the vector
                ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);

                std::vector<BlockInteractions<ParticleGroupClass>>* externalInteractions = &externalInteractionsLeafLevelMpi[idxGroup];

                #pragma omp task default(none) firstprivate(idxGroup, containers, externalInteractions)
                { // Can be a task(inout:iterCells)
                    std::vector<OutOfBlockInteraction> outsideInteractions;
                    const MortonIndex blockStartIdx = containers->getStartingIndex();
                    const MortonIndex blockEndIdx   = containers->getEndingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                        // ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(mindex);
                        if(containers->exists(mindex)){
                            MortonIndex interactionsIndexes[26];
                            int interactionsPosition[26];
                            FTreeCoordinate coord(mindex, tree->getHeight()-1);
                            int counter = coord.getNeighborsIndexes(tree->getHeight(),interactionsIndexes,interactionsPosition);

                            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                                if(interactionsIndexes[idxInter] < myFirstIndex ||
                                        myLastIndex <= interactionsIndexes[idxInter]){
                                    OutOfBlockInteraction property;
                                    property.insideIndex = mindex;
                                    property.outIndex    = interactionsIndexes[idxInter];
                                    property.outPosition = interactionsPosition[idxInter];
                                    outsideInteractions.push_back(property);
                                }
                            }
                        }
                    }

                    // Sort to match external order
                    FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                    int currentOutInteraction = 0;
                    for(int idxOtherGroup = 0 ; idxOtherGroup < int(processesBlockInfos[tree->getHeight()-1].size())
                                                && currentOutInteraction < int(outsideInteractions.size()) ; ++idxOtherGroup){
                        // Skip my blocks
                        if(idxOtherGroup == nbBlocksBeforeMinPerLevel[tree->getHeight()-1]){
                            idxOtherGroup += tree->getNbCellGroupAtLevel(tree->getHeight()-1);
                            if(idxOtherGroup == int(processesBlockInfos[tree->getHeight()-1].size())){
                                break;
                            }
                            FAssertLF(idxOtherGroup < int(processesBlockInfos[tree->getHeight()-1].size()));
                        }

                        const MortonIndex blockStartIdx = processesBlockInfos[tree->getHeight()-1][idxOtherGroup].firstIndex;
                        const MortonIndex blockEndIdx   = processesBlockInfos[tree->getHeight()-1][idxOtherGroup].lastIndex;

                        while(currentOutInteraction < int(outsideInteractions.size()) && outsideInteractions[currentOutInteraction].outIndex < blockStartIdx){
                            currentOutInteraction += 1;
                        }

                        int lastOutInteraction = currentOutInteraction;
                        while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdx){
                            lastOutInteraction += 1;
                        }

                        const int nbInteractionsBetweenBlocks = (lastOutInteraction-currentOutInteraction);
                        if(nbInteractionsBetweenBlocks){
                            if(remoteParticleGroupss[idxOtherGroup].ptrSymb == nullptr){
                                #pragma omp critical(CreateM2LRemotes)
                                {
                                    if(remoteParticleGroupss[idxOtherGroup].ptrSymb == nullptr){
                                        const int nbBytesInBlock = processesBlockInfos[tree->getHeight()-1][idxOtherGroup].leavesBufferSize;
                                        unsigned char* memoryBlock = (unsigned char*)FAlignedMemory::AllocateBytes<32>(nbBytesInBlock);
                                        remoteParticleGroupss[idxOtherGroup].ptrSymb = memoryBlock;
                                        starpu_variable_data_register(&remoteParticleGroupss[idxOtherGroup].handleSymb, 0,
                                                                      (uintptr_t)remoteParticleGroupss[idxOtherGroup].ptrSymb, nbBytesInBlock);
                                    }
                                }
                            }

                            externalInteractions->emplace_back();
                            BlockInteractions<ParticleGroupClass>* interactions = &externalInteractions->back();
                            //interactions->otherBlock = remoteParticleGroupss[idxOtherGroup].ptr;
                            interactions->otherBlockId = idxOtherGroup;
                            interactions->interactions.resize(nbInteractionsBetweenBlocks);
                            std::copy(outsideInteractions.begin() + currentOutInteraction,
                                      outsideInteractions.begin() + lastOutInteraction,
                                      interactions->interactions.begin());
                        }

                        currentOutInteraction = lastOutInteraction;
                    }
                }
            }
        }
    }

    struct MpiDependency{
        int src;
        int dest;
        int level;
        int globalBlockId;
    };

    std::vector<MpiDependency> toSend;

    void postRecvAllocatedBlocks(){
        std::vector<MpiDependency> toRecv;
        FAssertLF(tree->getHeight() == int(remoteCellGroups.size()));
        for(int idxLevel = 0 ; idxLevel < tree->getHeight() ; ++idxLevel){
            for(int idxHandle = 0 ; idxHandle < int(remoteCellGroups[idxLevel].size()) ; ++idxHandle){
                if(remoteCellGroups[idxLevel][idxHandle].ptrSymb){
                    FAssertLF(remoteCellGroups[idxLevel][idxHandle].ptrUp);
                    FLOG(FLog::Controller << "[SMpi] " << idxLevel << " Post a recv during M2L for Idx " << processesBlockInfos[idxLevel][idxHandle].firstIndex <<
                         " and dest is " << processesBlockInfos[idxLevel][idxHandle].owner << " tag " << getTag(idxLevel,processesBlockInfos[idxLevel][idxHandle].firstIndex, 0) << "\n");
                    FLOG(FLog::Controller << "[SMpi] " << idxLevel << " Post a recv during M2L for Idx " << processesBlockInfos[idxLevel][idxHandle].firstIndex <<
                         " and dest is " << processesBlockInfos[idxLevel][idxHandle].owner << " tag " << getTag(idxLevel,processesBlockInfos[idxLevel][idxHandle].firstIndex, 1) << "\n");

                    starpu_mpi_irecv_detached( remoteCellGroups[idxLevel][idxHandle].handleSymb,
                                                processesBlockInfos[idxLevel][idxHandle].owner,
                                                getTag(idxLevel,processesBlockInfos[idxLevel][idxHandle].firstIndex, 0),
                                                comm.getComm(), 0, 0 );
                    starpu_mpi_irecv_detached( remoteCellGroups[idxLevel][idxHandle].handleUp,
                                                processesBlockInfos[idxLevel][idxHandle].owner,
                                                getTag(idxLevel,processesBlockInfos[idxLevel][idxHandle].firstIndex, 1),
                                                comm.getComm(), 0, 0 );

                    toRecv.push_back({processesBlockInfos[idxLevel][idxHandle].owner,
                                        comm.processId(), idxLevel, idxHandle});
                }
            }
        }
        {
            for(int idxHandle = 0 ; idxHandle < int(remoteParticleGroupss.size()) ; ++idxHandle){
                if(remoteParticleGroupss[idxHandle].ptrSymb){
                    FLOG(FLog::Controller << "[SMpi] Post a recv during P2P for Idx " << processesBlockInfos[tree->getHeight()-1][idxHandle].firstIndex <<
                         " and dest is " << processesBlockInfos[tree->getHeight()-1][idxHandle].owner << " tag " << getTag(tree->getHeight(),processesBlockInfos[tree->getHeight()-1][idxHandle].firstIndex, 0) << "\n");

                    starpu_mpi_irecv_detached( remoteParticleGroupss[idxHandle].handleSymb,
                                                processesBlockInfos[tree->getHeight()-1][idxHandle].owner,
                                                getTag(tree->getHeight(),processesBlockInfos[tree->getHeight()-1][idxHandle].firstIndex, 0),
                                                comm.getComm(), 0, 0 );

                    toRecv.push_back({processesBlockInfos[tree->getHeight()-1][idxHandle].owner,
                                        comm.processId(), tree->getHeight(), idxHandle});
                }
            }
        }

        FQuickSort<MpiDependency, int>::QsSequential(toRecv.data(),int(toRecv.size()),[](const MpiDependency& d1, const MpiDependency& d2){
            return d1.src <= d2.src;
        });

        std::unique_ptr<int[]> nbBlocksToRecvFromEach(new int[comm.processCount()]);
        memset(nbBlocksToRecvFromEach.get(), 0, sizeof(int)*comm.processCount());
        for(int idxDep = 0 ; idxDep < int(toRecv.size()) ; ++idxDep){
            nbBlocksToRecvFromEach[toRecv[idxDep].src] += 1;
        }

        FAssertLF(nbBlocksToRecvFromEach[comm.processId()] == 0);
        int offset = 0;

        for(int idxProc = 0 ; idxProc < comm.processCount() ; ++idxProc){
            if(idxProc == comm.processId()){
                // How much to send to each
                std::unique_ptr<int[]> nbBlocksToSendToEach(new int[comm.processCount()]);
                FMpi::Assert(MPI_Gather(&nbBlocksToRecvFromEach[idxProc], 1,
                                 MPI_INT, nbBlocksToSendToEach.get(), 1,
                                 MPI_INT, idxProc, comm.getComm() ), __LINE__);

                std::unique_ptr<int[]> displs(new int[comm.processCount()]);
                displs[0] = 0;
                for(int idxProc = 1 ; idxProc < comm.processCount() ; ++idxProc){
                    displs[idxProc] = displs[idxProc-1] + nbBlocksToSendToEach[idxProc-1];
                }
                toSend.resize(displs[comm.processCount()-1] + nbBlocksToSendToEach[comm.processCount()-1]);

                // We work in bytes
                for(int idxProc = 0 ; idxProc < comm.processCount() ; ++idxProc){
                    nbBlocksToSendToEach[idxProc] *= sizeof(MpiDependency);
                    displs[idxProc] *= sizeof(MpiDependency);
                }

                FMpi::Assert(MPI_Gatherv( nullptr, 0, MPI_BYTE,
                                 toSend.data(),
                                 nbBlocksToSendToEach.get(), displs.get(),
                                 MPI_BYTE, idxProc, comm.getComm()), __LINE__);
            }
            else{
                FMpi::Assert(MPI_Gather(&nbBlocksToRecvFromEach[idxProc], 1,
                                 MPI_INT, 0, 0, MPI_INT, idxProc, comm.getComm() ), __LINE__);
                FMpi::Assert(MPI_Gatherv(
                                 &toRecv[offset], int(nbBlocksToRecvFromEach[idxProc]*sizeof(MpiDependency)), MPI_BYTE,
                                 0, 0, 0, MPI_BYTE, idxProc, comm.getComm() ), __LINE__);

                offset += nbBlocksToRecvFromEach[idxProc];
            }
        }
    }

    void insertParticlesSend(){
        for(int idxSd = 0 ; idxSd < int(toSend.size()) ; ++idxSd){
            const MpiDependency sd = toSend[idxSd];
            if(sd.level == tree->getHeight()){
                const int localId = sd.globalBlockId - nbBlocksBeforeMinPerLevel[tree->getHeight()-1];
                FAssertLF(sd.src == comm.processId());
                FAssertLF(0 <= localId);
                FAssertLF(localId < tree->getNbParticleGroup());

                FLOG(FLog::Controller << "[SMpi] Post a send during P2P for Idx " << tree->getParticleGroup(localId)->getStartingIndex() <<
                     " and dest is " << sd.dest << " tag " << getTag(tree->getHeight(),tree->getParticleGroup(localId)->getStartingIndex(), 0) <<  "\n");

                starpu_mpi_isend_detached( particleHandles[localId].symb, sd.dest,
                        getTag(tree->getHeight(),tree->getParticleGroup(localId)->getStartingIndex(), 0),
                        comm.getComm(), 0/*callback*/, 0/*arg*/ );
            }
        }
    }

    void insertCellsSend(){
        for(int idxSd = 0 ; idxSd < int(toSend.size()) ; ++idxSd){
            const MpiDependency sd = toSend[idxSd];
            if(sd.level != tree->getHeight()){
                const int localId = sd.globalBlockId - nbBlocksBeforeMinPerLevel[sd.level];
                FAssertLF(sd.src == comm.processId());
                FAssertLF(0 <= localId);
                FAssertLF(localId < tree->getNbCellGroupAtLevel(sd.level));

                FLOG(FLog::Controller << "[SMpi] " << sd.level << " Post a send during M2L for Idx " << tree->getCellGroup(sd.level, localId)->getStartingIndex() <<
                     " and dest is " << sd.dest << " tag " << getTag(sd.level,tree->getCellGroup(sd.level, localId)->getStartingIndex(), 0) << "\n");
                FLOG(FLog::Controller << "[SMpi] " << sd.level << " Post a send during M2L for Idx " << tree->getCellGroup(sd.level, localId)->getStartingIndex() <<
                     " and dest is " << sd.dest << " tag " << getTag(sd.level,tree->getCellGroup(sd.level, localId)->getStartingIndex(), 1) << "\n");

                starpu_mpi_isend_detached( cellHandles[sd.level][localId].symb, sd.dest,
                        getTag(sd.level,tree->getCellGroup(sd.level, localId)->getStartingIndex(), 0),
                        comm.getComm(), 0/*callback*/, 0/*arg*/ );
                starpu_mpi_isend_detached( cellHandles[sd.level][localId].up, sd.dest,
                        getTag(sd.level,tree->getCellGroup(sd.level, localId)->getStartingIndex(), 1),
                        comm.getComm(), 0/*callback*/, 0/*arg*/ );
            }
        }
    }

    void cleanHandleMpi(){
        for(int idxLevel = 0 ; idxLevel < int(remoteCellGroups.size()) ; ++idxLevel){
            for(int idxHandle = 0 ; idxHandle < int(remoteCellGroups[idxLevel].size()) ; ++idxHandle){
                if(remoteCellGroups[idxLevel][idxHandle].ptrSymb){
                    starpu_data_unregister(remoteCellGroups[idxLevel][idxHandle].handleSymb);
                    starpu_data_unregister(remoteCellGroups[idxLevel][idxHandle].handleUp);
                    FAlignedMemory::DeallocBytes(remoteCellGroups[idxLevel][idxHandle].ptrSymb);
                    FAlignedMemory::DeallocBytes(remoteCellGroups[idxLevel][idxHandle].ptrUp);

                    if(remoteCellGroups[idxLevel][idxHandle].ptrDown){
                        starpu_data_unregister(remoteCellGroups[idxLevel][idxHandle].handleDown);
                        FAlignedMemory::DeallocBytes(remoteCellGroups[idxLevel][idxHandle].ptrDown);
                    }
                }
            }
            remoteCellGroups[idxLevel].clear();
        }
        {
            for(int idxHandle = 0 ; idxHandle < int(remoteParticleGroupss.size()) ; ++idxHandle){
                if(remoteParticleGroupss[idxHandle].ptrSymb){
                    starpu_data_unregister(remoteParticleGroupss[idxHandle].handleSymb);
                    FAlignedMemory::DeallocBytes(remoteParticleGroupss[idxHandle].ptrSymb);
                }
            }
            remoteParticleGroupss.clear();
        }
    }

    ////////////////////////////////////////////////////////////////////////////

    /** Reset the handles array and create new ones to define
     * in a starpu way each block of data
     */
    void buildHandles(){
        cleanHandle();

        for(int idxLevel = 2 ; idxLevel < tree->getHeight() ; ++idxLevel){
            cellHandles[idxLevel].resize(tree->getNbCellGroupAtLevel(idxLevel));

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                const CellContainerClass* currentCells = tree->getCellGroup(idxLevel, idxGroup);
                starpu_variable_data_register(&cellHandles[idxLevel][idxGroup].symb, 0,
                                              (uintptr_t)currentCells->getRawBuffer(), currentCells->getBufferSizeInByte());
                starpu_variable_data_register(&cellHandles[idxLevel][idxGroup].up, 0,
                                              (uintptr_t)currentCells->getRawMultipoleBuffer(), currentCells->getMultipoleBufferSizeInByte());
                starpu_variable_data_register(&cellHandles[idxLevel][idxGroup].down, 0,
                                              (uintptr_t)currentCells->getRawLocalBuffer(), currentCells->getLocalBufferSizeInByte());
            }
        }
        {
            particleHandles.resize(tree->getNbParticleGroup());
            for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
                ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);
                starpu_variable_data_register(&particleHandles[idxGroup].symb, 0,
                                              (uintptr_t)containers->getRawBuffer(), containers->getBufferSizeInByte());
                starpu_variable_data_register(&particleHandles[idxGroup].down, 0,
                                              (uintptr_t)containers->getRawAttributesBuffer(), containers->getAttributesBufferSizeInByte());
            }
        }
    }

    /**
     * This function is creating the interactions vector between blocks.
     * It fills externalInteractionsAllLevel and externalInteractionsLeafLevel.
     * Warning, the omp task for now are using the class attributes!
     *
     */
    void buildExternalInteractionVecs(){
        FLOG( FTic timer; FTic leafTimer; FTic cellTimer; );
        // Reset interactions
        externalInteractionsAllLevel.clear();
        externalInteractionsLeafLevel.clear();
        // One per level + leaf level
        externalInteractionsAllLevel.resize(tree->getHeight());

        // First leaf level
        {
            // We create one big vector per block
            externalInteractionsLeafLevel.resize(tree->getNbParticleGroup());

            for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
                // Create the vector
                ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);

                std::vector<BlockInteractions<ParticleGroupClass>>* externalInteractions = &externalInteractionsLeafLevel[idxGroup];

                #pragma omp task default(none) firstprivate(idxGroup, containers, externalInteractions)
                { // Can be a task(inout:iterCells)
                    std::vector<OutOfBlockInteraction> outsideInteractions;
                    const MortonIndex blockStartIdx = containers->getStartingIndex();
                    const MortonIndex blockEndIdx   = containers->getEndingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                        if(containers->exists(mindex)){
                            MortonIndex interactionsIndexes[26];
                            int interactionsPosition[26];
                            FTreeCoordinate coord(mindex, tree->getHeight()-1);
                            int counter = coord.getNeighborsIndexes(tree->getHeight(),interactionsIndexes,interactionsPosition);

                            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                    // Inside block interaction, do nothing
                                }
                                else if(interactionsIndexes[idxInter] < mindex){
                                    OutOfBlockInteraction property;
                                    property.insideIndex = mindex;
                                    property.outIndex    = interactionsIndexes[idxInter];
                                    property.outPosition = interactionsPosition[idxInter];
                                    outsideInteractions.push_back(property);
                                }
                            }
                        }
                    }

                    // Sort to match external order
                    FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                    int currentOutInteraction = 0;
                    for(int idxLeftGroup = 0 ; idxLeftGroup < idxGroup && currentOutInteraction < int(outsideInteractions.size()) ; ++idxLeftGroup){
                        ParticleGroupClass* leftContainers = tree->getParticleGroup(idxLeftGroup);
                        const MortonIndex blockStartIdx    = leftContainers->getStartingIndex();
                        const MortonIndex blockEndIdx      = leftContainers->getEndingIndex();

                        while(currentOutInteraction < int(outsideInteractions.size()) && outsideInteractions[currentOutInteraction].outIndex < blockStartIdx){
                            currentOutInteraction += 1;
                        }

                        int lastOutInteraction = currentOutInteraction;
                        while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdx){
                            lastOutInteraction += 1;
                        }

                        const int nbInteractionsBetweenBlocks = (lastOutInteraction-currentOutInteraction);
                        if(nbInteractionsBetweenBlocks){
                            externalInteractions->emplace_back();
                            BlockInteractions<ParticleGroupClass>* interactions = &externalInteractions->back();
                            interactions->otherBlock = leftContainers;
                            interactions->otherBlockId = idxLeftGroup;
                            interactions->interactions.resize(nbInteractionsBetweenBlocks);
                            std::copy(outsideInteractions.begin() + currentOutInteraction,
                                      outsideInteractions.begin() + lastOutInteraction,
                                      interactions->interactions.begin());
                        }

                        currentOutInteraction = lastOutInteraction;
                    }
                }
            }
        }
        FLOG( leafTimer.tac(); );
        FLOG( cellTimer.tic(); );
        {
            for(int idxLevel = tree->getHeight()-1 ; idxLevel >= 2 ; --idxLevel){
                externalInteractionsAllLevel[idxLevel].resize(tree->getNbCellGroupAtLevel(idxLevel));

                for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                    CellContainerClass* currentCells = tree->getCellGroup(idxLevel, idxGroup);

                    std::vector<BlockInteractions<CellContainerClass>>* externalInteractions = &externalInteractionsAllLevel[idxLevel][idxGroup];

                    #pragma omp task default(none) firstprivate(idxGroup, currentCells, idxLevel, externalInteractions)
                    {
                        std::vector<OutOfBlockInteraction> outsideInteractions;
                        const MortonIndex blockStartIdx = currentCells->getStartingIndex();
                        const MortonIndex blockEndIdx   = currentCells->getEndingIndex();

                        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                            if(currentCells->exists(mindex)){
                                typename CellContainerClass::CompleteCellClass cell = currentCells->getCompleteCell(mindex);
                                FAssertLF(cell.getMortonIndex() == mindex);
                                MortonIndex interactionsIndexes[189];
                                int interactionsPosition[189];
                                const FTreeCoordinate coord(cell.getCoordinate());
                                int counter = coord.getInteractionNeighbors(idxLevel,interactionsIndexes,interactionsPosition);

                                for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                                    if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                        // Nothing to do
                                    }
                                    else if(interactionsIndexes[idxInter] < mindex){
                                        OutOfBlockInteraction property;
                                        property.insideIndex = mindex;
                                        property.outIndex    = interactionsIndexes[idxInter];
                                        property.outPosition = interactionsPosition[idxInter];
                                        outsideInteractions.push_back(property);
                                    }
                                }
                            }
                        }

                        // Manage outofblock interaction
                        FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                        int currentOutInteraction = 0;
                        for(int idxLeftGroup = 0 ; idxLeftGroup < idxGroup && currentOutInteraction < int(outsideInteractions.size()) ; ++idxLeftGroup){
                            CellContainerClass* leftCells   = tree->getCellGroup(idxLevel, idxLeftGroup);
                            const MortonIndex blockStartIdx = leftCells->getStartingIndex();
                            const MortonIndex blockEndIdx   = leftCells->getEndingIndex();

                            while(currentOutInteraction < int(outsideInteractions.size()) && outsideInteractions[currentOutInteraction].outIndex < blockStartIdx){
                                currentOutInteraction += 1;
                            }

                            int lastOutInteraction = currentOutInteraction;
                            while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdx){
                                lastOutInteraction += 1;
                            }

                            // Create interactions
                            const int nbInteractionsBetweenBlocks = (lastOutInteraction-currentOutInteraction);
                            if(nbInteractionsBetweenBlocks){
                                externalInteractions->emplace_back();
                                BlockInteractions<CellContainerClass>* interactions = &externalInteractions->back();
                                interactions->otherBlock = leftCells;
                                interactions->otherBlockId = idxLeftGroup;
                                interactions->interactions.resize(nbInteractionsBetweenBlocks);
                                std::copy(outsideInteractions.begin() + currentOutInteraction,
                                          outsideInteractions.begin() + lastOutInteraction,
                                          interactions->interactions.begin());
                            }

                            currentOutInteraction = lastOutInteraction;
                        }
                    }
                }
            }
        }
        FLOG( cellTimer.tac(); );

        #pragma omp taskwait

        FLOG( FLog::Controller << "\t\t Prepare in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t Prepare at leaf level in   " << leafTimer.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t Prepare at other levels in " << cellTimer.elapsed() << "s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Bottom Pass
    /////////////////////////////////////////////////////////////////////////////////////

    void bottomPass(){
        FLOG( FTic timer; );

        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            starpu_insert_task(&p2m_cl,
                    STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                    STARPU_R, cellHandles[tree->getHeight()-1][idxGroup].symb,
                    STARPU_RW, cellHandles[tree->getHeight()-1][idxGroup].up,
                    STARPU_R, particleHandles[idxGroup].symb,
                    0);
        }

        FLOG( FLog::Controller << "\t\t bottomPass in " << timer.tacAndElapsed() << "s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Upward Pass
    /////////////////////////////////////////////////////////////////////////////////////

    void upwardPass(){
        FLOG( FTic timer; );
        for(int idxLevel = tree->getHeight()-2 ; idxLevel >= 2 ; --idxLevel){
            int idxSubGroup = 0;

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel)
                    && idxSubGroup < tree->getNbCellGroupAtLevel(idxLevel+1) ; ++idxGroup){
                CellContainerClass*const currentCells = tree->getCellGroup(idxLevel, idxGroup);

                struct starpu_task* const task = starpu_task_create();
                task->dyn_handles = (starpu_data_handle_t*)malloc(sizeof(starpu_data_handle_t)*20);
                task->dyn_handles[0] = cellHandles[idxLevel][idxGroup].symb;
                task->dyn_handles[1] = cellHandles[idxLevel][idxGroup].up;

                // Skip current group if needed
                if( tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++idxSubGroup;
                    FAssertLF( idxSubGroup != tree->getNbCellGroupAtLevel(idxLevel+1) );
                    FAssertLF( (tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }

                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                task->dyn_handles[(nbSubCellGroups*2) + 2] = cellHandles[idxLevel+1][idxSubGroup].symb;
                task->dyn_handles[(nbSubCellGroups*2) + 3] = cellHandles[idxLevel+1][idxSubGroup].up;
                nbSubCellGroups += 1;

                while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (((currentCells->getEndingIndex()-1)<<3)+7)
                      && (idxSubGroup+1) != tree->getNbCellGroupAtLevel(idxLevel+1)
                      && tree->getCellGroup(idxLevel+1, idxSubGroup+1)->getStartingIndex() <= ((currentCells->getEndingIndex()-1)<<3)+7 ){
                    idxSubGroup += 1;
                    task->dyn_handles[(nbSubCellGroups*2) + 2] = cellHandles[idxLevel+1][idxSubGroup].symb;
                    task->dyn_handles[(nbSubCellGroups*2) + 3] = cellHandles[idxLevel+1][idxSubGroup].up;
                    nbSubCellGroups += 1;
                    FAssertLF( nbSubCellGroups <= 9 );
                }

                // put the right codelet
                task->cl = &m2m_cl[nbSubCellGroups-1];
                // put args values
                char *arg_buffer;
                size_t arg_buffer_size;
                starpu_codelet_pack_args((void**)&arg_buffer, &arg_buffer_size,
                                         STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                                         STARPU_VALUE, &nbSubCellGroups, sizeof(nbSubCellGroups),
                                         STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                                         0);
                task->cl_arg = arg_buffer;
                task->cl_arg_size = arg_buffer_size;
                FAssertLF(starpu_task_submit(task) == 0);
            }

            /////////////////////////////////////////////////////////////
            // Exchange for mpi
            /////////////////////////////////////////////////////////////
            // Manage the external operations
            // Find what to recv
            if(tree->getNbCellGroupAtLevel(idxLevel)){
                // Take last block at this level
                const CellContainerClass* currentCells = tree->getCellGroup(idxLevel, tree->getNbCellGroupAtLevel(idxLevel)-1);
                // Take the last cell index of the last block
                const MortonIndex myLastIdx = currentCells->getEndingIndex()-1;
                // Find the descriptor of the first block that belong to someone else at lower level
                const int firstOtherBlock = nbBlocksBeforeMinPerLevel[idxLevel+1] + tree->getNbCellGroupAtLevel(idxLevel+1);
                FAssertLF(processesBlockInfos[idxLevel+1][firstOtherBlock-1].owner == comm.processId());
                // Iterate while the block has our cell has parent
                int idxBlockToRecv = 0;
                while(firstOtherBlock + idxBlockToRecv < int(processesBlockInfos[idxLevel+1].size()) &&
                      myLastIdx == (processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].firstIndex >> 3)){

                    if(remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].ptrSymb == nullptr){
                        const int nbBytesInBlockSymb = processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].bufferSizeSymb;
                        unsigned char* memoryBlockSymb = (unsigned char*)FAlignedMemory::AllocateBytes<32>(nbBytesInBlockSymb);
                        remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].ptrSymb = memoryBlockSymb;
                        starpu_variable_data_register(&remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].handleSymb, 0,
                                                      (uintptr_t)remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].ptrSymb, nbBytesInBlockSymb);

                        const int nbBytesInBlockUp = processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].bufferSizeUp;
                        unsigned char* memoryBlockUp = (unsigned char*)FAlignedMemory::AllocateBytes<32>(nbBytesInBlockUp);
                        remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].ptrUp = memoryBlockUp;
                        starpu_variable_data_register(&remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].handleUp, 0,
                                                      (uintptr_t)remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].ptrUp, nbBytesInBlockUp);
                    }

                    FLOG(FLog::Controller << "[SMpi] Post a recv during M2M for Idx " << processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].firstIndex <<
                         " and owner is " << processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].owner << " tag " << getTag(idxLevel,processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].firstIndex, 0) << "\n");
                    FLOG(FLog::Controller << "[SMpi] Post a recv during M2M for Idx " << processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].firstIndex <<
                         " and owner is " << processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].owner << " tag " << getTag(idxLevel,processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].firstIndex, 1) << "\n");

                    starpu_mpi_irecv_detached ( remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].handleSymb,
                                                processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].owner,
                                                getTag(idxLevel,processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].firstIndex, 0),
                                                comm.getComm(), 0/*callback*/, 0/*arg*/ );
                    starpu_mpi_irecv_detached ( remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].handleUp,
                                                processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].owner,
                                                getTag(idxLevel,processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].firstIndex, 1),
                                                comm.getComm(), 0/*callback*/, 0/*arg*/ );


                    idxBlockToRecv += 1;
                }
                FAssertLF(idxBlockToRecv < 8);
                if(idxBlockToRecv){// Perform the work
                    struct starpu_task* const task = starpu_task_create();
                    task->dyn_handles = (starpu_data_handle_t*)malloc(sizeof(starpu_data_handle_t)*20);
                    task->dyn_handles[0] = cellHandles[idxLevel][tree->getNbCellGroupAtLevel(idxLevel)-1].symb;
                    task->dyn_handles[1] = cellHandles[idxLevel][tree->getNbCellGroupAtLevel(idxLevel)-1].up;

                    // Copy at max 8 groups
                    int nbSubCellGroups = 0;
                    while(nbSubCellGroups < idxBlockToRecv){
                        task->dyn_handles[(nbSubCellGroups*2) + 2] = remoteCellGroups[idxLevel+1][firstOtherBlock + nbSubCellGroups].handleSymb;
                        task->dyn_handles[(nbSubCellGroups*2) + 3] = remoteCellGroups[idxLevel+1][firstOtherBlock + nbSubCellGroups].handleUp;
                        nbSubCellGroups += 1;
                    }

                    // put the right codelet
                    task->cl = &m2m_cl[nbSubCellGroups-1];
                    // put args values
                    char *arg_buffer;
                    size_t arg_buffer_size;
                    starpu_codelet_pack_args((void**)&arg_buffer, &arg_buffer_size,
                                             STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                                             STARPU_VALUE, &nbSubCellGroups, sizeof(nbSubCellGroups),
                                             STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                                             0);
                    task->cl_arg = arg_buffer;
                    task->cl_arg_size = arg_buffer_size;
                    FAssertLF(starpu_task_submit(task) == 0);
                }
            }
            // Find what to send
            if(tree->getNbCellGroupAtLevel(idxLevel+1)
                    && nbBlocksBeforeMinPerLevel[idxLevel] != 0){
                // Take the first lower block
                const CellContainerClass* currentCells = tree->getCellGroup(idxLevel+1, 0);
                // Take its first index
                const MortonIndex myFirstChildIdx = currentCells->getStartingIndex();
                const MortonIndex missingParentIdx = (myFirstChildIdx>>3);
                // If no parent or the first parent is not the good one
                if(tree->getNbCellGroupAtLevel(idxLevel) == 0
                        || tree->getCellGroup(idxLevel, 0)->getStartingIndex() != missingParentIdx){
                    // Look if the parent is owned by another block
                    const int firstOtherBlock = nbBlocksBeforeMinPerLevel[idxLevel]-1;
                    FAssertLF(processesBlockInfos[idxLevel][firstOtherBlock].lastIndex-1 == missingParentIdx);
                    const int dest = processesBlockInfos[idxLevel][firstOtherBlock].owner;
                    int lowerIdxToSend = 0;
                    while(lowerIdxToSend != tree->getNbCellGroupAtLevel(idxLevel+1)
                          && missingParentIdx == (tree->getCellGroup(idxLevel+1, lowerIdxToSend)->getStartingIndex()>>3)){

                        FLOG(FLog::Controller << "[SMpi] Post a send during M2M for Idx " << tree->getCellGroup(idxLevel+1, lowerIdxToSend)->getStartingIndex() <<
                             " and dest is " << dest << " tag " << getTag(idxLevel,tree->getCellGroup(idxLevel+1, lowerIdxToSend)->getStartingIndex(), 0) << "\n");
                        FLOG(FLog::Controller << "[SMpi] Post a send during M2M for Idx " << tree->getCellGroup(idxLevel+1, lowerIdxToSend)->getStartingIndex() <<
                             " and dest is " << dest << " tag " << getTag(idxLevel,tree->getCellGroup(idxLevel+1, lowerIdxToSend)->getStartingIndex(), 1) << "\n");

                        starpu_mpi_isend_detached( cellHandles[idxLevel+1][lowerIdxToSend].symb, dest,
                                getTag(idxLevel,tree->getCellGroup(idxLevel+1, lowerIdxToSend)->getStartingIndex(), 0),
                                comm.getComm(), 0/*callback*/, 0/*arg*/ );
                        starpu_mpi_isend_detached( cellHandles[idxLevel+1][lowerIdxToSend].up, dest,
                                getTag(idxLevel,tree->getCellGroup(idxLevel+1, lowerIdxToSend)->getStartingIndex(), 1),
                                comm.getComm(), 0/*callback*/, 0/*arg*/ );
                        lowerIdxToSend += 1;
                    }

                }
            }
            /////////////////////////////////////////////////////////////
        }
        FLOG( FLog::Controller << "\t\t upwardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Transfer Pass Mpi
    /////////////////////////////////////////////////////////////////////////////////////

    void transferPassMpi(){
        FLOG( FTic timer; );
        for(int idxLevel = tree->getHeight()-1 ; idxLevel >= 2 ; --idxLevel){
            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                for(int idxInteraction = 0; idxInteraction < int(externalInteractionsAllLevelMpi[idxLevel][idxGroup].size()) ; ++idxInteraction){
                    const int interactionid = externalInteractionsAllLevelMpi[idxLevel][idxGroup][idxInteraction].otherBlockId;
                    const std::vector<OutOfBlockInteraction>* outsideInteractions = &externalInteractionsAllLevelMpi[idxLevel][idxGroup][idxInteraction].interactions;

                    starpu_insert_task(&m2l_cl_inout_mpi,
                            STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                            STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                            STARPU_VALUE, &outsideInteractions, sizeof(outsideInteractions),
                            STARPU_R, cellHandles[idxLevel][idxGroup].symb,
                            (STARPU_RW|STARPU_COMMUTE), cellHandles[idxLevel][idxGroup].down,
                            STARPU_R, remoteCellGroups[idxLevel][interactionid].handleSymb,
                            STARPU_R, remoteCellGroups[idxLevel][interactionid].handleUp,
                            0);
                }
            }
        }
        FLOG( FLog::Controller << "\t\t transferPassMpi in " << timer.tacAndElapsed() << "s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Transfer Pass
    /////////////////////////////////////////////////////////////////////////////////////

    void transferPass(){
        FLOG( FTic timer; );
        FLOG( FTic timerInBlock; FTic timerOutBlock; );
        for(int idxLevel = tree->getHeight()-1 ; idxLevel >= 2 ; --idxLevel){
            FLOG( timerInBlock.tic() );
            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                starpu_insert_task(&m2l_cl_in,
                        STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                        STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                                   STARPU_R, cellHandles[idxLevel][idxGroup].symb,
                                   STARPU_R, cellHandles[idxLevel][idxGroup].up,
                                   (STARPU_RW|STARPU_COMMUTE), cellHandles[idxLevel][idxGroup].down,
                        0);
            }
            FLOG( timerInBlock.tac() );

            FLOG( timerOutBlock.tic() );
            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                for(int idxInteraction = 0; idxInteraction < int(externalInteractionsAllLevel[idxLevel][idxGroup].size()) ; ++idxInteraction){
                    const int interactionid = externalInteractionsAllLevel[idxLevel][idxGroup][idxInteraction].otherBlockId;
                    const std::vector<OutOfBlockInteraction>* outsideInteractions = &externalInteractionsAllLevel[idxLevel][idxGroup][idxInteraction].interactions;

                    starpu_insert_task(&m2l_cl_inout,
                            STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                            STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                            STARPU_VALUE, &outsideInteractions, sizeof(outsideInteractions),
                               STARPU_R, cellHandles[idxLevel][idxGroup].symb,
                               STARPU_R, cellHandles[idxLevel][idxGroup].up,
                               (STARPU_RW|STARPU_COMMUTE), cellHandles[idxLevel][idxGroup].down,
                               STARPU_R, cellHandles[idxLevel][interactionid].symb,
                               STARPU_R, cellHandles[idxLevel][interactionid].up,
                               (STARPU_RW|STARPU_COMMUTE), cellHandles[idxLevel][interactionid].down,
                            0);
                }
            }
            FLOG( timerOutBlock.tac() );
        }
        FLOG( FLog::Controller << "\t\t transferPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock in  " << timerInBlock.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.elapsed() << "s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Downard Pass
    /////////////////////////////////////////////////////////////////////////////////////

    void downardPass(){
        FLOG( FTic timer; );
        for(int idxLevel = 2 ; idxLevel <= tree->getHeight()-2 ; ++idxLevel){
            /////////////////////////////////////////////////////////////
            // Exchange for MPI
            /////////////////////////////////////////////////////////////
            // Manage the external operations
            // Find what to recv
            if(tree->getNbCellGroupAtLevel(idxLevel)){
                // Take last block at this level
                const int idxLastBlock = tree->getNbCellGroupAtLevel(idxLevel)-1;
                const CellContainerClass* currentCells = tree->getCellGroup(idxLevel, idxLastBlock);
                // Take the last cell index of the last block
                const MortonIndex myLastIdx = currentCells->getEndingIndex()-1;
                // Find the descriptor of the first block that belong to someone else at lower level
                const int firstOtherBlock = nbBlocksBeforeMinPerLevel[idxLevel+1] + tree->getNbCellGroupAtLevel(idxLevel+1);
                FAssertLF(processesBlockInfos[idxLevel+1][firstOtherBlock-1].owner == comm.processId());
                // Iterate while the block has our cell has parent
                int idxBlockToSend = 0;
                int lastProcSend   = 0;
                while(firstOtherBlock + idxBlockToSend < int(processesBlockInfos[idxLevel+1].size()) &&
                      myLastIdx == (processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToSend].firstIndex >> 3)){

                    if(lastProcSend != processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToSend].owner){

                        FLOG(FLog::Controller << "[SMpi] Post a send during L2L for Idx " << tree->getCellGroup(idxLevel, idxLastBlock)->getStartingIndex() <<
                             " and dest is " << processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToSend].owner
                                << " size " << tree->getCellGroup(idxLevel, idxLastBlock)->getBufferSizeInByte()
                                << " tag " << getTag(idxLevel,tree->getCellGroup(idxLevel, idxLastBlock)->getStartingIndex(), 0) << "\n");
                        FLOG(FLog::Controller << "[SMpi] Post a send during L2L for Idx " << tree->getCellGroup(idxLevel, idxLastBlock)->getStartingIndex() <<
                             " and dest is " << processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToSend].owner
                                << " size " << tree->getCellGroup(idxLevel, idxLastBlock)->getLocalBufferSizeInByte()
                                << " tag " << getTag(idxLevel,tree->getCellGroup(idxLevel, idxLastBlock)->getStartingIndex(), 2) << "\n");

                        starpu_mpi_isend_detached( cellHandles[idxLevel][idxLastBlock].symb,
                                                   processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToSend].owner,
                                                   getTag(idxLevel,tree->getCellGroup(idxLevel, idxLastBlock)->getStartingIndex(), 0),
                                                   comm.getComm(), 0/*callback*/, 0/*arg*/ );
                        starpu_mpi_isend_detached( cellHandles[idxLevel][idxLastBlock].down,
                                                   processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToSend].owner,
                                                   getTag(idxLevel,tree->getCellGroup(idxLevel, idxLastBlock)->getStartingIndex(), 2),
                                                   comm.getComm(), 0/*callback*/, 0/*arg*/ );

                        lastProcSend = processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToSend].owner;
                    }
                    idxBlockToSend += 1;
                }
            }
            // Find what to send
            if(tree->getNbCellGroupAtLevel(idxLevel+1)
                    && nbBlocksBeforeMinPerLevel[idxLevel] != 0){
                // Take the first lower block
                const CellContainerClass* currentCells = tree->getCellGroup(idxLevel+1, 0);
                // Take its first index
                const MortonIndex myFirstChildIdx = currentCells->getStartingIndex();
                const MortonIndex missingParentIdx = (myFirstChildIdx>>3);
                // If no parent or the first parent is not the good one
                if(tree->getNbCellGroupAtLevel(idxLevel) == 0
                        || tree->getCellGroup(idxLevel, 0)->getStartingIndex() != missingParentIdx){

                    // Look if the parent is owned by another block
                    const int firstOtherBlock = nbBlocksBeforeMinPerLevel[idxLevel]-1;
                    FAssertLF(processesBlockInfos[idxLevel][firstOtherBlock].lastIndex-1 == missingParentIdx);

                    if(remoteCellGroups[idxLevel][firstOtherBlock].ptrSymb == nullptr){
                        const int nbBytesInBlock = processesBlockInfos[idxLevel][firstOtherBlock].bufferSizeSymb;
                        unsigned char* memoryBlock = (unsigned char*)FAlignedMemory::AllocateBytes<32>(nbBytesInBlock);
                        remoteCellGroups[idxLevel][firstOtherBlock].ptrSymb = memoryBlock;
                        starpu_variable_data_register(&remoteCellGroups[idxLevel][firstOtherBlock].handleSymb, 0,
                                                      (uintptr_t)remoteCellGroups[idxLevel][firstOtherBlock].ptrSymb, nbBytesInBlock);
                    }
                    if(remoteCellGroups[idxLevel][firstOtherBlock].ptrDown == nullptr){
                        const int nbBytesInBlock = processesBlockInfos[idxLevel][firstOtherBlock].bufferSizeDown;
                        unsigned char* memoryBlock = (unsigned char*)FAlignedMemory::AllocateBytes<32>(nbBytesInBlock);
                        remoteCellGroups[idxLevel][firstOtherBlock].ptrDown = memoryBlock;
                        starpu_variable_data_register(&remoteCellGroups[idxLevel][firstOtherBlock].handleDown, 0,
                                                      (uintptr_t)remoteCellGroups[idxLevel][firstOtherBlock].ptrDown, nbBytesInBlock);
                    }

                    FLOG(FLog::Controller << "[SMpi] Post a recv during L2L for Idx " << processesBlockInfos[idxLevel][firstOtherBlock].firstIndex <<
                         " and owner " << processesBlockInfos[idxLevel][firstOtherBlock].owner
                         << " size " << processesBlockInfos[idxLevel][firstOtherBlock].bufferSizeSymb
                         << " tag " << getTag(idxLevel,processesBlockInfos[idxLevel][firstOtherBlock].firstIndex, 0) << "\n");
                    FLOG(FLog::Controller << "[SMpi] Post a recv during L2L for Idx " << processesBlockInfos[idxLevel][firstOtherBlock].firstIndex <<
                         " and owner " << processesBlockInfos[idxLevel][firstOtherBlock].owner
                         << " size " << processesBlockInfos[idxLevel][firstOtherBlock].bufferSizeDown
                         << " tag " << getTag(idxLevel,processesBlockInfos[idxLevel][firstOtherBlock].firstIndex, 2) << "\n");

                    starpu_mpi_irecv_detached ( remoteCellGroups[idxLevel][firstOtherBlock].handleSymb,
                                                processesBlockInfos[idxLevel][firstOtherBlock].owner,
                                                getTag(idxLevel,processesBlockInfos[idxLevel][firstOtherBlock].firstIndex, 0),
                                                comm.getComm(), 0/*callback*/, 0/*arg*/ );
                    starpu_mpi_irecv_detached ( remoteCellGroups[idxLevel][firstOtherBlock].handleDown,
                                                processesBlockInfos[idxLevel][firstOtherBlock].owner,
                                                getTag(idxLevel,processesBlockInfos[idxLevel][firstOtherBlock].firstIndex, 2),
                                                comm.getComm(), 0/*callback*/, 0/*arg*/ );

                    {
                        struct starpu_task* const task = starpu_task_create();
                        task->dyn_handles = (starpu_data_handle_t*)malloc(sizeof(starpu_data_handle_t)*20);
                        task->dyn_handles[0] = remoteCellGroups[idxLevel][firstOtherBlock].handleSymb;
                        task->dyn_handles[1] = remoteCellGroups[idxLevel][firstOtherBlock].handleDown;

                        const MortonIndex parentStartingIdx = processesBlockInfos[idxLevel][firstOtherBlock].firstIndex;
                        const MortonIndex parentEndingIdx = processesBlockInfos[idxLevel][firstOtherBlock].lastIndex;

                        int idxSubGroup = 0;
                        // Skip current group if needed
                        if( tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (parentStartingIdx<<3) ){
                            ++idxSubGroup;
                            FAssertLF( idxSubGroup != tree->getNbCellGroupAtLevel(idxLevel+1) );
                            FAssertLF( (tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex()>>3) == parentStartingIdx );
                        }
                        // Copy at max 8 groups
                        int nbSubCellGroups = 0;
                        task->dyn_handles[(nbSubCellGroups*2) + 2] = cellHandles[idxLevel+1][idxSubGroup].symb;
                        task->dyn_handles[(nbSubCellGroups*2) + 3] = cellHandles[idxLevel+1][idxSubGroup].down;
                        nbSubCellGroups += 1;
                        while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= ((parentEndingIdx<<3)+7)
                              && (idxSubGroup + 1) != tree->getNbCellGroupAtLevel(idxLevel+1)
                              && tree->getCellGroup(idxLevel+1, idxSubGroup+1)->getStartingIndex() <= (parentEndingIdx<<3)+7 ){
                            idxSubGroup += 1;
                            task->dyn_handles[(nbSubCellGroups*2) + 2] = cellHandles[idxLevel+1][idxSubGroup].symb;
                            task->dyn_handles[(nbSubCellGroups*2) + 3] = cellHandles[idxLevel+1][idxSubGroup].down;
                            nbSubCellGroups += 1;
                            FAssertLF( nbSubCellGroups <= 9 );
                        }

                        // put the right codelet
                        task->cl = &l2l_cl[nbSubCellGroups-1];
                        // put args values
                        char *arg_buffer;
                        size_t arg_buffer_size;
                        starpu_codelet_pack_args((void**)&arg_buffer, &arg_buffer_size,
                                                 STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                                                 STARPU_VALUE, &nbSubCellGroups, sizeof(nbSubCellGroups),
                                                 STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                                                 0);
                        task->cl_arg = arg_buffer;
                        task->cl_arg_size = arg_buffer_size;
                        FAssertLF(starpu_task_submit(task) == 0);
                    }
                }
            }
            /////////////////////////////////////////////////////////////



            int idxSubGroup = 0;

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel)
                    && idxSubGroup < tree->getNbCellGroupAtLevel(idxLevel+1) ; ++idxGroup){
                CellContainerClass*const currentCells = tree->getCellGroup(idxLevel, idxGroup);

                struct starpu_task* const task = starpu_task_create();
                task->dyn_handles = (starpu_data_handle_t*)malloc(sizeof(starpu_data_handle_t)*20);
                task->dyn_handles[0] = cellHandles[idxLevel][idxGroup].symb;
                task->dyn_handles[1] = cellHandles[idxLevel][idxGroup].down;

                // Skip current group if needed
                if( tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++idxSubGroup;
                    FAssertLF( idxSubGroup != tree->getNbCellGroupAtLevel(idxLevel+1) );
                    FAssertLF( (tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                task->dyn_handles[(nbSubCellGroups*2) + 2] = cellHandles[idxLevel+1][idxSubGroup].symb;
                task->dyn_handles[(nbSubCellGroups*2) + 3] = cellHandles[idxLevel+1][idxSubGroup].down;
                nbSubCellGroups += 1;
                while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (((currentCells->getEndingIndex()-1)<<3)+7)
                      && (idxSubGroup+1) != tree->getNbCellGroupAtLevel(idxLevel+1)
                      && tree->getCellGroup(idxLevel+1, idxSubGroup+1)->getStartingIndex() <= ((currentCells->getEndingIndex()-1)<<3)+7 ){
                    idxSubGroup += 1;
                    task->dyn_handles[(nbSubCellGroups*2) + 2] = cellHandles[idxLevel+1][idxSubGroup].symb;
                    task->dyn_handles[(nbSubCellGroups*2) + 3] = cellHandles[idxLevel+1][idxSubGroup].down;
                    nbSubCellGroups += 1;
                    FAssertLF( nbSubCellGroups <= 9 );
                }

                // put the right codelet
                task->cl = &l2l_cl[nbSubCellGroups-1];
                // put args values
                char *arg_buffer;
                size_t arg_buffer_size;
                starpu_codelet_pack_args((void**)&arg_buffer, &arg_buffer_size,
                                         STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                                         STARPU_VALUE, &nbSubCellGroups, sizeof(nbSubCellGroups),
                                         STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                                         0);
                task->cl_arg = arg_buffer;
                task->cl_arg_size = arg_buffer_size;
                FAssertLF(starpu_task_submit(task) == 0);
            }
        }
        FLOG( FLog::Controller << "\t\t downardPass in " << timer.tacAndElapsed() << "s\n" );
    }
    /////////////////////////////////////////////////////////////////////////////////////
    /// Direct Pass MPI
    /////////////////////////////////////////////////////////////////////////////////////

    void directPassMpi(){
        FLOG( FTic timer; );

        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            for(int idxInteraction = 0; idxInteraction < int(externalInteractionsLeafLevelMpi[idxGroup].size()) ; ++idxInteraction){
                const int interactionid = externalInteractionsLeafLevelMpi[idxGroup][idxInteraction].otherBlockId;
                const std::vector<OutOfBlockInteraction>* outsideInteractions = &externalInteractionsLeafLevelMpi[idxGroup][idxInteraction].interactions;
                starpu_insert_task(&p2p_cl_inout_mpi,
                        STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                        STARPU_VALUE, &outsideInteractions, sizeof(outsideInteractions),
                        STARPU_R, particleHandles[idxGroup].symb,
                        (STARPU_RW|STARPU_COMMUTE), particleHandles[idxGroup].down,
                        STARPU_R, remoteParticleGroupss[interactionid].handleSymb,
                        0);
            }
        }

        FLOG( FLog::Controller << "\t\t directPass in MPI " << timer.tacAndElapsed() << "s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Direct Pass
    /////////////////////////////////////////////////////////////////////////////////////

    void directPass(){
        FLOG( FTic timer; );
        FLOG( FTic timerInBlock; FTic timerOutBlock; );

        FLOG( timerInBlock.tic() );
        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            starpu_insert_task(&p2p_cl_in,
                    STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                               STARPU_R, particleHandles[idxGroup].symb,
                               (STARPU_RW|STARPU_COMMUTE), particleHandles[idxGroup].down,
                    0);
        }
        FLOG( timerInBlock.tac() );
        FLOG( timerOutBlock.tic() );
        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            for(int idxInteraction = 0; idxInteraction < int(externalInteractionsLeafLevel[idxGroup].size()) ; ++idxInteraction){
                const int interactionid = externalInteractionsLeafLevel[idxGroup][idxInteraction].otherBlockId;
                const std::vector<OutOfBlockInteraction>* outsideInteractions = &externalInteractionsLeafLevel[idxGroup][idxInteraction].interactions;
                starpu_insert_task(&p2p_cl_inout,
                        STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                        STARPU_VALUE, &outsideInteractions, sizeof(outsideInteractions),
                        STARPU_R, particleHandles[idxGroup].symb,
                                   (STARPU_RW|STARPU_COMMUTE), particleHandles[idxGroup].down,
                        STARPU_R, particleHandles[interactionid].symb,
                                   (STARPU_RW|STARPU_COMMUTE), particleHandles[interactionid].down,
                        0);
            }
        }
        FLOG( timerOutBlock.tac() );

        FLOG( FLog::Controller << "\t\t directPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock  in " << timerInBlock.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.elapsed() << "s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Merge Pass
    /////////////////////////////////////////////////////////////////////////////////////

    void mergePass(){
        FLOG( FTic timer; );

        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            starpu_insert_task(&l2p_cl,
                    STARPU_VALUE, &wrapperptr, sizeof(wrapperptr),
                    STARPU_R, cellHandles[tree->getHeight()-1][idxGroup].symb,
                    STARPU_R, cellHandles[tree->getHeight()-1][idxGroup].down,
                    STARPU_R, particleHandles[idxGroup].symb,
                    (STARPU_RW|STARPU_COMMUTE), particleHandles[idxGroup].down,
                    0);
        }

        FLOG( FLog::Controller << "\t\t L2P in " << timer.tacAndElapsed() << "s\n" );
    }
};

#endif // FGROUPTASKSTARPUMPIALGORITHM_HPP
