// Keep in private GIT
// @SCALFMM_PRIVATE
#ifndef FGROUPTASKSTARPUMPIALGORITHM_HPP
#define FGROUPTASKSTARPUMPIALGORITHM_HPP

#include "../Utils/FGlobal.hpp"
#include "../Core/FCoreCommon.hpp"
#include "../Utils/FQuickSort.hpp"
#include "../Containers/FTreeCoordinate.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FAssert.hpp"
#include "../Utils/FAlignedMemory.hpp"
#include "../Utils/FAssert.hpp"

#include "../Utils/FMpi.hpp"

#include "FOutOfBlockInteraction.hpp"

#include <vector>
#include <memory>

#include <omp.h>

//extern "C"{
#include <starpu.h>
#include <starpu_mpi.h>
//}

#ifdef STARPU_USE_CPU
#include "FStarPUCpuWrapper.hpp"
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
#include "FStarPUCudaWrapper.hpp"
#include "Cuda/FCudaEmptyKernel.hpp"
#include "Cuda/FCudaGroupAttachedLeaf.hpp"
#include "Cuda/FCudaGroupOfParticles.hpp"
#include "Cuda/FCudaGroupOfCells.hpp"
#endif
#ifdef STARPU_USE_OPENCL
#include "FStarPUOpenClWrapper.hpp"
#endif
#include "FStarPUUtils.hpp"

template <class OctreeClass, class CellContainerClass, class CellClass, class KernelClass, class ParticleGroupClass, class ParticleContainerClass
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
    , class CudaCellContainerClass = FCudaGroupOfCells<0>, class CudaParticleGroupClass = FCudaGroupOfParticles<0, int>, class CudaParticleContainerClass = FCudaGroupAttachedLeaf<0, int>,
    class CudaKernelClass = FCudaEmptyKernel<>
#endif
          >
class FGroupTaskStarPUMpiAlgorithm {
protected:
    typedef FGroupTaskStarPUMpiAlgorithm<OctreeClass, CellContainerClass, CellClass, KernelClass, ParticleGroupClass, ParticleContainerClass
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        , CudaCellContainerClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass
#endif
    > ThisClass;

    int getTag(const int inLevel, const MortonIndex mindex) const{
        int shift = 0;
        int height = tree->getHeight();
        while(height) { shift += 1; height >>= 1; }
        return (mindex<<shift) + inLevel;
    }

    const FMpi::FComm& comm;

    template <class OtherBlockClass>
    struct BlockInteractions{
        OtherBlockClass* otherBlock;
        int otherBlockId;
        std::vector<OutOfBlockInteraction> interactions;
    };

    std::vector< std::vector< std::vector<BlockInteractions<CellContainerClass>>>> externalInteractionsAllLevel;
    std::vector< std::vector<BlockInteractions<ParticleGroupClass>>> externalInteractionsLeafLevel;

    int MaxThreads;         //< The number of threads
    OctreeClass*const tree;       //< The Tree
    KernelClass*const originalCpuKernel;

    std::vector<starpu_data_handle_t>* handles_up;
    std::vector<starpu_data_handle_t>* handles_down;

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
    typedef FStarPUCpuWrapper<CellContainerClass, CellClass, KernelClass, ParticleGroupClass, ParticleContainerClass> StarPUCpuWrapperClass;
    StarPUCpuWrapperClass cpuWrapper;
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
    typedef FStarPUCudaWrapper<KernelClass, CudaCellContainerClass, CudaParticleGroupClass, CudaParticleContainerClass, CudaKernelClass> StarPUCudaWrapperClass;
    StarPUCudaWrapperClass cudaWrapper;
#endif
#ifdef STARPU_USE_OPENCL
    typedef FStarPUOpenClWrapper<CellContainerClass, CellClass, KernelClass, ParticleGroupClass, ParticleContainerClass> StarPUOpenClWrapperClass;
    StarPUOpenClWrapperClass openclWrapper;
#endif

    FStarPUPtrInterface wrappers;
    FStarPUPtrInterface* wrapperptr;

public:
    FGroupTaskStarPUMpiAlgorithm(const FMpi::FComm& inComm, OctreeClass*const inTree, KernelClass* inKernels, const int inMaxThreads = -1)
        :   comm(inComm), MaxThreads(inMaxThreads), tree(inTree), originalCpuKernel(inKernels),
            handles_up(nullptr), handles_down(nullptr),
#ifdef STARPU_USE_CPU
            cpuWrapper(tree->getHeight()),
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
            cudaWrapper(tree->getHeight()),
#endif
#ifdef STARPU_USE_OPENCL
            openclWrapper(tree->getHeight()),
#endif
            wrapperptr(&wrappers){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");
        FAssertLF(MaxThreads <= STARPU_MAXCPUS, "number of threads to high");

        struct starpu_conf conf;
        FAssertLF(starpu_conf_init(&conf) == 0);
        conf.ncpus = MaxThreads;
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
#ifdef STARPU_USE_OPENCL
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

        handles_up = new std::vector<starpu_data_handle_t>[tree->getHeight()+1];
        handles_down = new std::vector<starpu_data_handle_t>[tree->getHeight()+1];

        initCodelet();
        initCodeletMpi();

        FLOG(FLog::Controller << "FGroupTaskStarPUMpiAlgorithm (Max Thread " << MaxThreads << ")\n");
    }

    ~FGroupTaskStarPUMpiAlgorithm(){
        cleanHandle();
        cleanHandleMpi();
        delete[] handles_up;
        delete[] handles_down;

        starpu_resume();
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
        postRecvAllocatedBlocks();

        starpu_resume();

        if( operationsToProceed & FFmmP2P ) insertParticlesSend();
        if(operationsToProceed & FFmmM2L) insertCellsSend();

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

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
#ifdef STARPU_USE_OPENCL
        if(originalCpuKernel->supportP2M(FSTARPU_OPENCL_IDX)){
            p2m_cl.opencl_funcs[0] = StarPUOpenClWrapperClass::bottomPassCallback;
            p2m_cl.where |= STARPU_OPENCL;
        }
#endif
        p2m_cl.nbuffers = 2;
        p2m_cl.modes[0] = STARPU_RW;
        p2m_cl.modes[1] = STARPU_R;
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
#ifdef STARPU_USE_OPENCL
            if(originalCpuKernel->supportM2M(FSTARPU_OPENCL_IDX)){
                m2m_cl[idx].opencl_funcs[0] = StarPUOpenClWrapperClass::upwardPassCallback;
                m2m_cl[idx].where |= STARPU_OPENCL;
            }
#endif
            m2m_cl[idx].nbuffers = idx+2;
            m2m_cl[idx].dyn_modes = (starpu_data_access_mode*)malloc((idx+2)*sizeof(starpu_data_access_mode));
            m2m_cl[idx].dyn_modes[0] = STARPU_RW;
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
#ifdef STARPU_USE_OPENCL
            if(originalCpuKernel->supportL2L(FSTARPU_OPENCL_IDX)){
                l2l_cl[idx].opencl_funcs[0] = StarPUOpenClWrapperClass::downardPassCallback;
                l2l_cl[idx].where |= STARPU_OPENCL;
            }
#endif
            l2l_cl[idx].nbuffers = idx+2;
            l2l_cl[idx].dyn_modes = (starpu_data_access_mode*)malloc((idx+2)*sizeof(starpu_data_access_mode));
            l2l_cl[idx].dyn_modes[0] = STARPU_R;
            l2l_cl[idx].name = "l2l_cl";

            for(int idxBuffer = 0 ; idxBuffer <= idx ; ++idxBuffer){
                m2m_cl[idx].dyn_modes[idxBuffer+1] = STARPU_R;
                l2l_cl[idx].dyn_modes[idxBuffer+1] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
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
#ifdef STARPU_USE_OPENCL
        if(originalCpuKernel->supportL2P(FSTARPU_OPENCL_IDX)){
            l2p_cl.opencl_funcs[0] = StarPUOpenClWrapperClass::mergePassCallback;
            l2p_cl.where |= STARPU_OPENCL;
        }
#endif
        l2p_cl.nbuffers = 2;
        l2p_cl.modes[0] = STARPU_R;
        l2p_cl.modes[1] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
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
#ifdef STARPU_USE_OPENCL
        if(originalCpuKernel->supportP2P(FSTARPU_OPENCL_IDX)){
            p2p_cl_in.opencl_funcs[0] = StarPUOpenClWrapperClass::directInPassCallback;
            p2p_cl_in.where |= STARPU_OPENCL;
        }
#endif
        p2p_cl_in.nbuffers = 1;
        p2p_cl_in.modes[0] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        p2p_cl_in.name = "p2p_cl_in";
        memset(&p2p_cl_inout, 0, sizeof(p2p_cl_inout));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportP2P(FSTARPU_CPU_IDX)){
            p2p_cl_inout.cpu_funcs[0] = StarPUCpuWrapperClass::directInoutPassCallback;
            p2p_cl_inout.where |= STARPU_CPU;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportP2P(FSTARPU_CUDA_IDX)){
            p2p_cl_inout.cuda_funcs[0] = StarPUCudaWrapperClass::directInoutPassCallback;
            p2p_cl_inout.where |= STARPU_CUDA;
        }
#endif
#ifdef STARPU_USE_OPENCL
        if(originalCpuKernel->supportP2P(FSTARPU_OPENCL_IDX)){
            p2p_cl_inout.opencl_funcs[0] = StarPUOpenClWrapperClass::directInoutPassCallback;
            p2p_cl_inout.where |= STARPU_OPENCL;
        }
#endif
        p2p_cl_inout.nbuffers = 2;
        p2p_cl_inout.modes[0] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        p2p_cl_inout.modes[1] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
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
#ifdef STARPU_USE_OPENCL
        if(originalCpuKernel->supportM2L(FSTARPU_OPENCL_IDX)){
            m2l_cl_in.opencl_funcs[0] = StarPUOpenClWrapperClass::transferInPassCallback;
            m2l_cl_in.where |= STARPU_OPENCL;
        }
#endif
        m2l_cl_in.nbuffers = 2;
        m2l_cl_in.modes[0] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        m2l_cl_in.modes[1] = STARPU_R;
        m2l_cl_in.name = "m2l_cl_in";
        memset(&m2l_cl_inout, 0, sizeof(m2l_cl_inout));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportM2L(FSTARPU_CPU_IDX)){
            m2l_cl_inout.cpu_funcs[0] = StarPUCpuWrapperClass::transferInoutPassCallback;
            m2l_cl_inout.where |= STARPU_CPU;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportM2L(FSTARPU_CUDA_IDX)){
            m2l_cl_inout.cuda_funcs[0] = StarPUCudaWrapperClass::transferInoutPassCallback;
            m2l_cl_inout.where |= STARPU_CUDA;
        }
#endif
#ifdef STARPU_USE_OPENCL
        if(originalCpuKernel->supportM2L(FSTARPU_OPENCL_IDX)){
            m2l_cl_inout.opencl_funcs[0] = StarPUOpenClWrapperClass::transferInoutPassCallback;
            m2l_cl_inout.where |= STARPU_OPENCL;
        }
#endif
        m2l_cl_inout.nbuffers = 4;
        m2l_cl_inout.modes[0] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        m2l_cl_inout.modes[1] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        m2l_cl_inout.modes[2] = STARPU_R;
        m2l_cl_inout.modes[3] = STARPU_R;
        m2l_cl_inout.name = "m2l_cl_inout";
    }

    /** dealloc in a starpu way all the defined handles */
    void cleanHandle(){
        for(int idxLevel = 0 ; idxLevel < tree->getHeight() ; ++idxLevel){
            for(int idxHandle = 0 ; idxHandle < int(handles_up[idxLevel].size()) ; ++idxHandle){
                starpu_data_unregister(handles_up[idxLevel][idxHandle]);
            }
            handles_up[idxLevel].clear();
            for(int idxHandle = 0 ; idxHandle < int(handles_down[idxLevel].size()) ; ++idxHandle){
                starpu_data_unregister(handles_down[idxLevel][idxHandle]);
            }
            handles_down[idxLevel].clear();
        }
        {
            const int idxLevel = tree->getHeight();
            for(int idxHandle = 0 ; idxHandle < int(handles_up[idxLevel].size()) ; ++idxHandle){
                starpu_data_unregister(handles_up[idxLevel][idxHandle]);
            }
            handles_up[idxLevel].clear();
            for(int idxHandle = 0 ; idxHandle < int(handles_down[idxLevel].size()) ; ++idxHandle){
                starpu_data_unregister(handles_down[idxLevel][idxHandle]);
            }
            handles_down[idxLevel].clear();
        }
    }

    ////////////////////////////////////////////////////////////////////////////

    void initCodeletMpi(){
        memset(&p2p_cl_inout_mpi, 0, sizeof(p2p_cl_inout_mpi));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportM2L(FSTARPU_CPU_IDX)){
            p2p_cl_inout_mpi.where |= STARPU_CPU;
            p2p_cl_inout_mpi.cpu_funcs[0] = StarPUCpuWrapperClass::directInoutPassCallbackMpi;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportM2L(FSTARPU_CUDA_IDX)){
            p2p_cl_inout_mpi.where |= STARPU_CUDA;
            p2p_cl_inout_mpi.cuda_funcs[0] = StarPUCudaWrapperClass::directInoutPassCallbackMpi;
        }
#endif
#ifdef STARPU_USE_OPENCL
        if(originalCpuKernel->supportM2L(FSTARPU_OPENCL_IDX)){
            p2p_cl_inout_mpi.where |= STARPU_OPENCL;
            p2p_cl_inout_mpi.opencl_funcs[0] = StarPUOpenClWrapperClass::directInoutPassCallbackMpi;
        }
#endif
        p2p_cl_inout_mpi.nbuffers = 2;
        p2p_cl_inout_mpi.modes[0] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        p2p_cl_inout_mpi.modes[1] = STARPU_R;
        p2p_cl_inout_mpi.name = "p2p_cl_inout_mpi";

        memset(&m2l_cl_inout_mpi, 0, sizeof(m2l_cl_inout_mpi));
#ifdef STARPU_USE_CPU
        if(originalCpuKernel->supportM2L(FSTARPU_CPU_IDX)){
            m2l_cl_inout_mpi.where |= STARPU_CPU;
            m2l_cl_inout_mpi.cpu_funcs[0] = StarPUCpuWrapperClass::transferInoutPassCallbackMpi;
        }
#endif
#ifdef ScalFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportM2L(FSTARPU_CUDA_IDX)){
            m2l_cl_inout_mpi.where |= STARPU_CUDA;
            m2l_cl_inout_mpi.cuda_funcs[0] = StarPUCudaWrapperClass::transferInoutPassCallbackMpi;
        }
#endif
#ifdef STARPU_USE_OPENCL
        if(originalCpuKernel->supportM2L(FSTARPU_OPENCL_IDX)){
            m2l_cl_inout_mpi.where |= STARPU_OPENCL;
            m2l_cl_inout_mpi.opencl_funcs[0] = StarPUOpenClWrapperClass::transferInoutPassCallbackMpi;
        }
#endif
        m2l_cl_inout_mpi.nbuffers = 2;
        m2l_cl_inout_mpi.modes[0] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        m2l_cl_inout_mpi.modes[1] = STARPU_R;
        m2l_cl_inout_mpi.name = "m2l_cl_inout_mpi";
    }

    std::vector<std::pair<MortonIndex,MortonIndex>> processesIntervalPerLevels;
    struct BlockDescriptor{
        MortonIndex firstIndex;
        MortonIndex lastIndex;
        int owner;
        int bufferSize;
        int leavesBufferSize;
    };
    std::vector<std::vector<BlockDescriptor>> processesBlockInfos;
    std::vector<int> nbBlocksPerLevelAll;
    std::vector<int> nbBlocksBeforeMinPerLevel;

    std::vector< std::vector< std::vector<BlockInteractions<CellContainerClass>>>> externalInteractionsAllLevelMpi;
    std::vector< std::vector<BlockInteractions<ParticleGroupClass>>> externalInteractionsLeafLevelMpi;

    struct RemoteHandle{
        RemoteHandle() : ptr(nullptr){
            memset(&handle, 0, sizeof(handle));
        }

        unsigned char * ptr;
        starpu_data_handle_t handle;
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
                myBlocksAtLevel[idxGroup].bufferSize = currentCells->getBufferSizeInByte();

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
                const CellContainerClass* currentCells = tree->getCellGroup(idxLevel, idxGroup);

                std::vector<BlockInteractions<CellContainerClass>>* externalInteractions = &externalInteractionsAllLevelMpi[idxLevel][idxGroup];

                #pragma omp task default(none) firstprivate(idxGroup, currentCells, idxLevel, externalInteractions)
                {
                    std::vector<OutOfBlockInteraction> outsideInteractions;
                    const MortonIndex blockStartIdx = currentCells->getStartingIndex();
                    const MortonIndex blockEndIdx   = currentCells->getEndingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                        const CellClass* cell = currentCells->getCell(mindex);
                        if(cell){
                            FAssertLF(cell->getMortonIndex() == mindex);
                            MortonIndex interactionsIndexes[189];
                            int interactionsPosition[189];
                            const FTreeCoordinate coord(cell->getCoordinate());
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
                            if(remoteCellGroups[idxLevel][idxOtherGroup].ptr == nullptr){
                                #pragma omp critical(CreateM2LRemotes)
                                {
                                    if(remoteCellGroups[idxLevel][idxOtherGroup].ptr == nullptr){
                                        const int nbBytesInBlock = processesBlockInfos[idxLevel][idxOtherGroup].bufferSize;
                                        unsigned char* memoryBlock = (unsigned char*)FAlignedMemory::Allocate32BAligned(nbBytesInBlock);
                                        remoteCellGroups[idxLevel][idxOtherGroup].ptr = memoryBlock;
                                        starpu_variable_data_register(&remoteCellGroups[idxLevel][idxOtherGroup].handle, 0,
                                                                      (uintptr_t)remoteCellGroups[idxLevel][idxOtherGroup].ptr, nbBytesInBlock);
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
                        ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(mindex);
                        if(particles.isAttachedToSomething()){
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
                            if(remoteParticleGroupss[idxOtherGroup].ptr == nullptr){
                                #pragma omp critical(CreateM2LRemotes)
                                {
                                    if(remoteParticleGroupss[idxOtherGroup].ptr == nullptr){
                                        const int nbBytesInBlock = processesBlockInfos[tree->getHeight()-1][idxOtherGroup].leavesBufferSize;
                                        unsigned char* memoryBlock = (unsigned char*)FAlignedMemory::Allocate32BAligned(nbBytesInBlock);
                                        remoteParticleGroupss[idxOtherGroup].ptr = memoryBlock;
                                        starpu_variable_data_register(&remoteParticleGroupss[idxOtherGroup].handle, 0,
                                                                      (uintptr_t)remoteParticleGroupss[idxOtherGroup].ptr, nbBytesInBlock);
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
                if(remoteCellGroups[idxLevel][idxHandle].ptr){
                    FLOG(FLog::Controller << "[SMpi] Post a recv during M2L for Idx " << processesBlockInfos[idxLevel][idxHandle].firstIndex <<
                         " and dest is " << processesBlockInfos[idxLevel][idxHandle].owner << "\n");

                    starpu_mpi_irecv_detached( remoteCellGroups[idxLevel][idxHandle].handle,
                                                processesBlockInfos[idxLevel][idxHandle].owner,
                                                getTag(idxLevel,processesBlockInfos[idxLevel][idxHandle].firstIndex),
                                                comm.getComm(), 0, 0 );

                    toRecv.push_back({processesBlockInfos[idxLevel][idxHandle].owner,
                                        comm.processId(), idxLevel, idxHandle});
                }
            }
        }
        {
            for(int idxHandle = 0 ; idxHandle < int(remoteParticleGroupss.size()) ; ++idxHandle){
                if(remoteParticleGroupss[idxHandle].ptr){
                    FLOG(FLog::Controller << "[SMpi] Post a recv during P2P for Idx " << processesBlockInfos[tree->getHeight()-1][idxHandle].firstIndex <<
                         " and dest is " << processesBlockInfos[tree->getHeight()-1][idxHandle].owner << "\n");

                    starpu_mpi_irecv_detached( remoteParticleGroupss[idxHandle].handle,
                                                processesBlockInfos[tree->getHeight()-1][idxHandle].owner,
                                                getTag(tree->getHeight(),processesBlockInfos[tree->getHeight()-1][idxHandle].firstIndex),
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
                     " and dest is " << sd.dest << "\n");

                starpu_mpi_isend_detached( handles_up[tree->getHeight()][localId], sd.dest,
                        getTag(tree->getHeight(),tree->getParticleGroup(localId)->getStartingIndex()),
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

                FLOG(FLog::Controller << "[SMpi] Post a send during M2L for Idx " << tree->getCellGroup(sd.level, localId)->getStartingIndex() <<
                     " and dest is " << sd.dest << "\n");

                starpu_mpi_isend_detached( handles_up[sd.level][localId], sd.dest,
                        getTag(sd.level,tree->getCellGroup(sd.level, localId)->getStartingIndex()),
                        comm.getComm(), 0/*callback*/, 0/*arg*/ );
            }
        }
    }

    void cleanHandleMpi(){
        for(int idxLevel = 0 ; idxLevel < int(remoteCellGroups.size()) ; ++idxLevel){
            for(int idxHandle = 0 ; idxHandle < int(remoteCellGroups[idxLevel].size()) ; ++idxHandle){
                if(remoteCellGroups[idxLevel][idxHandle].ptr){
                    starpu_data_unregister(remoteCellGroups[idxLevel][idxHandle].handle);
                    FAlignedMemory::Dealloc32BAligned(remoteCellGroups[idxLevel][idxHandle].ptr);
                }
            }
            remoteCellGroups[idxLevel].clear();
        }
        {
            for(int idxHandle = 0 ; idxHandle < int(remoteParticleGroupss.size()) ; ++idxHandle){
                if(remoteParticleGroupss[idxHandle].ptr){
                    starpu_data_unregister(remoteParticleGroupss[idxHandle].handle);
                    FAlignedMemory::Dealloc32BAligned(remoteParticleGroupss[idxHandle].ptr);
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
            handles_up[idxLevel].resize(tree->getNbCellGroupAtLevel(idxLevel));
            handles_down[idxLevel].resize(tree->getNbCellGroupAtLevel(idxLevel));

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                const CellContainerClass* currentCells = tree->getCellGroup(idxLevel, idxGroup);
                starpu_variable_data_register(&handles_up[idxLevel][idxGroup], 0,
                                              (uintptr_t)currentCells->getRawBuffer(), currentCells->getBufferSizeInByte());
                starpu_variable_data_register(&handles_down[idxLevel][idxGroup], 0,
                                              (uintptr_t)currentCells->getRawBuffer(), currentCells->getBufferSizeInByte());
            }
        }
        {
            const int idxLevel = tree->getHeight();
            handles_up[idxLevel].resize(tree->getNbParticleGroup());
            handles_down[idxLevel].resize(tree->getNbParticleGroup());
            for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
                ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);
                starpu_variable_data_register(&handles_up[idxLevel][idxGroup], 0,
                                              (uintptr_t)containers->getRawBuffer(), containers->getBufferSizeInByte());
                starpu_variable_data_register(&handles_down[idxLevel][idxGroup], 0,
                                              (uintptr_t)containers->getRawBuffer(), containers->getBufferSizeInByte());
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
                        ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(mindex);
                        if(particles.isAttachedToSomething()){
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
                    const CellContainerClass* currentCells = tree->getCellGroup(idxLevel, idxGroup);

                    std::vector<BlockInteractions<CellContainerClass>>* externalInteractions = &externalInteractionsAllLevel[idxLevel][idxGroup];

                    #pragma omp task default(none) firstprivate(idxGroup, currentCells, idxLevel, externalInteractions)
                    {
                        std::vector<OutOfBlockInteraction> outsideInteractions;
                        const MortonIndex blockStartIdx = currentCells->getStartingIndex();
                        const MortonIndex blockEndIdx   = currentCells->getEndingIndex();

                        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                            const CellClass* cell = currentCells->getCell(mindex);
                            if(cell){
                                FAssertLF(cell->getMortonIndex() == mindex);
                                MortonIndex interactionsIndexes[189];
                                int interactionsPosition[189];
                                const FTreeCoordinate coord(cell->getCoordinate());
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
                    STARPU_RW, handles_up[tree->getHeight()-1][idxGroup],
                    STARPU_R, handles_up[tree->getHeight()][idxGroup],
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
                task->dyn_handles = (starpu_data_handle_t*)malloc(sizeof(starpu_data_handle_t)*10);
                task->dyn_handles[0] = handles_up[idxLevel][idxGroup];

                // Skip current group if needed
                if( tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++idxSubGroup;
                    FAssertLF( idxSubGroup != tree->getNbCellGroupAtLevel(idxLevel+1) );
                    FAssertLF( (tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                task->dyn_handles[nbSubCellGroups + 1] = handles_up[idxLevel+1][idxSubGroup];
                nbSubCellGroups += 1;
                while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (((currentCells->getEndingIndex()-1)<<3)+7)
                      && (idxSubGroup+1) != tree->getNbCellGroupAtLevel(idxLevel+1)
                      && tree->getCellGroup(idxLevel+1, idxSubGroup+1)->getStartingIndex() <= ((currentCells->getEndingIndex()-1)<<3)+7 ){
                    idxSubGroup += 1;
                    task->dyn_handles[nbSubCellGroups + 1] = handles_up[idxLevel+1][idxSubGroup];
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

                    if(remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].ptr == nullptr){
                        const int nbBytesInBlock = processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].bufferSize;
                        unsigned char* memoryBlock = (unsigned char*)FAlignedMemory::Allocate32BAligned(nbBytesInBlock);
                        remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].ptr = memoryBlock;
                        starpu_variable_data_register(&remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].handle, 0,
                                                      (uintptr_t)remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].ptr, nbBytesInBlock);
                    }

                    FLOG(FLog::Controller << "[SMpi] Post a recv during M2M for Idx " << processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].firstIndex <<
                         " and owner is " << processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].owner << "\n");

                    starpu_mpi_irecv_detached ( remoteCellGroups[idxLevel+1][firstOtherBlock + idxBlockToRecv].handle,
                                                processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].owner,
                                                getTag(idxLevel,processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToRecv].firstIndex),
                                                comm.getComm(), 0/*callback*/, 0/*arg*/ );


                    idxBlockToRecv += 1;
                }
                FAssertLF(idxBlockToRecv < 8);
                if(idxBlockToRecv){// Perform the work
                    struct starpu_task* const task = starpu_task_create();
                    task->dyn_handles = (starpu_data_handle_t*)malloc(sizeof(starpu_data_handle_t)*10);
                    task->dyn_handles[0] = handles_up[idxLevel][tree->getNbCellGroupAtLevel(idxLevel)-1];

                    // Copy at max 8 groups
                    int nbSubCellGroups = 0;
                    while(nbSubCellGroups < idxBlockToRecv){
                        task->dyn_handles[nbSubCellGroups + 1] = remoteCellGroups[idxLevel+1][firstOtherBlock + nbSubCellGroups].handle;
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
                             " and dest is " << dest << "\n");

                        starpu_mpi_isend_detached( handles_up[idxLevel+1][lowerIdxToSend], dest,
                                getTag(idxLevel,tree->getCellGroup(idxLevel+1, lowerIdxToSend)->getStartingIndex()),
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
                            (STARPU_RW|STARPU_COMMUTE), handles_down[idxLevel][idxGroup],
                            STARPU_R, remoteCellGroups[idxLevel][interactionid].handle,
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
                                   (STARPU_RW|STARPU_COMMUTE), handles_down[idxLevel][idxGroup],
                                   STARPU_R, handles_up[idxLevel][idxGroup],
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
                            (STARPU_RW|STARPU_COMMUTE), handles_down[idxLevel][idxGroup],
                            (STARPU_RW|STARPU_COMMUTE), handles_down[idxLevel][interactionid],
                                       STARPU_R, handles_up[idxLevel][idxGroup],
                                       STARPU_R, handles_up[idxLevel][interactionid],
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
                             " and dest is " << processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToSend].owner << "\n");

                        starpu_mpi_isend_detached( handles_down[idxLevel][idxLastBlock],
                                                   processesBlockInfos[idxLevel+1][firstOtherBlock + idxBlockToSend].owner,
                                                   getTag(idxLevel,tree->getCellGroup(idxLevel, idxLastBlock)->getStartingIndex()),
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

                    if(remoteCellGroups[idxLevel][firstOtherBlock].ptr == nullptr){
                        const int nbBytesInBlock = processesBlockInfos[idxLevel][firstOtherBlock].bufferSize;
                        unsigned char* memoryBlock = (unsigned char*)FAlignedMemory::Allocate32BAligned(nbBytesInBlock);
                        remoteCellGroups[idxLevel][firstOtherBlock].ptr = memoryBlock;
                        starpu_variable_data_register(&remoteCellGroups[idxLevel][firstOtherBlock].handle, 0,
                                                      (uintptr_t)remoteCellGroups[idxLevel][firstOtherBlock].ptr, nbBytesInBlock);
                    }

                    FLOG(FLog::Controller << "[SMpi] Post a recv during L2L for Idx " << processesBlockInfos[idxLevel][firstOtherBlock].firstIndex <<
                         " and owner " << processesBlockInfos[idxLevel][firstOtherBlock].owner << "\n");

                    starpu_mpi_irecv_detached ( remoteCellGroups[idxLevel][firstOtherBlock].handle,
                                                processesBlockInfos[idxLevel][firstOtherBlock].owner,
                                                getTag(idxLevel,processesBlockInfos[idxLevel][firstOtherBlock].firstIndex),
                                                comm.getComm(), 0/*callback*/, 0/*arg*/ );

                    {
                        struct starpu_task* const task = starpu_task_create();
                        task->dyn_handles = (starpu_data_handle_t*)malloc(sizeof(starpu_data_handle_t)*10);
                        task->dyn_handles[0] = remoteCellGroups[idxLevel][firstOtherBlock].handle;

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
                        task->dyn_handles[nbSubCellGroups + 1] = handles_down[idxLevel+1][idxSubGroup];
                        nbSubCellGroups += 1;
                        while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= ((parentEndingIdx<<3)+7)
                              && (idxSubGroup + 1) != tree->getNbCellGroupAtLevel(idxLevel+1)
                              && tree->getCellGroup(idxLevel+1, idxSubGroup+1)->getStartingIndex() <= (parentEndingIdx<<3)+7 ){
                            idxSubGroup += 1;
                            task->dyn_handles[nbSubCellGroups + 1] = handles_down[idxLevel+1][idxSubGroup];
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
                task->dyn_handles = (starpu_data_handle_t*)malloc(sizeof(starpu_data_handle_t)*10);
                task->dyn_handles[0] = handles_down[idxLevel][idxGroup];

                // Skip current group if needed
                if( tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++idxSubGroup;
                    FAssertLF( idxSubGroup != tree->getNbCellGroupAtLevel(idxLevel+1) );
                    FAssertLF( (tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                task->dyn_handles[nbSubCellGroups + 1] = handles_down[idxLevel+1][idxSubGroup];
                nbSubCellGroups += 1;
                while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (((currentCells->getEndingIndex()-1)<<3)+7)
                      && (idxSubGroup+1) != tree->getNbCellGroupAtLevel(idxLevel+1)
                      && tree->getCellGroup(idxLevel+1, idxSubGroup+1)->getStartingIndex() <= ((currentCells->getEndingIndex()-1)<<3)+7 ){
                    idxSubGroup += 1;
                    task->dyn_handles[nbSubCellGroups + 1] = handles_down[idxLevel+1][idxSubGroup];
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
                        (STARPU_RW|STARPU_COMMUTE), handles_down[tree->getHeight()][idxGroup],
                        STARPU_R, remoteParticleGroupss[interactionid].handle,
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
                    (STARPU_RW|STARPU_COMMUTE), handles_down[tree->getHeight()][idxGroup],
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
                        (STARPU_RW|STARPU_COMMUTE), handles_down[tree->getHeight()][idxGroup],
                        (STARPU_RW|STARPU_COMMUTE), handles_down[tree->getHeight()][interactionid],
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
                    STARPU_R, handles_down[tree->getHeight()-1][idxGroup],
                    (STARPU_RW|STARPU_COMMUTE), handles_down[tree->getHeight()][idxGroup],
                    0);
        }

        FLOG( FLog::Controller << "\t\t L2P in " << timer.tacAndElapsed() << "s\n" );
    }
};

#endif // FGROUPTASKSTARPUMPIALGORITHM_HPP
