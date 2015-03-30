// Keep in private GIT
// @SCALFMM_PRIVATE
#ifndef FGROUPTASKSTARPUALGORITHM_HPP
#define FGROUPTASKSTARPUALGORITHM_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Core/FCoreCommon.hpp"
#include "../../Utils/FQuickSort.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../../Utils/FLog.hpp"
#include "../../Utils/FTic.hpp"
#include "../../Utils/FAssert.hpp"

#include "FOutOfBlockInteraction.hpp"

#include <vector>
#include <vector>

#include <omp.h>

#include <starpu.h>
#include "../StarPUUtils/FStarPUUtils.hpp"
#include "../StarPUUtils/FStarPUFmmPriorities.hpp"

#ifdef STARPU_USE_CPU
#include "../StarPUUtils/FStarPUCpuWrapper.hpp"
#endif
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
#include "../StarPUUtils/FStarPUCudaWrapper.hpp"
#include "../Cuda/FCudaEmptyKernel.hpp"
#include "../Cuda/FCudaGroupAttachedLeaf.hpp"
#include "../Cuda/FCudaGroupOfParticles.hpp"
#include "../Cuda/FCudaGroupOfCells.hpp"
#include "../Cuda/FCudaEmptyCellSymb.hpp"
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
#include "../StarPUUtils/FStarPUOpenClWrapper.hpp"
#include "../OpenCl/FOpenCLDeviceWrapper.hpp"
#endif


template <class OctreeClass, class CellContainerClass, class KernelClass, class ParticleGroupClass, class StarPUCpuWrapperClass
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
    , class StarPUCudaWrapperClass = FStarPUCudaWrapper<KernelClass, FCudaEmptyCellSymb, int, int, FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                                                        FCudaGroupOfParticles<int, 0, 0, int>, FCudaGroupAttachedLeaf<int, 0, 0, int>, FCudaEmptyKernel<int> >
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
    , class StarPUOpenClWrapperClass = FStarPUOpenClWrapper<KernelClass, FOpenCLDeviceWrapper<KernelClass>>
#endif
          >
class FGroupTaskStarPUAlgorithm {
protected:
    typedef FGroupTaskStarPUAlgorithm<OctreeClass, CellContainerClass, KernelClass, ParticleGroupClass, StarPUCpuWrapperClass
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
        , StarPUCudaWrapperClass
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
    , StarPUOpenClWrapperClass
#endif
    > ThisClass;

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
        int intervalSize;
    };

    struct ParticleHandles{
        starpu_data_handle_t symb;
        starpu_data_handle_t down;
        int intervalSize;
    };

    std::vector< std::vector< std::vector<BlockInteractions<CellContainerClass>>>> externalInteractionsAllLevel;
    std::vector< std::vector<BlockInteractions<ParticleGroupClass>>> externalInteractionsLeafLevel;

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

    starpu_codelet p2p_cl_in;
    starpu_codelet p2p_cl_inout;

#ifdef STARPU_USE_CPU
    StarPUCpuWrapperClass cpuWrapper;
#endif
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
    StarPUCudaWrapperClass cudaWrapper;
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
    StarPUOpenClWrapperClass openclWrapper;
#endif

    FStarPUPtrInterface wrappers;
    FStarPUPtrInterface* wrapperptr;

public:
    FGroupTaskStarPUAlgorithm(OctreeClass*const inTree, KernelClass* inKernels)
        : tree(inTree), originalCpuKernel(inKernels),
          cellHandles(nullptr),
#ifdef STARPU_USE_CPU
            cpuWrapper(tree->getHeight()),
#endif
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
            cudaWrapper(tree->getHeight()),
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
            openclWrapper(tree->getHeight()),
#endif
            wrapperptr(&wrappers){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");

        struct starpu_conf conf;
        FAssertLF(starpu_conf_init(&conf) == 0);
        FStarPUFmmPriorities::Controller().init(&conf, tree->getHeight(), inKernels);
        FAssertLF(starpu_init(&conf) == 0);

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
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
        FStarPUUtils::ExecOnWorkers(STARPU_CUDA, [&](){
            starpu_pthread_mutex_lock(&initMutex);
            cudaWrapper.initKernel(starpu_worker_get_id(), inKernels);
            starpu_pthread_mutex_unlock(&initMutex);
        });
        wrappers.set(FSTARPU_CUDA_IDX, &cudaWrapper);
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
        FStarPUUtils::ExecOnWorkers(STARPU_OPENCL, [&](){
            starpu_pthread_mutex_lock(&initMutex);
            openclWrapper.initKernel(starpu_worker_get_id(), inKernels);
            starpu_pthread_mutex_unlock(&initMutex);
        });
        wrappers.set(FSTARPU_OPENCL_IDX, &openclWrapper);
#endif
        starpu_pthread_mutex_destroy(&initMutex);

        starpu_pause();

        cellHandles   = new std::vector<CellHandles>[tree->getHeight()];

        initCodelet();

        FLOG(FLog::Controller << "FGroupTaskStarPUAlgorithm (Max Worker " << starpu_worker_get_count() << ")\n");
#ifdef STARPU_USE_CPU
        FLOG(FLog::Controller << "FGroupTaskStarPUAlgorithm (Max CPU " << starpu_cpu_worker_get_count() << ")\n");
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
        FLOG(FLog::Controller << "FGroupTaskStarPUAlgorithm (Max OpenCL " << starpu_opencl_worker_get_count() << ")\n");
#endif
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
        FLOG(FLog::Controller << "FGroupTaskStarPUAlgorithm (Max CUDA " << starpu_cuda_worker_get_count() << ")\n");
#endif
    }

    ~FGroupTaskStarPUAlgorithm(){
        starpu_resume();

        cleanHandle();
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
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
        FStarPUUtils::ExecOnWorkers(STARPU_CUDA, [&](){
            starpu_pthread_mutex_lock(&releaseMutex);
            cudaWrapper.releaseKernel(starpu_worker_get_id());
            starpu_pthread_mutex_unlock(&releaseMutex);
        });
        wrappers.set(FSTARPU_CUDA_IDX, &cudaWrapper);
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
        FStarPUUtils::ExecOnWorkers(STARPU_OPENCL, [&](){
            starpu_pthread_mutex_lock(&releaseMutex);
            openclWrapper.releaseKernel(starpu_worker_get_id());
            starpu_pthread_mutex_unlock(&releaseMutex);
        });
        wrappers.set(FSTARPU_OPENCL_IDX, &openclWrapper);
#endif
        starpu_pthread_mutex_destroy(&releaseMutex);

        starpu_shutdown();
    }

    void execute(const unsigned operationsToProceed = FFmmNearAndFarFields){
        FLOG( FLog::Controller << "\tStart FGroupTaskStarPUAlgorithm\n" );

        #pragma omp parallel
        #pragma omp single
        buildExternalInteractionVecs();

        buildHandles();

        starpu_resume();

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if( operationsToProceed & FFmmP2P ) directPass();

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
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportP2M(FSTARPU_CUDA_IDX)){
            p2m_cl.cuda_funcs[0] = StarPUCudaWrapperClass::bottomPassCallback;
            p2m_cl.where |= STARPU_CUDA;
        }
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
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
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
            if(originalCpuKernel->supportM2M(FSTARPU_CUDA_IDX)){
                m2m_cl[idx].cuda_funcs[0] = StarPUCudaWrapperClass::upwardPassCallback;
                m2m_cl[idx].where |= STARPU_CUDA;
            }
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
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
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
            if(originalCpuKernel->supportL2L(FSTARPU_CUDA_IDX)){
                l2l_cl[idx].cuda_funcs[0] = StarPUCudaWrapperClass::downardPassCallback;
                l2l_cl[idx].where |= STARPU_CUDA;
            }
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
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
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportL2P(FSTARPU_CUDA_IDX)){
            l2p_cl.cuda_funcs[0] = StarPUCudaWrapperClass::mergePassCallback;
            l2p_cl.where |= STARPU_CUDA;
        }
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
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
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportP2P(FSTARPU_CUDA_IDX)){
            p2p_cl_in.cuda_funcs[0] = StarPUCudaWrapperClass::directInPassCallback;
            p2p_cl_in.where |= STARPU_CUDA;
        }
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
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
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportP2PExtern(FSTARPU_CUDA_IDX)){
            p2p_cl_inout.cuda_funcs[0] = StarPUCudaWrapperClass::directInoutPassCallback;
            p2p_cl_inout.where |= STARPU_CUDA;
        }
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
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
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportM2L(FSTARPU_CUDA_IDX)){
            m2l_cl_in.cuda_funcs[0] = StarPUCudaWrapperClass::transferInPassCallback;
            m2l_cl_in.where |= STARPU_CUDA;
        }
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
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
#ifdef SCALFMM_ENABLE_CUDA_KERNEL
        if(originalCpuKernel->supportM2LExtern(FSTARPU_CUDA_IDX)){
            m2l_cl_inout.cuda_funcs[0] = StarPUCudaWrapperClass::transferInoutPassCallback;
            m2l_cl_inout.where |= STARPU_CUDA;
        }
#endif
#ifdef SCALFMM_ENABLE_OPENCL_KERNEL
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
                cellHandles[idxLevel][idxGroup].intervalSize = int(currentCells->getEndingIndex() - currentCells->getStartingIndex());
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
                particleHandles[idxGroup].intervalSize = int(containers->getEndingIndex() - containers->getStartingIndex());
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
                        const MortonIndex blockStartIdxOther    = leftContainers->getStartingIndex();
                        const MortonIndex blockEndIdxOther      = leftContainers->getEndingIndex();

                        while(currentOutInteraction < int(outsideInteractions.size()) && outsideInteractions[currentOutInteraction].outIndex < blockStartIdxOther){
                            currentOutInteraction += 1;
                        }

                        int lastOutInteraction = currentOutInteraction;
                        while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdxOther){
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
                            const MortonIndex blockStartIdxOther = leftCells->getStartingIndex();
                            const MortonIndex blockEndIdxOther   = leftCells->getEndingIndex();

                            while(currentOutInteraction < int(outsideInteractions.size()) && outsideInteractions[currentOutInteraction].outIndex < blockStartIdxOther){
                                currentOutInteraction += 1;
                            }

                            int lastOutInteraction = currentOutInteraction;
                            while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdxOther){
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
                    STARPU_VALUE, &cellHandles[tree->getHeight()-1][idxGroup].intervalSize, sizeof(int),
                    STARPU_PRIORITY, FStarPUFmmPriorities::Controller().getPrioP2M(),
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

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
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
                                         STARPU_VALUE, &cellHandles[idxLevel][idxGroup].intervalSize, sizeof(int),
                                         0);
                task->cl_arg = arg_buffer;
                task->cl_arg_size = arg_buffer_size;
                task->priority = FStarPUFmmPriorities::Controller().getPrioM2M(idxLevel);
                FAssertLF(starpu_task_submit(task) == 0);
            }
        }
        FLOG( FLog::Controller << "\t\t upwardPass in " << timer.tacAndElapsed() << "s\n" );
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
                        STARPU_VALUE, &cellHandles[idxLevel][idxGroup].intervalSize, sizeof(int),
                                   STARPU_PRIORITY, FStarPUFmmPriorities::Controller().getPrioM2L(idxLevel),
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
                            STARPU_VALUE, &cellHandles[idxLevel][idxGroup].intervalSize, sizeof(int),
                                       STARPU_PRIORITY, FStarPUFmmPriorities::Controller().getPrioM2LExtern(idxLevel),
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
            int idxSubGroup = 0;

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
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
                                         STARPU_VALUE, &cellHandles[idxLevel][idxGroup].intervalSize, sizeof(int),
                                         0);
                task->cl_arg = arg_buffer;
                task->cl_arg_size = arg_buffer_size;
                task->priority = FStarPUFmmPriorities::Controller().getPrioL2L(idxLevel);
                FAssertLF(starpu_task_submit(task) == 0);
            }
        }
        FLOG( FLog::Controller << "\t\t downardPass in " << timer.tacAndElapsed() << "s\n" );
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
                    STARPU_VALUE, &particleHandles[idxGroup].intervalSize, sizeof(int),
                               STARPU_PRIORITY, FStarPUFmmPriorities::Controller().getPrioP2P(),
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
                        STARPU_VALUE, &particleHandles[idxGroup].intervalSize, sizeof(int),
                                   STARPU_PRIORITY, FStarPUFmmPriorities::Controller().getPrioP2PExtern(),
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
                    STARPU_VALUE, &cellHandles[tree->getHeight()-1][idxGroup].intervalSize, sizeof(int),
                               STARPU_PRIORITY, FStarPUFmmPriorities::Controller().getPrioL2P(),
                    STARPU_R, cellHandles[tree->getHeight()-1][idxGroup].symb,
                    STARPU_R, cellHandles[tree->getHeight()-1][idxGroup].down,
                    STARPU_R, particleHandles[idxGroup].symb,
                    (STARPU_RW|STARPU_COMMUTE), particleHandles[idxGroup].down,
                    0);
        }

        FLOG( FLog::Controller << "\t\t L2P in " << timer.tacAndElapsed() << "s\n" );
    }
};

#endif // FGROUPTASKSTARPUALGORITHM_HPP
