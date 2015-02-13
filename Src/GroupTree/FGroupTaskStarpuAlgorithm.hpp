
// Keep in private GIT
// @SCALFMM_PRIVATE
#ifndef FGROUPTASKSTARPUALGORITHM_HPP
#define FGROUPTASKSTARPUALGORITHM_HPP

#include "../Utils/FGlobal.hpp"
#include "../Core/FCoreCommon.hpp"
#include "../Utils/FQuickSort.hpp"
#include "../Containers/FTreeCoordinate.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FAssert.hpp"

#include <vector>
#include <vector>

#include <omp.h>

extern "C"{
#include <starpu.h>
}

#include "FStarPUCpuWrapper.hpp"
#include "FStarPUUtils.hpp"

template <class OctreeClass, class CellContainerClass, class CellClass, class KernelClass, class ParticleGroupClass, class ParticleContainerClass>
class FGroupTaskStarPUAlgorithm {
protected:
    typedef FGroupTaskStarPUAlgorithm<OctreeClass, CellContainerClass, CellClass, KernelClass, ParticleGroupClass, ParticleContainerClass> ThisClass;

    struct OutOfBlockInteraction{
        MortonIndex outIndex;
        MortonIndex insideIndex;
        int outPosition;
        // To sort
        bool operator <=(const OutOfBlockInteraction& other) const{
            return outIndex <= other.outIndex;
        }
    };

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

    std::vector<starpu_data_handle_t>* handles_up;
    std::vector<starpu_data_handle_t>* handles_down;

    starpu_codelet p2m_cl;
    starpu_codelet m2m_cl[9];
    starpu_codelet l2l_cl[9];
    starpu_codelet l2p_cl;

    starpu_codelet m2l_cl_in;
    starpu_codelet m2l_cl_inout;

    starpu_codelet p2p_cl_in;
    starpu_codelet p2p_cl_inout;

    typedef FStarPUCpuWrapper<CellContainerClass, CellClass, KernelClass, ParticleGroupClass, ParticleContainerClass> StarPUCpuWrapperClass;
    StarPUCpuWrapperClass cpuWrapper;

    FStarPUPtrInterface wrappers;
    FStarPUPtrInterface* wrapperptr;

public:
    FGroupTaskStarPUAlgorithm(OctreeClass*const inTree, KernelClass* inKernels, const int inMaxThreads = -1)
        : MaxThreads(inMaxThreads), tree(inTree),
          handles_up(nullptr), handles_down(nullptr),
            cpuWrapper(tree->getHeight()), wrapperptr(&wrappers){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");
        FAssertLF(MaxThreads <= STARPU_MAXCPUS, "number of threads to high");

        struct starpu_conf conf;
        FAssertLF(starpu_conf_init(&conf) == 0);
        conf.ncpus = MaxThreads;
        FAssertLF(starpu_init(&conf) == 0);

        starpu_pthread_mutex_t initMutex;
        starpu_pthread_mutex_init(&initMutex, NULL);
        FStarPUUtils::ExecOnWorkers(STARPU_CPU, [&](){
            starpu_pthread_mutex_lock(&initMutex);
            cpuWrapper.initKernel(starpu_worker_get_id(), inKernels);
            starpu_pthread_mutex_unlock(&initMutex);
        });
        wrappers.set(FSTARPU_CPU_IDX, &cpuWrapper);

        starpu_pthread_mutex_destroy(&initMutex);

        starpu_pause();

        MaxThreads = starpu_worker_get_count();//starpu_cpu_worker_get_count();

        handles_up   = new std::vector<starpu_data_handle_t>[tree->getHeight()+1];
        handles_down = new std::vector<starpu_data_handle_t>[tree->getHeight()+1];

        initCodelet();

        FLOG(FLog::Controller << "FGroupTaskStarPUAlgorithm (Max Thread " << MaxThreads << ")\n");
    }

    ~FGroupTaskStarPUAlgorithm(){
        cleanHandle();
        delete[] handles_up;
        delete[] handles_down;

        starpu_resume();
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
        p2m_cl.where = STARPU_CPU;
        p2m_cl.cpu_funcs[0] = StarPUCpuWrapperClass::bottomPassCallback;
        p2m_cl.nbuffers = 2;
        p2m_cl.modes[0] = STARPU_RW;
        p2m_cl.modes[1] = STARPU_R;
        p2m_cl.name = "p2m_cl";

        memset(m2m_cl, 0, sizeof(m2m_cl[0])*9);
        memset(l2l_cl, 0, sizeof(l2l_cl[0])*9);
        for(int idx = 0 ; idx < 9 ; ++idx){
            m2m_cl[idx].where = STARPU_CPU;
            m2m_cl[idx].cpu_funcs[0] = StarPUCpuWrapperClass::upwardPassCallback;
            m2m_cl[idx].nbuffers = idx+2;
            m2m_cl[idx].dyn_modes = (starpu_data_access_mode*)malloc((idx+2)*sizeof(starpu_data_access_mode));
            m2m_cl[idx].dyn_modes[0] = STARPU_RW;
            m2m_cl[idx].name = "m2m_cl";

            l2l_cl[idx].where = STARPU_CPU;
            l2l_cl[idx].cpu_funcs[0] = StarPUCpuWrapperClass::downardPassCallback;
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
        l2p_cl.where = STARPU_CPU;
        l2p_cl.cpu_funcs[0] = StarPUCpuWrapperClass::mergePassCallback;
        l2p_cl.nbuffers = 2;
        l2p_cl.modes[0] = STARPU_R;
        l2p_cl.modes[1] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        l2p_cl.name = "l2p_cl";

        memset(&p2p_cl_in, 0, sizeof(p2p_cl_in));
        p2p_cl_in.where = STARPU_CPU;
        p2p_cl_in.cpu_funcs[0] = StarPUCpuWrapperClass::directInPassCallback;
        p2p_cl_in.nbuffers = 1;
        p2p_cl_in.modes[0] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        p2p_cl_in.name = "p2p_cl_in";
        memset(&p2p_cl_inout, 0, sizeof(p2p_cl_inout));
        p2p_cl_inout.where = STARPU_CPU;
        p2p_cl_inout.cpu_funcs[0] = StarPUCpuWrapperClass::directInoutPassCallback;
        p2p_cl_inout.nbuffers = 2;
        p2p_cl_inout.modes[0] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        p2p_cl_inout.modes[1] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        p2p_cl_inout.name = "p2p_cl_inout";

        memset(&m2l_cl_in, 0, sizeof(m2l_cl_in));
        m2l_cl_in.where = STARPU_CPU;
        m2l_cl_in.cpu_funcs[0] = StarPUCpuWrapperClass::transferInPassCallback;
        m2l_cl_in.nbuffers = 2;
        m2l_cl_in.modes[0] = starpu_data_access_mode(STARPU_RW|STARPU_COMMUTE);
        m2l_cl_in.modes[1] = STARPU_R;
        m2l_cl_in.name = "m2l_cl_in";
        memset(&m2l_cl_inout, 0, sizeof(m2l_cl_inout));
        m2l_cl_inout.where = STARPU_CPU;
        m2l_cl_inout.cpu_funcs[0] = StarPUCpuWrapperClass::transferInoutPassCallback;
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

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
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
            int idxSubGroup = 0;

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
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

#endif // FGROUPTASKSTARPUALGORITHM_HPP
