#ifndef FGROUPTASKSTARPUALGORITHM_HPP
#define FGROUPTASKSTARPUALGORITHM_HPP

#include "../Utils/FGlobal.hpp"
#include "../Core/FCoreCommon.hpp"
#include "../Utils/FQuickSort.hpp"
#include "../Containers/FTreeCoordinate.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FTic.hpp"

#include <vector>
#include <vector>

#include <omp.h>

extern "C"{
#include <starpu.h>
}

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

    const int MaxThreads;         //< The number of threads
    OctreeClass*const tree;       //< The Tree
    KernelClass** kernels;        //< The kernels

    std::vector<starpu_data_handle_t>* handles;

    starpu_codelet p2m_cl;
    starpu_codelet m2m_cl[9];
    starpu_codelet l2l_cl[9];
    starpu_codelet l2p_cl;

    starpu_codelet m2l_cl_in;
    starpu_codelet m2l_cl_inout;

    starpu_codelet p2p_cl_in;
    starpu_codelet p2p_cl_inout;

public:
    FGroupTaskStarPUAlgorithm(OctreeClass*const inTree, KernelClass* inKernels, const int inMaxThreads = -1)
        : MaxThreads(inMaxThreads==-1?omp_get_max_threads():inMaxThreads), tree(inTree), kernels(nullptr),
          handles(nullptr){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");

        kernels = new KernelClass*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass(*inKernels);
        }

        handles = new std::vector<starpu_data_handle_t>[tree->getHeight()+1];

        FLOG(FLog::Controller << "FGroupTaskStarPUAlgorithm (Max Thread " << MaxThreads << ")\n");

        FAssertLF(starpu_init(NULL) == 0);
        starpu_pause();
    }

    ~FGroupTaskStarPUAlgorithm(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete[] kernels;

        cleanHandle();
        delete[] handles;

        starpu_shutdown();
    }

    void execute(const unsigned operationsToProceed = FFmmNearAndFarFields){
        FLOG( FLog::Controller << "\tStart FGroupTaskStarPUAlgorithm\n" );

        #pragma omp parallel
        #pragma omp single
        buildExternalInteractionVecs();

        buildHandles();

        starpu_resume();

        if( operationsToProceed & FFmmP2P ) directPass();

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if( operationsToProceed & FFmmL2P ) mergePass();

        starpu_task_wait_for_all();
        starpu_pause();
    }

protected:
    void initCodelet(){
        memset(&p2m_cl, 0, sizeof(p2m_cl));
        p2m_cl.where = STARPU_CPU;
        p2m_cl.cpu_funcs[0] = bottomPassCallback;
        p2m_cl.nbuffers = 2;
        p2m_cl.modes[0] = STARPU_RW;
        p2m_cl.modes[1] = STARPU_R;

        memset(m2m_cl, 0, sizeof(m2m_cl[0])*9);
        memset(l2l_cl, 0, sizeof(l2l_cl[0])*9);
        for(int idx = 0 ; idx < 9 ; ++idx){
            m2m_cl[idx].where = STARPU_CPU;
            m2m_cl[idx].cpu_funcs[0] = upwardPassCallback;
            m2m_cl[idx].nbuffers = idx+1;
            m2m_cl[idx].modes[0] = STARPU_RW;

            l2l_cl[idx].where = STARPU_CPU;
            l2l_cl[idx].cpu_funcs[0] = downardPassCallback;
            l2l_cl[idx].nbuffers = idx+1;
            l2l_cl[idx].modes[0] = STARPU_R;

            for(int idxBuffer = 0 ; idxBuffer <= idx ; ++idxBuffer){
                m2m_cl[idx].modes[idxBuffer+1] = STARPU_R;
                l2l_cl[idx].modes[idxBuffer+1] = STARPU_RW;
            }
        }

        memset(&l2p_cl, 0, sizeof(l2p_cl));
        l2p_cl.where = STARPU_CPU;
        l2p_cl.cpu_funcs[0] = mergePassCallback;
        l2p_cl.nbuffers = 2;
        l2p_cl.modes[0] = STARPU_R;
        l2p_cl.modes[1] = STARPU_RW;

        memset(&p2p_cl_in, 0, sizeof(p2p_cl_in));
        p2p_cl_in.where = STARPU_CPU;
        p2p_cl_in.cpu_funcs[0] = directInPassCallback;
        p2p_cl_in.nbuffers = 1;
        p2p_cl_in.modes[0] = STARPU_RW;
        memset(&p2p_cl_inout, 0, sizeof(p2p_cl_inout));
        p2p_cl_inout.where = STARPU_CPU;
        p2p_cl_inout.cpu_funcs[0] = directInoutPassCallback;
        p2p_cl_inout.nbuffers = 2;
        p2p_cl_inout.modes[0] = STARPU_RW;
        p2p_cl_inout.modes[1] = STARPU_RW;

        memset(&m2l_cl_in, 0, sizeof(m2l_cl_in));
        m2l_cl_in.where = STARPU_CPU;
        m2l_cl_in.cpu_funcs[0] = transferInPassCallback;
        m2l_cl_in.nbuffers = 1;
        m2l_cl_in.modes[0] = STARPU_RW;
        memset(&m2l_cl_inout, 0, sizeof(m2l_cl_inout));
        m2l_cl_inout.where = STARPU_CPU;
        m2l_cl_inout.cpu_funcs[0] = transferInoutPassCallback;
        m2l_cl_inout.nbuffers = 2;
        m2l_cl_inout.modes[0] = STARPU_RW;
        m2l_cl_inout.modes[1] = STARPU_RW;
    }

    /** dealloc in a starpu way all the defined handles */
    void cleanHandle(){
        for(int idxLevel = 0 ; idxLevel < tree->getHeight() ; ++idxLevel){
            for(int idxHandle = 0 ; idxHandle < int(handles[idxLevel].size()) ; ++idxHandle){
                starpu_data_unregister(handles[idxLevel][idxHandle]);
            }
            handles[idxLevel].clear();
        }
        {
            const int idxLevel = tree->getHeight();
            for(int idxHandle = 0 ; idxHandle < int(handles[idxLevel].size()) ; ++idxHandle){
                starpu_data_unregister(handles[idxLevel][idxHandle]);
            }
            handles[idxLevel].clear();
        }
    }

    /** Reset the handles array and create new ones to define
     * in a starpu way each block of data
     */
    void buildHandles(){
        cleanHandle();

        for(int idxLevel = 2 ; idxLevel < tree->getHeight() ; ++idxLevel){
            handles[idxLevel].resize(tree->getNbCellGroupAtLevel(idxLevel));

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                const CellContainerClass* currentCells = tree->getCellGroup(idxLevel, idxGroup);
                starpu_variable_data_register(&handles[idxLevel][idxGroup], 0, (uintptr_t)currentCells, 0); // TODO
            }
        }
        {
            const int idxLevel = tree->getHeight();
            for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
                ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);
                starpu_variable_data_register(&handles[idxLevel][idxGroup], 0, (uintptr_t)containers, 0); // TODO
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
                    STARPU_VALUE, this, sizeof(ThisClass*),
                    STARPU_RW, handles[tree->getHeight()-1][idxGroup],
                    STARPU_R, handles[tree->getHeight()][idxGroup],
                    0);
        }

        FLOG( FLog::Controller << "\t\t bottomPass in " << timer.tacAndElapsed() << "s\n" );
    }

    static void bottomPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass* leafCells  = STARPU_VARIABLE_GET_PTR(buffers[0]);
        ParticleGroupClass* containers = STARPU_VARIABLE_GET_PTR(buffers[1]);

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);
        worker->bottomPassPerform(leafCells, containers);
    }

    void bottomPassPerform(CellContainerClass* leafCells, ParticleGroupClass* containers){
        const MortonIndex blockStartIdx = leafCells->getStartingIndex();
        const MortonIndex blockEndIdx = leafCells->getEndingIndex();
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
            CellClass* cell = leafCells->getCell(mindex);
            if(cell){
                FAssertLF(cell->getMortonIndex() == mindex);
                ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(mindex);
                FAssertLF(particles.isAttachedToSomething());
                kernel->P2M(cell, &particles);
            }
        }
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
                task->handles[0] = handles[idxLevel][idxGroup];

                // Skip current group if needed
                if( tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++idxSubGroup;
                    FAssertLF( idxSubGroup != tree->getNbCellGroupAtLevel(idxLevel+1) );
                    FAssertLF( (tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                task->handles[nbSubCellGroups + 1] = handles[idxLevel+1][idxSubGroup];
                nbSubCellGroups += 1;
                while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= ((currentCells->getEndingIndex()<<3)+7)
                      && (++idxSubGroup) != tree->getNbCellGroupAtLevel(idxLevel+1)
                      && tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex() <= (currentCells->getEndingIndex()<<3)+7 ){
                    task->handles[nbSubCellGroups + 1] = handles[idxLevel+1][idxSubGroup];
                    nbSubCellGroups += 1;
                    FAssertLF( nbSubCellGroups <= 9 );
                }

                // put the right codelet
                task->cl = &m2m_cl[nbSubCellGroups-1];
                // put args values
                char *arg_buffer;
                size_t arg_buffer_size;
                starpu_codelet_pack_args((void**)&arg_buffer, &arg_buffer_size,
                                         STARPU_VALUE, this, sizeof(ThisClass*),
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

    static void upwardPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass*const currentCells = STARPU_VARIABLE_GET_PTR(buffers[0]);

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);
        int nbSubCellGroups = 0;
        starpu_codelet_unpack_args(cl_arg, &nbSubCellGroups);
        int idxLevel = 0;
        starpu_codelet_unpack_args(cl_arg, &idxLevel);

        CellContainerClass* subCellGroups[9];
        memset(subCellGroups, 0, 9*sizeof(CellContainerClass*));
        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            subCellGroups[idxSubGroup] = STARPU_VARIABLE_GET_PTR(buffers[idxSubGroup+1]);
        }

        worker->upwardPassPerform(currentCells, subCellGroups, nbSubCellGroups, idxLevel);
    }

    void upwardPassPerform(CellContainerClass*const currentCells,
                           CellContainerClass* subCellGroups[9],
                            const int nbSubCellGroups, const int idxLevel){
        const MortonIndex blockStartIdx = currentCells->getStartingIndex();
        const MortonIndex blockEndIdx   = currentCells->getEndingIndex();
        KernelClass*const kernel = kernels[starpu_worker_get_id()];
        int idxSubCellGroup = 0;

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
            CellClass* cell = currentCells->getCell(mindex);
            if(cell){
                FAssertLF(cell->getMortonIndex() == mindex);
                CellClass* child[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if( subCellGroups[idxSubCellGroup]->getEndingIndex() <= ((mindex<<3)+idxChild) ){
                        idxSubCellGroup += 1;
                    }
                    if( idxSubCellGroup == nbSubCellGroups ){
                        break;
                    }
                    child[idxChild] = subCellGroups[idxSubCellGroup]->getCell((mindex<<3)+idxChild);
                    FAssertLF(child[idxChild] == nullptr || child[idxChild]->getMortonIndex() == ((mindex<<3)+idxChild));
                }

                kernel->M2M(cell, child, idxLevel);
            }
        }
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
                        STARPU_VALUE, this, sizeof(ThisClass*),
                        STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                        STARPU_RW, handles[idxLevel][idxGroup],
                        0);
            }
            FLOG( timerInBlock.tac() );

            FLOG( timerOutBlock.tic() );
            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                for(int idxInteraction = 0; idxInteraction < int(externalInteractionsAllLevel[idxLevel][idxGroup].size()) ; ++idxInteraction){
                    const int interactionid = externalInteractionsAllLevel[idxLevel][idxGroup][idxInteraction].otherBlockId;
                    const std::vector<OutOfBlockInteraction>* outsideInteractions = &externalInteractionsAllLevel[idxLevel][idxGroup][idxInteraction].interactions;

                    starpu_insert_task(&m2l_cl_inout,
                            STARPU_VALUE, this, sizeof(ThisClass*),
                            STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                            STARPU_VALUE, outsideInteractions, sizeof(outsideInteractions),
                            STARPU_RW, handles[idxLevel][idxGroup],
                            STARPU_RW, handles[idxLevel][interactionid],
                            0);
                }
            }
            FLOG( timerOutBlock.tac() );
        }
        FLOG( FLog::Controller << "\t\t transferPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock in  " << timerInBlock.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.elapsed() << "s\n" );
    }

    static void transferInPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass*const currentCells = STARPU_VARIABLE_GET_PTR(buffers[0]);

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);
        int idxLevel = 0;
        starpu_codelet_unpack_args(cl_arg, &idxLevel);

        worker->transferInPassPerform(currentCells, idxLevel);
    }

    void transferInPassPerform(CellContainerClass*const currentCells, const int idxLevel){
        const MortonIndex blockStartIdx = currentCells->getStartingIndex();
        const MortonIndex blockEndIdx = currentCells->getEndingIndex();
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
            CellClass* cell = currentCells->getCell(mindex);
            if(cell){
                FAssertLF(cell->getMortonIndex() == mindex);
                MortonIndex interactionsIndexes[189];
                int interactionsPosition[189];
                const FTreeCoordinate coord(cell->getCoordinate());
                int counter = coord.getInteractionNeighbors(idxLevel,interactionsIndexes,interactionsPosition);

                const CellClass* interactions[343];
                memset(interactions, 0, 343*sizeof(CellClass*));
                int counterExistingCell = 0;

                for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                    if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                        CellClass* interCell = currentCells->getCell(interactionsIndexes[idxInter]);
                        if(interCell){
                            FAssertLF(interCell->getMortonIndex() == interactionsIndexes[idxInter]);
                            FAssertLF(interactions[interactionsPosition[idxInter]] == nullptr);
                            interactions[interactionsPosition[idxInter]] = interCell;
                            counterExistingCell += 1;
                        }
                    }
                }

                kernel->M2L( cell , interactions, counterExistingCell, idxLevel);
            }
        }
    }

    static void transferInoutPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass*const currentCells = STARPU_VARIABLE_GET_PTR(buffers[0]);
        CellContainerClass*const externalCells = STARPU_VARIABLE_GET_PTR(buffers[0]);

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);
        int idxLevel = 0;
        starpu_codelet_unpack_args(cl_arg, &idxLevel);
        const std::vector<OutOfBlockInteraction>* outsideInteractions;
        starpu_codelet_unpack_args(cl_arg, &outsideInteractions);

        worker->transferInoutPassPerform(currentCells, externalCells, idxLevel, outsideInteractions);
    }


    void transferInoutPassPerform(CellContainerClass*const currentCells,
                                  CellContainerClass*const cellsOther,
                                  const int idxLevel,
                                  const std::vector<OutOfBlockInteraction>* outsideInteractions){
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(int outInterIdx = 0 ; outInterIdx < int(outsideInteractions->size()) ; ++outInterIdx){
            CellClass* interCell = cellsOther->getCell((*outsideInteractions)[outInterIdx].outIndex);
            if(interCell){
                FAssertLF(interCell->getMortonIndex() == (*outsideInteractions)[outInterIdx].outIndex);
                CellClass* cell = currentCells->getCell((*outsideInteractions)[outInterIdx].insideIndex);
                FAssertLF(cell);
                FAssertLF(cell->getMortonIndex() == (*outsideInteractions)[outInterIdx].insideIndex);

                const CellClass* interactions[343];
                memset(interactions, 0, 343*sizeof(CellClass*));
                interactions[(*outsideInteractions)[outInterIdx].outPosition] = interCell;
                const int counter = 1;
                kernel->M2L( cell , interactions, counter, idxLevel);

                interactions[(*outsideInteractions)[outInterIdx].outPosition] = nullptr;
                interactions[getOppositeInterIndex((*outsideInteractions)[outInterIdx].outPosition)] = cell;
                kernel->M2L( interCell , interactions, counter, idxLevel);
            }
        }
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
                task->handles[0] = handles[idxLevel][idxGroup];

                // Skip current group if needed
                if( tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++idxSubGroup;
                    FAssertLF( idxSubGroup != tree->getNbCellGroupAtLevel(idxLevel+1) );
                    FAssertLF( (tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                task->handles[nbSubCellGroups + 1] = handles[idxLevel+1][idxSubGroup];
                nbSubCellGroups += 1;
                while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= ((currentCells->getEndingIndex()<<3)+7)
                      && (++idxSubGroup) != tree->getNbCellGroupAtLevel(idxLevel+1)
                      && tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex() <= (currentCells->getEndingIndex()<<3)+7 ){
                    task->handles[nbSubCellGroups + 1] = handles[idxLevel+1][idxSubGroup];
                    nbSubCellGroups += 1;
                    FAssertLF( nbSubCellGroups <= 9 );
                }

                // put the right codelet
                task->cl = &l2l_cl[nbSubCellGroups-1];
                // put args values
                char *arg_buffer;
                size_t arg_buffer_size;
                starpu_codelet_pack_args((void**)&arg_buffer, &arg_buffer_size,
                                         STARPU_VALUE, this, sizeof(ThisClass*),
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

    static void downardPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass*const currentCells = STARPU_VARIABLE_GET_PTR(buffers[0]);

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);
        int nbSubCellGroups = 0;
        starpu_codelet_unpack_args(cl_arg, &nbSubCellGroups);
        int idxLevel = 0;
        starpu_codelet_unpack_args(cl_arg, &idxLevel);

        CellContainerClass* subCellGroups[9];
        memset(subCellGroups, 0, 9*sizeof(CellContainerClass*));
        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            subCellGroups[idxSubGroup] = STARPU_VARIABLE_GET_PTR(buffers[idxSubGroup+1]);
        }

        worker->downardPassPerform(currentCells, subCellGroups, nbSubCellGroups, idxLevel);
    }

    void downardPassPerform(CellContainerClass*const currentCells,
                            CellContainerClass* subCellGroups[9],
                             const int nbSubCellGroups, const int idxLevel){
        const MortonIndex blockStartIdx = currentCells->getStartingIndex();
        const MortonIndex blockEndIdx = currentCells->getEndingIndex();
        KernelClass*const kernel = kernels[starpu_worker_get_id()];
        int idxSubCellGroup = 0;

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && idxSubCellGroup != nbSubCellGroups; ++mindex){
            CellClass* cell = currentCells->getCell(mindex);
            if(cell){
                FAssertLF(cell->getMortonIndex() == mindex);
                CellClass* child[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if( subCellGroups[idxSubCellGroup]->getEndingIndex() <= ((mindex<<3)+idxChild) ){
                        idxSubCellGroup += 1;
                    }
                    if( idxSubCellGroup == nbSubCellGroups ){
                        break;
                    }
                    child[idxChild] = subCellGroups[idxSubCellGroup]->getCell((mindex<<3)+idxChild);
                    FAssertLF(child[idxChild] == nullptr || child[idxChild]->getMortonIndex() == ((mindex<<3)+idxChild));
                }

                kernel->L2L(cell, child, idxLevel);
            }
        }
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
                    STARPU_VALUE, this, sizeof(ThisClass*),
                    STARPU_RW, handles[tree->getHeight()][idxGroup],
                    0);
        }
        FLOG( timerInBlock.tac() );
        FLOG( timerOutBlock.tic() );
        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            for(int idxInteraction = 0; idxInteraction < int(externalInteractionsLeafLevel.size()) ; ++idxInteraction){
                const int interactionid = externalInteractionsLeafLevel[idxGroup][idxInteraction].otherBlockId;
                const std::vector<OutOfBlockInteraction>* outsideInteractions = &externalInteractionsLeafLevel[idxGroup][idxInteraction].interactions;
                starpu_insert_task(&p2p_cl_in,
                        STARPU_VALUE, this, sizeof(ThisClass*),
                        STARPU_VALUE, outsideInteractions, sizeof(outsideInteractions),
                        STARPU_RW, handles[tree->getHeight()][idxGroup],
                        STARPU_RW, handles[tree->getHeight()][interactionid],
                        0);
            }
        }
        FLOG( timerOutBlock.tac() );

        FLOG( FLog::Controller << "\t\t directPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock  in " << timerInBlock.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.elapsed() << "s\n" );
    }

    static void directInPassCallback(void *buffers[], void *cl_arg){
        ParticleGroupClass* containers = STARPU_VARIABLE_GET_PTR(buffers[0]);

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);

        worker->directInPassPerform(containers);
    }

    void directInPassPerform(ParticleGroupClass* containers){
        const MortonIndex blockStartIdx = containers->getStartingIndex();
        const MortonIndex blockEndIdx = containers->getEndingIndex();
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
            ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(mindex);
            if(particles.isAttachedToSomething()){
                MortonIndex interactionsIndexes[26];
                int interactionsPosition[26];
                FTreeCoordinate coord(mindex, tree->getHeight()-1);
                int counter = coord.getNeighborsIndexes(tree->getHeight(),interactionsIndexes,interactionsPosition);

                ParticleContainerClass interactionsObjects[27];
                ParticleContainerClass* interactions[27];
                memset(interactions, 0, 27*sizeof(ParticleContainerClass*));
                int counterExistingCell = 0;

                for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                    if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                        interactionsObjects[counterExistingCell] = containers->template getLeaf<ParticleContainerClass>(interactionsIndexes[idxInter]);
                        if(interactionsObjects[counterExistingCell].isAttachedToSomething()){
                            FAssertLF(interactions[interactionsPosition[idxInter]] == nullptr);
                            interactions[interactionsPosition[idxInter]] = &interactionsObjects[counterExistingCell];
                            counterExistingCell += 1;
                        }
                    }
                }

                kernel->P2P( coord, &particles, &particles , interactions, counterExistingCell);
            }
        }
    }

    static void directInoutPassCallback(void *buffers[], void *cl_arg){
        ParticleGroupClass* containers = STARPU_VARIABLE_GET_PTR(buffers[0]);
        ParticleGroupClass* externalContainers = STARPU_VARIABLE_GET_PTR(buffers[1]);

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);
        const std::vector<OutOfBlockInteraction>* outsideInteractions = nullptr;
        starpu_codelet_unpack_args(cl_arg, &outsideInteractions);

        worker->directInPassPerform(containers, externalContainers, outsideInteractions);
    }

    void directInoutPassPerform(ParticleGroupClass* containers, ParticleGroupClass* containersOther,
                                const std::vector<OutOfBlockInteraction>* outsideInteractions){
        KernelClass*const kernel = kernels[starpu_worker_get_id()];
        for(int outInterIdx = 0 ; outInterIdx < int(outsideInteractions->size()) ; ++outInterIdx){
            ParticleContainerClass interParticles = containersOther->template getLeaf<ParticleContainerClass>((*outsideInteractions)[outInterIdx].outIndex);
            if(interParticles.isAttachedToSomething()){
                ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>((*outsideInteractions)[outInterIdx].insideIndex);
                FAssertLF(particles.isAttachedToSomething());
                ParticleContainerClass* interactions[27];
                memset(interactions, 0, 27*sizeof(ParticleContainerClass*));
                interactions[(*outsideInteractions)[outInterIdx].outPosition] = &interParticles;
                const int counter = 1;
                kernel->P2PRemote( FTreeCoordinate((*outsideInteractions)[outInterIdx].insideIndex, tree->getHeight()-1), &particles, &particles , interactions, counter);

                interactions[(*outsideInteractions)[outInterIdx].outPosition] = nullptr;
                interactions[getOppositeNeighIndex((*outsideInteractions)[outInterIdx].outPosition)] = &particles;
                kernel->P2PRemote( FTreeCoordinate((*outsideInteractions)[outInterIdx].outIndex, tree->getHeight()-1), &interParticles, &interParticles , interactions, counter);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////
    /// Merge Pass
    /////////////////////////////////////////////////////////////////////////////////////

    void mergePass(){
        FLOG( FTic timer; );

        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            starpu_insert_task(&l2p_cl,
                    STARPU_VALUE, this, sizeof(ThisClass*),
                    STARPU_R, handles[tree->getHeight()-1][idxGroup],
                    STARPU_RW, handles[tree->getHeight()][idxGroup],
                    0);
        }

        FLOG( FLog::Controller << "\t\t L2P in " << timer.tacAndElapsed() << "s\n" );
    }

    static void mergePassCallback(void *buffers[], void *cl_arg){
        CellContainerClass* leafCells  = STARPU_VARIABLE_GET_PTR(buffers[0]);
        ParticleGroupClass* containers = STARPU_VARIABLE_GET_PTR(buffers[1]);

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);
        worker->mergePassPerform(leafCells, containers);
    }

    void mergePassPerform(CellContainerClass* leafCells, ParticleGroupClass* containers){
        const MortonIndex blockStartIdx = leafCells->getStartingIndex();
        const MortonIndex blockEndIdx = leafCells->getEndingIndex();
        KernelClass*const kernel = kernels[starpu_worker_get_id()];

        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
            CellClass* cell = leafCells->getCell(mindex);
            if(cell){
                FAssertLF(cell->getMortonIndex() == mindex);
                ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(mindex);
                FAssertLF(particles.isAttachedToSomething());
                kernel->L2P(cell, &particles);
            }
        }
    }


    /////////////////////////////////////////////////////////////////////////////////////
    /// Utils Pass
    /////////////////////////////////////////////////////////////////////////////////////

    int getOppositeNeighIndex(const int index) const {
        // ((idxX+1)*3 + (idxY+1)) * 3 + (idxZ+1)
        return 27-index-1;
    }

    int getOppositeInterIndex(const int index) const {
        // ((( (xdiff+3) * 7) + (ydiff+3))) * 7 + zdiff + 3
        return 343-index-1;
    }
};

#endif // FGROUPTASKSTARPUALGORITHM_HPP
