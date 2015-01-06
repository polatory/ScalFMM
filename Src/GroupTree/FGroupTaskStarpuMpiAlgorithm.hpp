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

#include "../Utils/FMpi.hpp"

#include <vector>
#include <memory>

#include <omp.h>

extern "C"{
#include <starpu.h>
#include <starpu_mpi.h>
}

template <class OctreeClass, class CellContainerClass, class CellClass, class KernelClass, class ParticleGroupClass, class ParticleContainerClass>
class FGroupTaskStarPUMpiAlgorithm {
protected:
    typedef FGroupTaskStarPUMpiAlgorithm<OctreeClass, CellContainerClass, CellClass, KernelClass, ParticleGroupClass, ParticleContainerClass> ThisClass;

    int getTag(const int inLevel, const MortonIndex mindex) const{
        int shift = 0;
        int height = tree->getHeight();
        while(height) { shift += 1; height >>= 1; }
        return (mindex<<shift) + inLevel;
    }

    const FMpi::FComm& comm;

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
    KernelClass** kernels;        //< The kernels
    ThisClass* thisptr;

    std::vector<starpu_data_handle_t>* handles;

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

public:
    FGroupTaskStarPUMpiAlgorithm(const FMpi::FComm& inComm, OctreeClass*const inTree, KernelClass* inKernels, const int inMaxThreads = -1)
        : comm(inComm), MaxThreads(inMaxThreads), tree(inTree), kernels(nullptr),
          thisptr(this), handles(nullptr){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");
        FAssertLF(MaxThreads <= STARPU_MAXCPUS, "number of threads to high");

        struct starpu_conf conf;
        FAssertLF(starpu_conf_init(&conf) == 0);
        conf.ncpus = MaxThreads;
        FAssertLF(starpu_init(&conf) == 0);
        FAssertLF(starpu_mpi_init ( 0, 0, 0 ) == 0);
        starpu_pause();

        MaxThreads = starpu_worker_get_count();//starpu_cpu_worker_get_count();

        handles = new std::vector<starpu_data_handle_t>[tree->getHeight()+1];
        kernels = new KernelClass*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass(*inKernels);
        }

        initCodelet();
        initCodeletMpi();

        FLOG(FLog::Controller << "FGroupTaskStarPUMpiAlgorithm (Max Thread " << MaxThreads << ")\n");
    }

    ~FGroupTaskStarPUMpiAlgorithm(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete[] kernels;

        cleanHandle();
        cleanHandleMpi();
        delete[] handles;

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
        if( operationsToProceed & FFmmP2P ) directPass();
        if( operationsToProceed & FFmmP2P ) directPassMpi();

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();
        if(operationsToProceed & FFmmM2L) insertCellsSend();

        if(operationsToProceed & FFmmM2L) transferPass();
        if(operationsToProceed & FFmmM2L) transferPassMpi();

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
            m2m_cl[idx].nbuffers = idx+2;
            m2m_cl[idx].dyn_modes = (starpu_data_access_mode*)malloc((idx+2)*sizeof(starpu_data_access_mode));
            m2m_cl[idx].dyn_modes[0] = STARPU_RW;

            l2l_cl[idx].where = STARPU_CPU;
            l2l_cl[idx].cpu_funcs[0] = downardPassCallback;
            l2l_cl[idx].nbuffers = idx+2;
            l2l_cl[idx].dyn_modes = (starpu_data_access_mode*)malloc((idx+2)*sizeof(starpu_data_access_mode));
            l2l_cl[idx].dyn_modes[0] = STARPU_R;

            for(int idxBuffer = 0 ; idxBuffer <= idx ; ++idxBuffer){
                m2m_cl[idx].dyn_modes[idxBuffer+1] = STARPU_R;
                l2l_cl[idx].dyn_modes[idxBuffer+1] = STARPU_RW;
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

    ////////////////////////////////////////////////////////////////////////////

    void initCodeletMpi(){
        memset(&p2p_cl_inout_mpi, 0, sizeof(p2p_cl_inout_mpi));
        p2p_cl_inout_mpi.where = STARPU_CPU;
        p2p_cl_inout_mpi.cpu_funcs[0] = directInoutPassCallbackMpi;
        p2p_cl_inout_mpi.nbuffers = 2;
        p2p_cl_inout_mpi.modes[0] = STARPU_RW;
        p2p_cl_inout_mpi.modes[1] = STARPU_R;

        memset(&m2l_cl_inout_mpi, 0, sizeof(m2l_cl_inout_mpi));
        m2l_cl_inout_mpi.where = STARPU_CPU;
        m2l_cl_inout_mpi.cpu_funcs[0] = transferInoutPassCallbackMpi;
        m2l_cl_inout_mpi.nbuffers = 2;
        m2l_cl_inout_mpi.modes[0] = STARPU_RW;
        m2l_cl_inout_mpi.modes[1] = STARPU_R;
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

                starpu_mpi_isend_detached( handles[tree->getHeight()][localId], sd.dest,
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

                starpu_mpi_isend_detached( handles[sd.level][localId], sd.dest,
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
            handles[idxLevel].resize(tree->getNbCellGroupAtLevel(idxLevel));

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
                const CellContainerClass* currentCells = tree->getCellGroup(idxLevel, idxGroup);
                starpu_variable_data_register(&handles[idxLevel][idxGroup], 0,
                                              (uintptr_t)currentCells->getRawBuffer(), currentCells->getBufferSizeInByte());
            }
        }
        {
            const int idxLevel = tree->getHeight();
            handles[idxLevel].resize(tree->getNbParticleGroup());
            for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
                ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);
                starpu_variable_data_register(&handles[idxLevel][idxGroup], 0,
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
                    STARPU_VALUE, &thisptr, sizeof(ThisClass*),
                    STARPU_RW, handles[tree->getHeight()-1][idxGroup],
                    STARPU_R, handles[tree->getHeight()][idxGroup],
                    0);
        }

        FLOG( FLog::Controller << "\t\t bottomPass in " << timer.tacAndElapsed() << "s\n" );
    }

    static void bottomPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass leafCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                            STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        ParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                            STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);
        worker->bottomPassPerform(&leafCells, &containers);
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

            for(int idxGroup = 0 ; idxGroup < tree->getNbCellGroupAtLevel(idxLevel)
                    && idxSubGroup < tree->getNbCellGroupAtLevel(idxLevel+1) ; ++idxGroup){
                CellContainerClass*const currentCells = tree->getCellGroup(idxLevel, idxGroup);

                struct starpu_task* const task = starpu_task_create();
                task->dyn_handles = (starpu_data_handle_t*)malloc(sizeof(starpu_data_handle_t)*10);
                task->dyn_handles[0] = handles[idxLevel][idxGroup];

                // Skip current group if needed
                if( tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++idxSubGroup;
                    FAssertLF( idxSubGroup != tree->getNbCellGroupAtLevel(idxLevel+1) );
                    FAssertLF( (tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                task->dyn_handles[nbSubCellGroups + 1] = handles[idxLevel+1][idxSubGroup];
                nbSubCellGroups += 1;
                while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= ((currentCells->getEndingIndex()<<3)+7)
                      && (idxSubGroup+1) != tree->getNbCellGroupAtLevel(idxLevel+1)
                      && tree->getCellGroup(idxLevel+1, idxSubGroup+1)->getStartingIndex() <= (currentCells->getEndingIndex()<<3)+7 ){
                    idxSubGroup += 1;
                    task->dyn_handles[nbSubCellGroups + 1] = handles[idxLevel+1][idxSubGroup];
                    nbSubCellGroups += 1;
                    FAssertLF( nbSubCellGroups <= 9 );
                }

                // put the right codelet
                task->cl = &m2m_cl[nbSubCellGroups-1];
                // put args values
                char *arg_buffer;
                size_t arg_buffer_size;
                starpu_codelet_pack_args((void**)&arg_buffer, &arg_buffer_size,
                                         STARPU_VALUE, &thisptr, sizeof(ThisClass*),
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
                    task->dyn_handles[0] = handles[idxLevel][tree->getNbCellGroupAtLevel(idxLevel)-1];

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
                                             STARPU_VALUE, &thisptr, sizeof(ThisClass*),
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

                        starpu_mpi_isend_detached( handles[idxLevel+1][lowerIdxToSend], dest,
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

    static void upwardPassCallback(void *buffers[], void *cl_arg){
        CellContainerClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));

        ThisClass* worker = nullptr;
        int nbSubCellGroups = 0;
        int idxLevel = 0;
        starpu_codelet_unpack_args(cl_arg, &worker, &nbSubCellGroups, &idxLevel);

        CellContainerClass* subCellGroups[9];
        memset(subCellGroups, 0, 9*sizeof(CellContainerClass*));
        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            subCellGroups[idxSubGroup] = new CellContainerClass((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[idxSubGroup+1]),
                    STARPU_VARIABLE_GET_ELEMSIZE(buffers[idxSubGroup+1]));
        }

        worker->upwardPassPerform(&currentCells, subCellGroups, nbSubCellGroups, idxLevel);

        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            delete subCellGroups[idxSubGroup];
        }
    }

    void upwardPassPerform(CellContainerClass*const currentCells,
                           CellContainerClass* subCellGroups[9],
                            const int nbSubCellGroups, const int idxLevel){
        FAssertLF(nbSubCellGroups != 0);
        const MortonIndex blockStartIdx = FMath::Max(currentCells->getStartingIndex(),
                                              subCellGroups[0]->getStartingIndex()>>3);
        const MortonIndex blockEndIdx   = FMath::Min(currentCells->getEndingIndex(),
                                              ((subCellGroups[nbSubCellGroups-1]->getEndingIndex()-1)>>3)+1);
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
                            STARPU_VALUE, &thisptr, sizeof(ThisClass*),
                            STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                            STARPU_VALUE, &outsideInteractions, sizeof(outsideInteractions),
                            STARPU_RW, handles[idxLevel][idxGroup],
                            STARPU_R, remoteCellGroups[idxLevel][interactionid].handle,
                            0);
                }
            }
        }
        FLOG( FLog::Controller << "\t\t transferPassMpi in " << timer.tacAndElapsed() << "s\n" );
    }

    static void transferInoutPassCallbackMpi(void *buffers[], void *cl_arg){
        CellContainerClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        CellContainerClass externalCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        ThisClass* worker = nullptr;
        int idxLevel = 0;
        const std::vector<OutOfBlockInteraction>* outsideInteractions;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel, &outsideInteractions);

        worker->transferInoutPassPerformMpi(&currentCells, &externalCells, idxLevel, outsideInteractions);
    }


    void transferInoutPassPerformMpi(CellContainerClass*const currentCells,
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
                        STARPU_VALUE, &thisptr, sizeof(ThisClass*),
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
                            STARPU_VALUE, &thisptr, sizeof(ThisClass*),
                            STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                            STARPU_VALUE, &outsideInteractions, sizeof(outsideInteractions),
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
        CellContainerClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));

        ThisClass* worker = nullptr;
        int idxLevel = 0;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel);

        worker->transferInPassPerform(&currentCells, idxLevel);
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
        CellContainerClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        CellContainerClass externalCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        ThisClass* worker = nullptr;
        int idxLevel = 0;
        const std::vector<OutOfBlockInteraction>* outsideInteractions;
        starpu_codelet_unpack_args(cl_arg, &worker, &idxLevel, &outsideInteractions);

        worker->transferInoutPassPerform(&currentCells, &externalCells, idxLevel, outsideInteractions);
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

                        starpu_mpi_isend_detached( handles[idxLevel][idxLastBlock],
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
                        task->dyn_handles[nbSubCellGroups + 1] = handles[idxLevel+1][idxSubGroup];
                        nbSubCellGroups += 1;
                        while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= ((parentEndingIdx<<3)+7)
                              && (idxSubGroup + 1) != tree->getNbCellGroupAtLevel(idxLevel+1)
                              && tree->getCellGroup(idxLevel+1, idxSubGroup+1)->getStartingIndex() <= (parentEndingIdx<<3)+7 ){
                            idxSubGroup += 1;
                            task->dyn_handles[nbSubCellGroups + 1] = handles[idxLevel+1][idxSubGroup];
                            nbSubCellGroups += 1;
                            FAssertLF( nbSubCellGroups <= 9 );
                        }

                        // put the right codelet
                        task->cl = &l2l_cl[nbSubCellGroups-1];
                        // put args values
                        char *arg_buffer;
                        size_t arg_buffer_size;
                        starpu_codelet_pack_args((void**)&arg_buffer, &arg_buffer_size,
                                                 STARPU_VALUE, &thisptr, sizeof(ThisClass*),
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
                task->dyn_handles[0] = handles[idxLevel][idxGroup];

                // Skip current group if needed
                if( tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++idxSubGroup;
                    FAssertLF( idxSubGroup != tree->getNbCellGroupAtLevel(idxLevel+1) );
                    FAssertLF( (tree->getCellGroup(idxLevel+1, idxSubGroup)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                task->dyn_handles[nbSubCellGroups + 1] = handles[idxLevel+1][idxSubGroup];
                nbSubCellGroups += 1;
                while(tree->getCellGroup(idxLevel+1, idxSubGroup)->getEndingIndex() <= ((currentCells->getEndingIndex()<<3)+7)
                      && (idxSubGroup+1) != tree->getNbCellGroupAtLevel(idxLevel+1)
                      && tree->getCellGroup(idxLevel+1, idxSubGroup+1)->getStartingIndex() <= (currentCells->getEndingIndex()<<3)+7 ){
                    idxSubGroup += 1;
                    task->dyn_handles[nbSubCellGroups + 1] = handles[idxLevel+1][idxSubGroup];
                    nbSubCellGroups += 1;
                    FAssertLF( nbSubCellGroups <= 9 );
                }

                // put the right codelet
                task->cl = &l2l_cl[nbSubCellGroups-1];
                // put args values
                char *arg_buffer;
                size_t arg_buffer_size;
                starpu_codelet_pack_args((void**)&arg_buffer, &arg_buffer_size,
                                         STARPU_VALUE, &thisptr, sizeof(ThisClass*),
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
        CellContainerClass currentCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                        STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));

        ThisClass* worker = nullptr;
        int nbSubCellGroups = 0;
        int idxLevel = 0;
        starpu_codelet_unpack_args(cl_arg, &worker, &nbSubCellGroups, &idxLevel);

        CellContainerClass* subCellGroups[9];
        memset(subCellGroups, 0, 9*sizeof(CellContainerClass*));
        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            subCellGroups[idxSubGroup] = new CellContainerClass((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[idxSubGroup+1]),
                    STARPU_VARIABLE_GET_ELEMSIZE(buffers[idxSubGroup+1]));
        }

        worker->downardPassPerform(&currentCells, subCellGroups, nbSubCellGroups, idxLevel);

        for(int idxSubGroup = 0; idxSubGroup < nbSubCellGroups ; ++idxSubGroup){
            delete subCellGroups[idxSubGroup];
        }
    }

    void downardPassPerform(CellContainerClass*const currentCells,
                            CellContainerClass* subCellGroups[9],
                             const int nbSubCellGroups, const int idxLevel){
        FAssertLF(nbSubCellGroups != 0);
        const MortonIndex blockStartIdx = FMath::Max(currentCells->getStartingIndex(),
                                              subCellGroups[0]->getStartingIndex()>>3);
        const MortonIndex blockEndIdx   = FMath::Min(currentCells->getEndingIndex(),
                                              ((subCellGroups[nbSubCellGroups-1]->getEndingIndex()-1)>>3)+1);
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
    /// Direct Pass MPI
    /////////////////////////////////////////////////////////////////////////////////////

    void directPassMpi(){
        FLOG( FTic timer; );

        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            for(int idxInteraction = 0; idxInteraction < int(externalInteractionsLeafLevelMpi[idxGroup].size()) ; ++idxInteraction){
                const int interactionid = externalInteractionsLeafLevelMpi[idxGroup][idxInteraction].otherBlockId;
                const std::vector<OutOfBlockInteraction>* outsideInteractions = &externalInteractionsLeafLevelMpi[idxGroup][idxInteraction].interactions;
                starpu_insert_task(&p2p_cl_inout_mpi,
                        STARPU_VALUE, &thisptr, sizeof(ThisClass*),
                        STARPU_VALUE, &outsideInteractions, sizeof(outsideInteractions),
                        STARPU_RW, handles[tree->getHeight()][idxGroup],
                        STARPU_R, remoteParticleGroupss[interactionid].handle,
                        0);
            }
        }

        FLOG( FLog::Controller << "\t\t directPass in MPI " << timer.tacAndElapsed() << "s\n" );
    }

    static void directInoutPassCallbackMpi(void *buffers[], void *cl_arg){
        ParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                      STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        ParticleGroupClass externalContainers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                                      STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        ThisClass* worker = nullptr;
        const std::vector<OutOfBlockInteraction>* outsideInteractions = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker, &outsideInteractions);

        worker->directInoutPassPerform(&containers, &externalContainers, outsideInteractions);
    }

    void directInoutPassPerformMpi(ParticleGroupClass* containers, ParticleGroupClass* containersOther,
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
                    STARPU_VALUE, &thisptr, sizeof(ThisClass*),
                    STARPU_RW, handles[tree->getHeight()][idxGroup],
                    0);
        }
        FLOG( timerInBlock.tac() );
        FLOG( timerOutBlock.tic() );
        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            for(int idxInteraction = 0; idxInteraction < int(externalInteractionsLeafLevel[idxGroup].size()) ; ++idxInteraction){
                const int interactionid = externalInteractionsLeafLevel[idxGroup][idxInteraction].otherBlockId;
                const std::vector<OutOfBlockInteraction>* outsideInteractions = &externalInteractionsLeafLevel[idxGroup][idxInteraction].interactions;
                starpu_insert_task(&p2p_cl_inout,
                        STARPU_VALUE, &thisptr, sizeof(ThisClass*),
                        STARPU_VALUE, &outsideInteractions, sizeof(outsideInteractions),
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
        ParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                      STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);

        worker->directInPassPerform(&containers);
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
        ParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                      STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        ParticleGroupClass externalContainers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                                      STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        ThisClass* worker = nullptr;
        const std::vector<OutOfBlockInteraction>* outsideInteractions = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker, &outsideInteractions);

        worker->directInoutPassPerform(&containers, &externalContainers, outsideInteractions);
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
                    STARPU_VALUE, &thisptr, sizeof(ThisClass*),
                    STARPU_R, handles[tree->getHeight()-1][idxGroup],
                    STARPU_RW, handles[tree->getHeight()][idxGroup],
                    0);
        }

        FLOG( FLog::Controller << "\t\t L2P in " << timer.tacAndElapsed() << "s\n" );
    }

    static void mergePassCallback(void *buffers[], void *cl_arg){
        CellContainerClass leafCells((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[0]),
                                     STARPU_VARIABLE_GET_ELEMSIZE(buffers[0]));
        ParticleGroupClass containers((unsigned char*)STARPU_VARIABLE_GET_PTR(buffers[1]),
                                     STARPU_VARIABLE_GET_ELEMSIZE(buffers[1]));

        ThisClass* worker = nullptr;
        starpu_codelet_unpack_args(cl_arg, &worker);
        worker->mergePassPerform(&leafCells, &containers);
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

#endif // FGROUPTASKSTARPUMPIALGORITHM_HPP
