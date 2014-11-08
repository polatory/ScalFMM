#ifndef FGROUPTASKDEPALGORITHM_HPP
#define FGROUPTASKDEPALGORITHM_HPP


#include "../Utils/FGlobal.hpp"
#include "../Core/FCoreCommon.hpp"
#include "../Utils/FQuickSort.hpp"
#include "../Containers/FTreeCoordinate.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FTic.hpp"

#include <vector>
#include <vector>

#include <omp.h>

template <class OctreeClass, class CellContainerClass, class CellClass, class KernelClass, class ParticleGroupClass, class ParticleContainerClass>
class FGroupTaskDepAlgorithm {
protected:
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
        std::vector<OutOfBlockInteraction> interactions;
    };

    std::vector< std::vector< std::vector<BlockInteractions<CellContainerClass>>>> externalInteractionsAllLevel;
    std::vector< std::vector<BlockInteractions<ParticleGroupClass>>> externalInteractionsLeafLevel;

    const int MaxThreads;         //< The number of threads
    OctreeClass*const tree;       //< The Tree
    KernelClass** kernels;        //< The kernels

public:
    FGroupTaskDepAlgorithm(OctreeClass*const inTree, KernelClass* inKernels, const int inMaxThreads = -1)
        : MaxThreads(inMaxThreads==-1?omp_get_max_threads():inMaxThreads), tree(inTree), kernels(nullptr){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");

        kernels = new KernelClass*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass(*inKernels);
        }

        FLOG(FLog::Controller << "FGroupTaskDepAlgorithm (Max Thread " << MaxThreads << ")\n");
    }

    ~FGroupTaskDepAlgorithm(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete[] kernels;
    }

    void execute(const unsigned operationsToProceed = FFmmNearAndFarFields){
        FLOG( FLog::Controller << "\tStart FGroupTaskDepAlgorithm\n" );

        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                // For now rebuild all external interaction
                buildExternalInteractionVecs();
            }
        }

        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                if( operationsToProceed & FFmmP2P ) directPass();

                if(operationsToProceed & FFmmP2M) bottomPass();

                if(operationsToProceed & FFmmM2M) upwardPass();

                if(operationsToProceed & FFmmM2L) transferPass();

                if(operationsToProceed & FFmmL2L) downardPass();

                if( operationsToProceed & FFmmL2P ) mergePass();

                #pragma omp taskwait
            }
        }
    }

protected:
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


    void bottomPass(){
        FLOG( FTic timer; );

        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            CellContainerClass* leafCells  = tree->getCellGroup(tree->getHeight()-1, idxGroup);
            ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);
            #pragma omp task default(shared) firstprivate(leafCells, containers) depend(inout: leafCells[0]) depend(in: containers[0])
            {
                const MortonIndex blockStartIdx = leafCells->getStartingIndex();
                const MortonIndex blockEndIdx = leafCells->getEndingIndex();
                KernelClass*const kernel = kernels[omp_get_thread_num()];

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
        }

        FLOG( FLog::Controller << "\t\t bottomPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void upwardPass(){
        FLOG( FTic timer; );
        for(int idxLevel = tree->getHeight()-2 ; idxLevel >= 2 ; --idxLevel){
            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

            typename OctreeClass::CellGroupIterator iterChildCells = tree->cellsBegin(idxLevel+1);
            const typename OctreeClass::CellGroupIterator endChildCells = tree->cellsEnd(idxLevel+1);

            while(iterCells != endCells && iterChildCells != endChildCells){
                CellContainerClass*const currentCells = (*iterCells);

                CellContainerClass* subCellGroups[9];
                memset(subCellGroups, 0, sizeof(CellContainerClass*) * 9);

                // Skip current group if needed
                if( (*iterChildCells)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++iterChildCells;
                    FAssertLF( iterChildCells != endChildCells );
                    FAssertLF( ((*iterChildCells)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                subCellGroups[nbSubCellGroups] = (*iterChildCells);
                nbSubCellGroups += 1;
                while((*iterChildCells)->getEndingIndex() <= ((currentCells->getEndingIndex()<<3)+7)
                      && (++iterChildCells) != endChildCells
                      && (*iterChildCells)->getStartingIndex() <= (currentCells->getEndingIndex()<<3)+7 ){
                    subCellGroups[nbSubCellGroups] = (*iterChildCells);
                    nbSubCellGroups += 1;
                    FAssertLF( nbSubCellGroups <= 9 );
                }

                #pragma omp task default(none) firstprivate(idxLevel, currentCells, subCellGroups, nbSubCellGroups) depend(inout: currentCells[0]) depend(in: subCellGroups[0][0], subCellGroups[1][0], subCellGroups[2][0], subCellGroups[3][0], subCellGroups[4][0], subCellGroups[5][0], subCellGroups[6][0], subCellGroups[7][0], subCellGroups[8][0])
                {
                    const MortonIndex blockStartIdx = currentCells->getStartingIndex();
                    const MortonIndex blockEndIdx   = currentCells->getEndingIndex();
                    KernelClass*const kernel = kernels[omp_get_thread_num()];
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

                ++iterCells;
            }

            FAssertLF(iterCells == endCells);
            FAssertLF((iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
            FAssertLF(iterCells == endCells && (iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
        }
        FLOG( FLog::Controller << "\t\t upwardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void transferPass(){
        FLOG( FTic timer; );
        FLOG( FTic timerInBlock; FTic timerOutBlock; );
        for(int idxLevel = tree->getHeight()-1 ; idxLevel >= 2 ; --idxLevel){
            FLOG( timerInBlock.tic() );
            {
                typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
                const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

                while(iterCells != endCells){
                    CellContainerClass* currentCells = (*iterCells);

                    #pragma omp task default(none) firstprivate(currentCells, idxLevel) depend(inout: currentCells[0])
                    {
                        const MortonIndex blockStartIdx = currentCells->getStartingIndex();
                        const MortonIndex blockEndIdx = currentCells->getEndingIndex();
                        KernelClass*const kernel = kernels[omp_get_thread_num()];

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
                    ++iterCells;
                }
            }
            FLOG( timerInBlock.tac() );
            FLOG( timerOutBlock.tic() );
            {
                typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
                const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

                typename std::vector<std::vector<BlockInteractions<CellContainerClass>>>::iterator externalInteractionsIter = externalInteractionsAllLevel[idxLevel].begin();

                while(iterCells != endCells){
                    CellContainerClass* currentCells = (*iterCells);

                    typename std::vector<BlockInteractions<CellContainerClass>>::iterator currentInteractions = (*externalInteractionsIter).begin();
                    const typename std::vector<BlockInteractions<CellContainerClass>>::iterator currentInteractionsEnd = (*externalInteractionsIter).end();

                    while(currentInteractions != currentInteractionsEnd){
                        CellContainerClass* cellsOther = (*currentInteractions).otherBlock;
                        const std::vector<OutOfBlockInteraction>* outsideInteractions = &(*currentInteractions).interactions;

                        #pragma omp task default(none) firstprivate(currentCells, outsideInteractions, cellsOther, idxLevel) depend(inout: currentCells[0], cellsOther[0])
                        {
                            KernelClass*const kernel = kernels[omp_get_thread_num()];

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

                        ++currentInteractions;
                    }

                    ++iterCells;
                    ++externalInteractionsIter;
                }
            }
            FLOG( timerOutBlock.tac() );
        }
        FLOG( FLog::Controller << "\t\t transferPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock in  " << timerInBlock.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.elapsed() << "s\n" );
    }

    void downardPass(){
        FLOG( FTic timer; );
        for(int idxLevel = 2 ; idxLevel <= tree->getHeight()-2 ; ++idxLevel){
            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

            typename OctreeClass::CellGroupIterator iterChildCells = tree->cellsBegin(idxLevel+1);
            const typename OctreeClass::CellGroupIterator endChildCells = tree->cellsEnd(idxLevel+1);

            while(iterCells != endCells && iterChildCells != endChildCells){
                CellContainerClass*const currentCells = (*iterCells);

                CellContainerClass* subCellGroups[9];
                memset(subCellGroups, 0, sizeof(CellContainerClass*) * 9);

                // Skip current group if needed
                if( (*iterChildCells)->getEndingIndex() <= (currentCells->getStartingIndex()<<3) ){
                    ++iterChildCells;
                    FAssertLF( iterChildCells != endChildCells );
                    FAssertLF( ((*iterChildCells)->getStartingIndex()>>3) == currentCells->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                subCellGroups[nbSubCellGroups] = (*iterChildCells);
                nbSubCellGroups += 1;
                while((*iterChildCells)->getEndingIndex() <= ((currentCells->getEndingIndex()<<3)+7)
                      && (++iterChildCells) != endChildCells
                      && (*iterChildCells)->getStartingIndex() <= (currentCells->getEndingIndex()<<3)+7 ){
                    subCellGroups[nbSubCellGroups] = (*iterChildCells);
                    nbSubCellGroups += 1;
                    FAssertLF( nbSubCellGroups <= 9 );
                }

                #pragma omp task default(none) firstprivate(idxLevel, currentCells, subCellGroups, nbSubCellGroups) depend(inout: subCellGroups[0][0], subCellGroups[1][0], subCellGroups[2][0], subCellGroups[3][0], subCellGroups[4][0], subCellGroups[5][0], subCellGroups[6][0], subCellGroups[7][0], subCellGroups[8][0]) depend(in: currentCells[0])
                {
                    const MortonIndex blockStartIdx = currentCells->getStartingIndex();
                    const MortonIndex blockEndIdx = currentCells->getEndingIndex();
                    KernelClass*const kernel = kernels[omp_get_thread_num()];
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

                ++iterCells;
            }

            FAssertLF(iterCells == endCells && (iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
        }
        FLOG( FLog::Controller << "\t\t downardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void directPass(){
        FLOG( FTic timer; );
        FLOG( FTic timerInBlock; FTic timerOutBlock; );

        FLOG( timerInBlock.tic() );
        {
            typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
            const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

            while(iterParticles != endParticles){
                ParticleGroupClass* containers = (*iterParticles);

                #pragma omp task default(none) firstprivate(containers) depend(inout: containers[0])
                {
                    const MortonIndex blockStartIdx = containers->getStartingIndex();
                    const MortonIndex blockEndIdx = containers->getEndingIndex();
                    KernelClass*const kernel = kernels[omp_get_thread_num()];

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
                ++iterParticles;
            }
        }
        FLOG( timerInBlock.tac() );
        FLOG( timerOutBlock.tic() );
        {
            typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
            const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

            typename std::vector<std::vector<BlockInteractions<ParticleGroupClass>>>::iterator externalInteractionsIter = externalInteractionsLeafLevel.begin();

            while(iterParticles != endParticles){
                typename std::vector<BlockInteractions<ParticleGroupClass>>::iterator currentInteractions = (*externalInteractionsIter).begin();
                const typename std::vector<BlockInteractions<ParticleGroupClass>>::iterator currentInteractionsEnd = (*externalInteractionsIter).end();

                ParticleGroupClass* containers = (*iterParticles);

                while(currentInteractions != currentInteractionsEnd){
                    ParticleGroupClass* containersOther = (*currentInteractions).otherBlock;
                    const std::vector<OutOfBlockInteraction>* outsideInteractions = &(*currentInteractions).interactions;

                    #pragma omp task default(none) firstprivate(containers, containersOther, outsideInteractions) depend(inout: containersOther[0], containers[0])
                    {
                        KernelClass*const kernel = kernels[omp_get_thread_num()];
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

                    ++currentInteractions;
                }

                ++iterParticles;
                ++externalInteractionsIter;
            }
        }
        FLOG( timerOutBlock.tac() );

        FLOG( FLog::Controller << "\t\t directPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock  in " << timerInBlock.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.elapsed() << "s\n" );
    }

    void mergePass(){
        FLOG( FTic timer; );

        for(int idxGroup = 0 ; idxGroup < tree->getNbParticleGroup() ; ++idxGroup){
            CellContainerClass* leafCells  = tree->getCellGroup(tree->getHeight()-1, idxGroup);
            ParticleGroupClass* containers = tree->getParticleGroup(idxGroup);
            #pragma omp task default(shared) firstprivate(leafCells, containers) depend(inout: containers[0]) depend(in: leafCells[0])
            {
                const MortonIndex blockStartIdx = leafCells->getStartingIndex();
                const MortonIndex blockEndIdx = leafCells->getEndingIndex();
                KernelClass*const kernel = kernels[omp_get_thread_num()];

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
        }

        FLOG( FLog::Controller << "\t\t L2P in " << timer.tacAndElapsed() << "s\n" );
    }

    int getOppositeNeighIndex(const int index) const {
        // ((idxX+1)*3 + (idxY+1)) * 3 + (idxZ+1)
        return 27-index-1;
    }

    int getOppositeInterIndex(const int index) const {
        // ((( (xdiff+3) * 7) + (ydiff+3))) * 7 + zdiff + 3
        return 343-index-1;
    }
};

#endif // FGROUPTASKDEPALGORITHM_HPP
