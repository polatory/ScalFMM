#ifndef FGROUPTASKALGORITHM_HPP
#define FGROUPTASKALGORITHM_HPP

#include "../Utils/FGlobal.hpp"
#include "../Core/FCoreCommon.hpp"
#include "../Utils/FQuickSort.hpp"
#include "../Containers/FTreeCoordinate.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FTic.hpp"

#include <list>
#include <vector>

#include <omp.h>

template <class OctreeClass, class CellContainerClass, class CellClass, class KernelClass, class ParticleGroupClass, class ParticleContainerClass>
class FGroupTaskAlgorithm {
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

    std::vector< std::list< std::list<BlockInteractions<CellContainerClass>>>> externalInteractionsAllLevel;
    std::list< std::list<BlockInteractions<ParticleGroupClass>>> externalInteractionsLeafLevel;

    const int MaxThreads;         //< The number of threads
    OctreeClass*const tree;       //< The Tree
    KernelClass** kernels;        //< The kernels

public:
    FGroupTaskAlgorithm(OctreeClass*const inTree, KernelClass* inKernels, const int inMaxThreads = -1)
        : MaxThreads(inMaxThreads==-1?omp_get_max_threads():inMaxThreads), tree(inTree), kernels(nullptr){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");

        kernels = new KernelClass*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass(*inKernels);
        }

        FLOG(FLog::Controller << "FGroupTaskAlgorithm (Max Thread " << MaxThreads << ")\n");
    }

    ~FGroupTaskAlgorithm(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete[] kernels;
    }

    void execute(const unsigned operationsToProceed = FFmmNearAndFarFields){
        FLOG( FLog::Controller << "\tStart FGroupTaskAlgorithm\n" );

        // For now rebuild all external interaction
        buildExternalInteractionVecs();

#pragma omp parallel
        {
#pragma omp single nowait
            {
                if(operationsToProceed & FFmmP2M) bottomPass();

                if(operationsToProceed & FFmmM2M) upwardPass();

                if(operationsToProceed & FFmmM2L) transferPass();

                if(operationsToProceed & FFmmL2L) downardPass();

            }

#pragma omp single nowait
            {
                if( operationsToProceed & FFmmP2P ) directPass();
            }

#pragma omp barrier

#pragma omp single nowait
            {
                if( operationsToProceed & FFmmL2P ) mergePass();
            }
        }
    }

protected:
    void buildExternalInteractionVecs(){
        // Reset interactions
        externalInteractionsAllLevel.clear();
        externalInteractionsLeafLevel.clear();
        // One per level + leaf level
        externalInteractionsAllLevel.resize(tree->getHeight());

        // First leaf level
        {
            // We create one big vector per block
            typename std::list< std::vector<OutOfBlockInteraction> > allOutsideInteractions;

            {
                typename std::list<ParticleGroupClass*>::iterator iterParticles = tree->leavesBegin();
                const typename std::list<ParticleGroupClass*>::iterator endParticles = tree->leavesEnd();

                // We iterate on blocks
                while(iterParticles != endParticles){
                    // Create the vector
                    allOutsideInteractions.push_back( std::vector<OutOfBlockInteraction>() );
                    typename std::vector<OutOfBlockInteraction>* outsideInteractions = &allOutsideInteractions.back();
                    ParticleGroupClass* containers = (*iterParticles);

                    #pragma omp task default(none) firstprivate(containers, outsideInteractions)
                    { // Can be a task(inout:iterCells, out:outsideInteractions)
                        const MortonIndex blockStartIdx = containers->getStartingIndex();
                        const MortonIndex blockEndIdx = containers->getEndingIndex();

                        for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                            ParticleContainerClass particles = containers->template getLeaf<ParticleContainerClass>(mindex);
                            if(particles.isAttachedToSomething()){
                                MortonIndex interactionsIndexes[26];
                                int interactionsPosition[26];
                                FTreeCoordinate coord(mindex, tree->getHeight()-1);
                                int counter = coord.getNeighborsIndexes(tree->getHeight(),interactionsIndexes,interactionsPosition);

                                ParticleContainerClass* interactions[27];
                                memset(interactions, 0, 27*sizeof(ParticleContainerClass*));

                                for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                                    if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                        // Inside block interaction, do nothing
                                    }
                                    else if(interactionsIndexes[idxInter] < mindex){
                                        OutOfBlockInteraction property;
                                        property.insideIndex = mindex;
                                        property.outIndex    = interactionsIndexes[idxInter];
                                        property.outPosition = interactionsPosition[idxInter];
                                        (*outsideInteractions).push_back(property);
                                    }
                                }
                            }
                        }

                        // Sort to match external order
                        FQuickSort<OutOfBlockInteraction, int>::QsSequential((*outsideInteractions).data(),int((*outsideInteractions).size()));
                    }
                    ++iterParticles;
                }
                #pragma omp taskwait
            }
            {
                // Now we have to split the vector of interaction based on the existing blocks
                typename std::list<ParticleGroupClass*>::iterator iterParticles = tree->leavesBegin();
                const typename std::list<ParticleGroupClass*>::iterator endParticles = tree->leavesEnd();

                typename std::list<std::vector<OutOfBlockInteraction> >::iterator iterInteractions = allOutsideInteractions.begin();

                // For each block (can be moved to a task or into the previous one)
                while(iterParticles != endParticles){
                    typename std::vector<OutOfBlockInteraction>* outsideInteractions = &(*iterInteractions);

                    externalInteractionsLeafLevel.emplace_back();
                    std::list<BlockInteractions<ParticleGroupClass>>* externalInteractions = &externalInteractionsLeafLevel.back();

                    typename std::list<ParticleGroupClass*>::iterator iterLeftParticles = tree->leavesBegin();
                    int currentOutInteraction = 0;
                    while(iterLeftParticles != iterParticles && currentOutInteraction < int((*outsideInteractions).size())){
                        ParticleGroupClass* leftContainers = (*iterLeftParticles);
                        const MortonIndex blockStartIdx = leftContainers->getStartingIndex();
                        const MortonIndex blockEndIdx = leftContainers->getEndingIndex();

                        while(currentOutInteraction < int((*outsideInteractions).size()) && (*outsideInteractions)[currentOutInteraction].outIndex < blockStartIdx){
                            currentOutInteraction += 1;
                        }

                        int lastOutInteraction = currentOutInteraction;
                        while(lastOutInteraction < int((*outsideInteractions).size()) && (*outsideInteractions)[lastOutInteraction].outIndex < blockEndIdx){
                            lastOutInteraction += 1;
                        }

                        const int nbInteractionsBetweenBlocks = (lastOutInteraction-currentOutInteraction);
                        if(nbInteractionsBetweenBlocks){
                            externalInteractions->emplace_back();
                            BlockInteractions<ParticleGroupClass>* interactions = &externalInteractions->back();
                            interactions->otherBlock = (*iterLeftParticles);
                            interactions->interactions.resize(nbInteractionsBetweenBlocks);
                            std::copy((*outsideInteractions).begin() + currentOutInteraction,
                                      (*outsideInteractions).begin() + lastOutInteraction,
                                      interactions->interactions.begin());
                        }

                        currentOutInteraction = lastOutInteraction;
                        ++iterLeftParticles;
                    }

                    ++iterParticles;
                    ++iterInteractions;
                }
            }
        }
    }


    void bottomPass(){
        FLOG( FTic timer; );
        typename std::list<ParticleGroupClass*>::iterator iterParticles = tree->leavesBegin();
        const typename std::list<ParticleGroupClass*>::iterator endParticles = tree->leavesEnd();

        typename std::list<CellContainerClass*>::iterator iterCells = tree->cellsBegin(tree->getHeight()-1);
        const typename std::list<CellContainerClass*>::iterator endCells = tree->cellsEnd(tree->getHeight()-1);

        while(iterParticles != endParticles && iterCells != endCells){
            CellContainerClass* leafCells = (*iterCells);
            ParticleGroupClass* containers = (*iterParticles);
#pragma omp task default(shared) firstprivate(leafCells, containers)
            { // Can be a task(in:iterParticles, out:iterCells)
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

            ++iterParticles;
            ++iterCells;
        }
        // Wait for task to complete
#pragma omp taskwait

        FAssertLF(iterParticles == endParticles && iterCells == endCells);
        FLOG( FLog::Controller << "\t\t bottomPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void upwardPass(){
        FLOG( FTic timer; );
        for(int idxLevel = tree->getHeight()-2 ; idxLevel >= 2 ; --idxLevel){
            typename std::list<CellContainerClass*>::iterator iterCells = tree->cellsBegin(idxLevel);
            const typename std::list<CellContainerClass*>::iterator endCells = tree->cellsEnd(idxLevel);

            typename std::list<CellContainerClass*>::iterator iterChildCells = tree->cellsBegin(idxLevel+1);
            const typename std::list<CellContainerClass*>::iterator endChildCells = tree->cellsEnd(idxLevel+1);

            while(iterCells != endCells && iterChildCells != endChildCells){
                CellContainerClass* currentCells = (*iterCells);

                CellContainerClass* subCellGroups[9];
                memset(subCellGroups, 0, sizeof(CellContainerClass*) * 9);

                // Skip current group if needed
                if( (*iterChildCells)->getEndingIndex() <= ((*iterCells)->getStartingIndex()<<3) ){
                    ++iterChildCells;
                    FAssertLF( iterChildCells != endChildCells );
                    FAssertLF( ((*iterChildCells)->getStartingIndex()>>3) == (*iterCells)->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                subCellGroups[nbSubCellGroups] = (*iterChildCells);
                nbSubCellGroups += 1;
                while((*iterChildCells)->getEndingIndex() <= (((*iterCells)->getEndingIndex()<<3)+7)
                      && (++iterChildCells) != endChildCells
                      && (*iterChildCells)->getStartingIndex() <= ((*iterCells)->getEndingIndex()<<3)+7 ){
                    subCellGroups[nbSubCellGroups] = (*iterChildCells);
                    nbSubCellGroups += 1;
                    FAssertLF( nbSubCellGroups <= 9 );
                }

#pragma omp task default(none) firstprivate(idxLevel, currentCells, subCellGroups, nbSubCellGroups)
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

                            kernel->M2M(cell, child, idxLevel);
                        }
                    }
                }

                ++iterCells;
            }
            // Wait this level before the next one
#pragma omp taskwait

            FAssertLF(iterCells == endCells);
            FAssertLF((iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
            FAssertLF(iterCells == endCells && (iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
        }
        FLOG( FLog::Controller << "\t\t upwardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void transferPass(){
        FLOG( FTic timer; );
        FLOG( FTic timerInBlock; FTic timerSort; FTic timerOutBlock; );
        for(int idxLevel = tree->getHeight()-1 ; idxLevel >= 2 ; --idxLevel){
            std::list<std::vector<OutOfBlockInteraction> > allOutsideInteractions;
            FLOG( timerInBlock.tic() );
            {
                typename std::list<CellContainerClass*>::iterator iterCells = tree->cellsBegin(idxLevel);
                const typename std::list<CellContainerClass*>::iterator endCells = tree->cellsEnd(idxLevel);

                while(iterCells != endCells){
                    allOutsideInteractions.push_back(std::vector<OutOfBlockInteraction>());
                    std::vector<OutOfBlockInteraction>* outsideInteractions = &allOutsideInteractions.back();

                    CellContainerClass* currentCells = (*iterCells);

                    #pragma omp task default(none) firstprivate(currentCells, outsideInteractions, idxLevel)
                    { // Can be a task(inout:iterCells, out:outsideInteractions)
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
                                    else if(interactionsIndexes[idxInter] < mindex){
                                        OutOfBlockInteraction property;
                                        property.insideIndex = mindex;
                                        property.outIndex    = interactionsIndexes[idxInter];
                                        property.outPosition = interactionsPosition[idxInter];
                                        (*outsideInteractions).push_back(property);
                                    }
                                }

                                kernel->M2L( cell , interactions, counterExistingCell, idxLevel);
                            }
                        }
                    }
                    ++iterCells;
                }
                #pragma omp taskwait
            }
            FLOG( timerInBlock.tac() );
            FLOG( timerSort.tic() );
            {
                typename std::list<CellContainerClass*>::iterator iterCells = tree->cellsBegin(idxLevel);
                const typename std::list<CellContainerClass*>::iterator endCells = tree->cellsEnd(idxLevel);

                typename std::list<std::vector<OutOfBlockInteraction> >::iterator iterInteractions = allOutsideInteractions.begin();

                while(iterCells != endCells){
                    typename std::vector<OutOfBlockInteraction>* outsideInteractions = &(*iterInteractions);
#pragma omp task default(none) firstprivate(outsideInteractions)
                    {
                        // Manage outofblock interaction
                        FQuickSort<OutOfBlockInteraction, int>::QsSequential((*outsideInteractions).data(),int((*outsideInteractions).size()));
                    }
                    ++iterCells;
                    ++iterInteractions;
                }

#pragma omp taskwait
            }
            FLOG( timerSort.tac() );
            FLOG( timerOutBlock.tic() );
            {
                typename std::list<CellContainerClass*>::iterator iterCells = tree->cellsBegin(idxLevel);
                const typename std::list<CellContainerClass*>::iterator endCells = tree->cellsEnd(idxLevel);

                typename std::list<std::vector<OutOfBlockInteraction> >::iterator iterInteractions = allOutsideInteractions.begin();

                while(iterCells != endCells){
                    CellContainerClass* currentCells = (*iterCells);

                    typename std::vector<OutOfBlockInteraction>* outsideInteractions = &(*iterInteractions);
                    typename std::list<CellContainerClass*>::iterator iterLeftCells = tree->cellsBegin(idxLevel);
                    int currentOutInteraction = 0;
                    while(iterLeftCells != iterCells && currentOutInteraction < int((*outsideInteractions).size())){
                        CellContainerClass* leftCells = (*iterLeftCells);
                        const MortonIndex blockStartIdx = (*iterLeftCells)->getStartingIndex();
                        const MortonIndex blockEndIdx = (*iterLeftCells)->getEndingIndex();

                        while(currentOutInteraction < int((*outsideInteractions).size()) && (*outsideInteractions)[currentOutInteraction].outIndex < blockStartIdx){
                            currentOutInteraction += 1;
                        }

                        int lastOutInteraction = currentOutInteraction;
                        while(lastOutInteraction < int((*outsideInteractions).size()) && (*outsideInteractions)[lastOutInteraction].outIndex < blockEndIdx){
                            lastOutInteraction += 1;
                        }

                        #pragma omp task default(none) firstprivate(currentCells, leftCells, outsideInteractions, currentOutInteraction, lastOutInteraction, idxLevel)
                        { // Can be a task(in:currentOutInteraction, in:outsideInteractions, in:lastOutInteraction, inout:iterLeftCells, inout:iterCells)
                            KernelClass*const kernel = kernels[omp_get_thread_num()];

                            for(int outInterIdx = currentOutInteraction ; outInterIdx < lastOutInteraction ; ++outInterIdx){
                                CellClass* interCell = leftCells->getCell((*outsideInteractions)[outInterIdx].outIndex);
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

                        #pragma omp taskwait

                        currentOutInteraction = lastOutInteraction;
                        ++iterLeftCells;
                    }

                    ++iterCells;
                    ++iterInteractions;
                }
            }
            FLOG( timerOutBlock.tac() );
        }
        FLOG( FLog::Controller << "\t\t transferPass in " << timer.tacAndElapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t inblock in " << timerInBlock.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t sort in " << timerSort.elapsed() << "s\n" );
        FLOG( FLog::Controller << "\t\t\t outblock in " << timerOutBlock.elapsed() << "s\n" );
    }

    void downardPass(){
        FLOG( FTic timer; );
        for(int idxLevel = 2 ; idxLevel <= tree->getHeight()-2 ; ++idxLevel){
            typename std::list<CellContainerClass*>::iterator iterCells = tree->cellsBegin(idxLevel);
            const typename std::list<CellContainerClass*>::iterator endCells = tree->cellsEnd(idxLevel);

            typename std::list<CellContainerClass*>::iterator iterChildCells = tree->cellsBegin(idxLevel+1);
            const typename std::list<CellContainerClass*>::iterator endChildCells = tree->cellsEnd(idxLevel+1);

            while(iterCells != endCells && iterChildCells != endChildCells){
                CellContainerClass* currentCells = (*iterCells);

                CellContainerClass* subCellGroups[9];
                memset(subCellGroups, 0, sizeof(CellContainerClass*) * 9);

                // Skip current group if needed
                if( (*iterChildCells)->getEndingIndex() <= ((*iterCells)->getStartingIndex()<<3) ){
                    ++iterChildCells;
                    FAssertLF( iterChildCells != endChildCells );
                    FAssertLF( ((*iterChildCells)->getStartingIndex()>>3) == (*iterCells)->getStartingIndex() );
                }
                // Copy at max 8 groups
                int nbSubCellGroups = 0;
                subCellGroups[nbSubCellGroups] = (*iterChildCells);
                nbSubCellGroups += 1;
                while((*iterChildCells)->getEndingIndex() <= (((*iterCells)->getEndingIndex()<<3)+7)
                      && (++iterChildCells) != endChildCells
                      && (*iterChildCells)->getStartingIndex() <= ((*iterCells)->getEndingIndex()<<3)+7 ){
                    subCellGroups[nbSubCellGroups] = (*iterChildCells);
                    nbSubCellGroups += 1;
                    FAssertLF( nbSubCellGroups <= 9 );
                }

#pragma omp task default(none) firstprivate(idxLevel, currentCells, subCellGroups, nbSubCellGroups)
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

#pragma omp taskwait

            FAssertLF(iterCells == endCells && (iterChildCells == endChildCells || (++iterChildCells) == endChildCells));
        }
        FLOG( FLog::Controller << "\t\t downardPass in " << timer.tacAndElapsed() << "s\n" );
    }

    void directPass(){
        FLOG( FTic timer; );
        FLOG( FTic timerInBlock; FTic timerOutBlock; );

        FLOG( timerInBlock.tic() );
        {
            typename std::list<ParticleGroupClass*>::iterator iterParticles = tree->leavesBegin();
            const typename std::list<ParticleGroupClass*>::iterator endParticles = tree->leavesEnd();

            while(iterParticles != endParticles){
                ParticleGroupClass* containers = (*iterParticles);

                #pragma omp task default(none) firstprivate(containers)
                { // Can be a task(inout:iterCells, out:outsideInteractions)
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
            #pragma omp taskwait
        }
        FLOG( timerInBlock.tac() );
        FLOG( timerOutBlock.tic() );
        {
            typename std::list<ParticleGroupClass*>::iterator iterParticles = tree->leavesBegin();
            const typename std::list<ParticleGroupClass*>::iterator endParticles = tree->leavesEnd();

            typename std::list<std::list<BlockInteractions<ParticleGroupClass>>>::iterator externalInteractionsIter = externalInteractionsLeafLevel.begin();

            while(iterParticles != endParticles){
                typename std::list<BlockInteractions<ParticleGroupClass>>::iterator currentInteractions = (*externalInteractionsIter).begin();
                const typename std::list<BlockInteractions<ParticleGroupClass>>::iterator currentInteractionsEnd = (*externalInteractionsIter).end();

                ParticleGroupClass* containers = (*iterParticles);

                while(currentInteractions != currentInteractionsEnd){
                    ParticleGroupClass* containersOther = (*currentInteractions).otherBlock;
                    const std::vector<OutOfBlockInteraction>* outsideInteractions = &(*currentInteractions).interactions;

                    #pragma omp task default(none) firstprivate(containers, containersOther, outsideInteractions)
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
                    // only one task but need to wait for it
                    #pragma omp taskwait

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
        {
            typename std::list<ParticleGroupClass*>::iterator iterParticles = tree->leavesBegin();
            const typename std::list<ParticleGroupClass*>::iterator endParticles = tree->leavesEnd();

            typename std::list<CellContainerClass*>::iterator iterCells = tree->cellsBegin(tree->getHeight()-1);
            const typename std::list<CellContainerClass*>::iterator endCells = tree->cellsEnd(tree->getHeight()-1);

            while(iterParticles != endParticles && iterCells != endCells){
                CellContainerClass* leafCells = (*iterCells);
                ParticleGroupClass* containers = (*iterParticles);
#pragma omp task default(shared) firstprivate(leafCells, containers)
                { // Can be a task(out:iterParticles, in:iterCells)
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

                ++iterParticles;
                ++iterCells;
            }
            // Wait for task to complete
#pragma omp taskwait

            FAssertLF(iterParticles == endParticles && iterCells == endCells);
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

#endif // FGROUPTASKALGORITHM_HPP
