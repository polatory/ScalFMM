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

    const int MaxThreads;         //< The number of threads
    OctreeClass*const tree;       //< The Tree
    KernelClass** kernels;    //< The kernels

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
        for(int idxLevel = tree->getHeight()-1 ; idxLevel >= 2 ; --idxLevel){
            typename std::list<CellContainerClass*>::iterator iterCells = tree->cellsBegin(idxLevel);
            const typename std::list<CellContainerClass*>::iterator endCells = tree->cellsEnd(idxLevel);

            while(iterCells != endCells){
                std::vector<OutOfBlockInteraction> outsideInteractions;

                { // Can be a task(inout:iterCells, out:outsideInteractions)
                    const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterCells)->getEndingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                        CellClass* cell = (*iterCells)->getCell(mindex);
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
                                    CellClass* interCell = (*iterCells)->getCell(interactionsIndexes[idxInter]);
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
                                    outsideInteractions.push_back(property);
                                }
                            }

                            kernels[0]->M2L( cell , interactions, counterExistingCell, idxLevel);
                        }
                    }
                }


                // Manage outofblock interaction
                FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                typename std::list<CellContainerClass*>::iterator iterLeftCells = tree->cellsBegin(idxLevel);
                int currentOutInteraction = 0;
                while(iterLeftCells != iterCells && currentOutInteraction < int(outsideInteractions.size())){
                    const MortonIndex blockStartIdx = (*iterLeftCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterLeftCells)->getEndingIndex();

                    while(currentOutInteraction < int(outsideInteractions.size()) && outsideInteractions[currentOutInteraction].outIndex < blockStartIdx){
                        currentOutInteraction += 1;
                    }

                    int lastOutInteraction = currentOutInteraction;
                    while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdx){
                        lastOutInteraction += 1;
                    }

                    { // Can be a task(in:currentOutInteraction, in:outsideInteractions, in:lastOutInteraction, inout:iterLeftCells, inout:iterCells)
                        for(int outInterIdx = currentOutInteraction ; outInterIdx < lastOutInteraction ; ++outInterIdx){
                            CellClass* interCell = (*iterLeftCells)->getCell(outsideInteractions[outInterIdx].outIndex);
                            if(interCell){
                                FAssertLF(interCell->getMortonIndex() == outsideInteractions[outInterIdx].outIndex);
                                CellClass* cell = (*iterCells)->getCell(outsideInteractions[outInterIdx].insideIndex);
                                FAssertLF(cell);
                                FAssertLF(cell->getMortonIndex() == outsideInteractions[outInterIdx].insideIndex);

                                const CellClass* interactions[343];
                                memset(interactions, 0, 343*sizeof(CellClass*));
                                interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
                                const int counter = 1;
                                kernels[0]->M2L( cell , interactions, counter, idxLevel);

                                interactions[outsideInteractions[outInterIdx].outPosition] = nullptr;
                                interactions[getOppositeInterIndex(outsideInteractions[outInterIdx].outPosition)] = cell;
                                kernels[0]->M2L( interCell , interactions, counter, idxLevel);
                            }
                        }
                    }

                    currentOutInteraction = lastOutInteraction;
                    ++iterLeftCells;
                }

                ++iterCells;
            }

        }
        FLOG( FLog::Controller << "\t\t transferPass in " << timer.tacAndElapsed() << "s\n" );
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
        {
            typename std::list<ParticleGroupClass*>::iterator iterParticles = tree->leavesBegin();
            const typename std::list<ParticleGroupClass*>::iterator endParticles = tree->leavesEnd();

            while(iterParticles != endParticles){
                typename std::vector<OutOfBlockInteraction> outsideInteractions;

                { // Can be a task(inout:iterCells, out:outsideInteractions)
                    const MortonIndex blockStartIdx = (*iterParticles)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterParticles)->getEndingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                        ParticleContainerClass particles = (*iterParticles)->template getLeaf<ParticleContainerClass>(mindex);
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
                                    interactionsObjects[counterExistingCell] = (*iterParticles)->template getLeaf<ParticleContainerClass>(interactionsIndexes[idxInter]);
                                    if(interactionsObjects[counterExistingCell].isAttachedToSomething()){
                                        FAssertLF(interactions[interactionsPosition[idxInter]] == nullptr);
                                        interactions[interactionsPosition[idxInter]] = &interactionsObjects[counterExistingCell];
                                        counterExistingCell += 1;
                                    }
                                }
                                else if(interactionsIndexes[idxInter] < mindex){
                                    OutOfBlockInteraction property;
                                    property.insideIndex = mindex;
                                    property.outIndex    = interactionsIndexes[idxInter];
                                    property.outPosition = interactionsPosition[idxInter];
                                    outsideInteractions.push_back(property);
                                }
                            }

                            kernels[0]->P2P( coord, &particles, &particles , interactions, counterExistingCell);
                        }
                    }
                }


                // Manage outofblock interaction
                FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                typename std::list<ParticleGroupClass*>::iterator iterLeftParticles = tree->leavesBegin();
                int currentOutInteraction = 0;
                while(iterLeftParticles != iterParticles && currentOutInteraction < int(outsideInteractions.size())){
                    const MortonIndex blockStartIdx = (*iterLeftParticles)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterLeftParticles)->getEndingIndex();

                    while(currentOutInteraction < int(outsideInteractions.size()) && outsideInteractions[currentOutInteraction].outIndex < blockStartIdx){
                        currentOutInteraction += 1;
                    }

                    int lastOutInteraction = currentOutInteraction;
                    while(lastOutInteraction < int(outsideInteractions.size()) && outsideInteractions[lastOutInteraction].outIndex < blockEndIdx){
                        lastOutInteraction += 1;
                    }

                    { // Can be a task(in:currentOutInteraction, in:outsideInteractions, in:lastOutInteraction, inout:iterLeftParticles, inout:iterParticles)
                        for(int outInterIdx = currentOutInteraction ; outInterIdx < lastOutInteraction ; ++outInterIdx){
                            ParticleContainerClass interParticles = (*iterLeftParticles)->template getLeaf<ParticleContainerClass>(outsideInteractions[outInterIdx].outIndex);
                            if(interParticles.isAttachedToSomething()){
                                ParticleContainerClass particles = (*iterParticles)->template getLeaf<ParticleContainerClass>(outsideInteractions[outInterIdx].insideIndex);
                                FAssertLF(particles.isAttachedToSomething());
                                ParticleContainerClass* interactions[27];
                                memset(interactions, 0, 27*sizeof(ParticleContainerClass*));
                                interactions[outsideInteractions[outInterIdx].outPosition] = &interParticles;
                                const int counter = 1;
                                kernels[0]->P2PRemote( FTreeCoordinate(outsideInteractions[outInterIdx].insideIndex, tree->getHeight()-1), &particles, &particles , interactions, counter);

                                interactions[outsideInteractions[outInterIdx].outPosition] = nullptr;
                                interactions[getOppositeNeighIndex(outsideInteractions[outInterIdx].outPosition)] = &particles;
                                kernels[0]->P2PRemote( FTreeCoordinate(outsideInteractions[outInterIdx].outIndex, tree->getHeight()-1), &interParticles, &interParticles , interactions, counter);
                            }
                        }
                    }

                    currentOutInteraction = lastOutInteraction;
                    ++iterLeftParticles;
                }

                ++iterParticles;
            }
        }
        FLOG( FLog::Controller << "\t\t directPass in " << timer.tacAndElapsed() << "s\n" );
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
