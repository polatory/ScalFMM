#ifndef FGROUPSEQALGORITHM_HPP
#define FGROUPSEQALGORITHM_HPP

#include "../Utils/FGlobal.hpp"
#include "../Core/FCoreCommon.hpp"
#include "../Utils/FQuickSort.hpp"
#include "../Containers/FTreeCoordinate.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FTic.hpp"

#include <list>
#include <vector>

template <class OctreeClass, class CellContainerClass, class CellClass, class KernelClass, class ParticleGroupClass, class ParticleContainerClass>
class FGroupSeqAlgorithm {
protected:
    struct OutOfBlockInteraction{
        MortonIndex outIndex;
        MortonIndex insideIndex;
        int outPosition;

        operator long long() const{
            return static_cast<long long>(outIndex);
        }
    };

    const int MaxThreads;         //< The number of threads
    OctreeClass*const tree;       //< The Tree
    KernelClass*const kernels;    //< The kernels

public:
    FGroupSeqAlgorithm(OctreeClass*const inTree, KernelClass* inKernels) : MaxThreads(1), tree(inTree), kernels(inKernels){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(kernels, "kernels cannot be null");

        FLOG(FLog::Controller << "FGroupSeqAlgorithm (Max Thread " << MaxThreads << ")\n");
    }

    ~FGroupSeqAlgorithm(){
    }

    void execute(const unsigned operationsToProceed = FFmmNearAndFarFields){
        FLOG( FLog::Controller << "\tStart FGroupSeqAlgorithm\n" );

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if( (operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P) ) directPass();
    }

protected:
    void bottomPass(){
        FLOG( FTic timer; );
        typename std::list<ParticleGroupClass*>::iterator iterParticles = tree->leavesBegin();
        const typename std::list<ParticleGroupClass*>::iterator endParticles = tree->leavesEnd();

        typename std::list<CellContainerClass*>::iterator iterCells = tree->cellsBegin(tree->getHeight()-1);
        const typename std::list<CellContainerClass*>::iterator endCells = tree->cellsEnd(tree->getHeight()-1);

        while(iterParticles != endParticles && iterCells != endCells){
            { // Can be a task(in:iterParticles, out:iterCells)
                const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                const MortonIndex blockEndIdx = (*iterCells)->getEndingIndex();

                for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                    CellClass* cell = (*iterCells)->getCell(mindex);
                    if(cell){
                        ParticleContainerClass particles = (*iterParticles)->template getLeaf<ParticleContainerClass>(mindex);
                        FAssertLF(particles.isAttachedToSomething());
                        kernels->P2M(cell, &particles);
                    }
                }
            }

            ++iterParticles;
            ++iterCells;
        }

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
                { // Can be a task(in:iterParticles, out:iterChildCells ...)
                    const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterCells)->getEndingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && iterChildCells != endChildCells; ++mindex){
                        CellClass* cell = (*iterCells)->getCell(mindex);
                        if(cell){
                            CellClass* child[8] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

                            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                                while(iterChildCells != endChildCells && (*iterChildCells)->getEndingIndex() < ((mindex<<3)+idxChild) ){
                                    ++iterChildCells;
                                }
                                if( iterChildCells == endChildCells ){
                                    break;
                                }
                                child[idxChild] = (*iterChildCells)->getCell((mindex<<3)+idxChild);
                            }

                            kernels->M2M(cell, child, idxLevel);
                        }
                    }
                }

                ++iterCells;
            }

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
                            MortonIndex interactionsIndexes[189];
                            int interactionsPosition[189];
                            FTreeCoordinate coord(mindex, idxLevel);
                            int counter = coord.getInteractionNeighbors(idxLevel,interactionsIndexes,interactionsPosition);

                            const CellClass* interactions[343];
                            memset(interactions, 0, 343*sizeof(CellClass*));
                            int counterExistingCell = 0;

                            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                    CellClass* interCell = (*iterCells)->getCell(interactionsIndexes[idxInter]);
                                    if(interCell){
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

                            kernels->M2L( cell , interactions, counterExistingCell, idxLevel);
                        }
                    }
                }


                // Manage outofblock interaction
                FQuickSort<OutOfBlockInteraction, long long, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                typename std::list<CellContainerClass*>::iterator iterLeftCells = tree->cellsBegin(idxLevel);
                int currentOutInteraction = 0;
                while(iterLeftCells != iterCells && currentOutInteraction < int(outsideInteractions.size())){
                    const MortonIndex blockStartIdx = (*iterLeftCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterLeftCells)->getEndingIndex();

                    while(outsideInteractions[currentOutInteraction].outIndex < blockStartIdx){
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
                                CellClass* cell = (*iterCells)->getCell(outsideInteractions[outInterIdx].insideIndex);
                                FAssertLF(cell);
                                const CellClass* interactions[343];
                                memset(interactions, 0, 343*sizeof(CellClass*));
                                interactions[outsideInteractions[outInterIdx].outPosition] = interCell;
                                const int counter = 1;
                                kernels->M2L( cell , interactions, counter, idxLevel);

                                interactions[outsideInteractions[outInterIdx].outPosition] = NULL;
                                interactions[getOppositeInterIndex(outsideInteractions[outInterIdx].outPosition)] = cell;
                                kernels->M2L( interCell , interactions, counter, idxLevel);
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
                { // Can be a task(in:iterParticles, inout:iterChildCells ...)
                    const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterCells)->getEndingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && iterChildCells != endChildCells; ++mindex){
                        CellClass* cell = (*iterCells)->getCell(mindex);
                        if(cell){
                            CellClass* child[8] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

                            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                                while(iterChildCells != endChildCells && (*iterChildCells)->getEndingIndex() < ((mindex<<3)+idxChild) ){
                                    ++iterChildCells;
                                }
                                if( iterChildCells == endChildCells ){
                                    break;
                                }
                                child[idxChild] = (*iterChildCells)->getCell((mindex<<3)+idxChild);
                            }

                            kernels->L2L(cell, child, idxLevel);
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
        {
            typename std::list<ParticleGroupClass*>::iterator iterParticles = tree->leavesBegin();
            const typename std::list<ParticleGroupClass*>::iterator endParticles = tree->leavesEnd();

            typename std::list<CellContainerClass*>::iterator iterCells = tree->cellsBegin(tree->getHeight()-1);
            const typename std::list<CellContainerClass*>::iterator endCells = tree->cellsEnd(tree->getHeight()-1);

            while(iterParticles != endParticles && iterCells != endCells){
                { // Can be a task(in:iterCells, inout:iterParticles)
                    const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterCells)->getStartingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                        CellClass* cell = (*iterCells)->getCell(mindex);
                        if(cell){
                            ParticleContainerClass particles = (*iterParticles)->template getLeaf<ParticleContainerClass>(mindex);
                            FAssertLF(particles.isAttachedToSomething());
                            kernels->P2M(cell, &particles);
                        }
                    }
                }

                ++iterParticles;
                ++iterCells;
            }

            FAssertLF(iterParticles == endParticles && iterCells == endCells);
        }
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

                            kernels->P2P( coord, &particles, &particles , interactions, counterExistingCell);
                        }
                    }
                }


                // Manage outofblock interaction
                FQuickSort<OutOfBlockInteraction, long long, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                typename std::list<ParticleGroupClass*>::iterator iterLeftParticles = tree->leavesBegin();
                int currentOutInteraction = 0;
                while(iterLeftParticles != iterParticles && currentOutInteraction < int(outsideInteractions.size())){
                    const MortonIndex blockStartIdx = (*iterLeftParticles)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterLeftParticles)->getEndingIndex();

                    while(outsideInteractions[currentOutInteraction].outIndex < blockStartIdx){
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
                                kernels->P2PRemote( FTreeCoordinate(outsideInteractions[outInterIdx].insideIndex, tree->getHeight()-1), &particles, &particles , interactions, counter);

                                interactions[outsideInteractions[outInterIdx].outPosition] = NULL;
                                interactions[getOppositeNeighIndex(outsideInteractions[outInterIdx].outPosition)] = &particles;
                                kernels->P2PRemote( FTreeCoordinate(outsideInteractions[outInterIdx].outIndex, tree->getHeight()-1), &interParticles, &interParticles , interactions, counter);
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

    int getOppositeNeighIndex(const int index) const {
        // ((idxX+1)*3 + (idxY+1)) * 3 + (idxZ+1)
        return 27-index;
    }

    int getOppositeInterIndex(const int index) const {
        // ((( (xdiff+3) * 7) + (ydiff+3))) * 7 + zdiff + 3
        return 343-index;
    }
};


#endif // FGROUPSEQALGORITHM_HPP
