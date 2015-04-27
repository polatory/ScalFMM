
// Keep in private GIT
// @SCALFMM_PRIVATE
#ifndef FGROUPSEQALGORITHM_HPP
#define FGROUPSEQALGORITHM_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Core/FCoreCommon.hpp"
#include "../../Utils/FQuickSort.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../../Utils/FLog.hpp"
#include "../../Utils/FTic.hpp"

#include "FOutOfBlockInteraction.hpp"

#include <vector>
#include <vector>

template <class OctreeClass, class CellContainerClass, class CellClass, class KernelClass, class ParticleGroupClass, class ParticleContainerClass>
class FGroupSeqAlgorithm : public FAbstractAlgorithm {
protected:
    const int MaxThreads;         //< The number of threads
    OctreeClass*const tree;       //< The Tree
    KernelClass*const kernels;    //< The kernels

public:
    FGroupSeqAlgorithm(OctreeClass*const inTree, KernelClass* inKernels) : MaxThreads(1), tree(inTree), kernels(inKernels){
        FAssertLF(tree, "tree cannot be null");
        FAssertLF(kernels, "kernels cannot be null");

        FAbstractAlgorithm::setNbLevelsInTree(tree->getHeight());

        FLOG(FLog::Controller << "FGroupSeqAlgorithm (Max Thread " << MaxThreads << ")\n");
    }

    ~FGroupSeqAlgorithm(){
    }

protected:
    /**
      * Runs the complete algorithm.
      */
    void executeCore(const unsigned operationsToProceed) override {
        FLOG( FLog::Controller << "\tStart FGroupSeqAlgorithm\n" );

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if( (operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P) ) directPass();
    }

    void bottomPass(){
        FLOG( FTic timer; );
        typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
        const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

        typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(tree->getHeight()-1);
        const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(tree->getHeight()-1);

        while(iterParticles != endParticles && iterCells != endCells){
            { // Can be a task(in:iterParticles, out:iterCells)
                const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                const MortonIndex blockEndIdx = (*iterCells)->getEndingIndex();

                for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                    if((*iterCells)->exists(mindex)){
                        CellClass cell = (*iterCells)->getUpCell(mindex);
                        FAssertLF(cell.getMortonIndex() == mindex);
                        ParticleContainerClass particles = (*iterParticles)->template getLeaf<ParticleContainerClass>(mindex);
                        FAssertLF(particles.isAttachedToSomething());
                        kernels->P2M(&cell, &particles);
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
        for(int idxLevel = FMath::Min(tree->getHeight() - 2, FAbstractAlgorithm::lowerWorkingLevel - 1) ; idxLevel >= FAbstractAlgorithm::upperWorkingLevel ; --idxLevel){
            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

            typename OctreeClass::CellGroupIterator iterChildCells = tree->cellsBegin(idxLevel+1);
            const typename OctreeClass::CellGroupIterator endChildCells = tree->cellsEnd(idxLevel+1);

            while(iterCells != endCells && iterChildCells != endChildCells){
                { // Can be a task(in:iterParticles, out:iterChildCells ...)
                    const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterCells)->getEndingIndex();
                    CellClass childData[8];

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && iterChildCells != endChildCells; ++mindex){
                        if((*iterCells)->exists(mindex)){
                            CellClass cell = (*iterCells)->getUpCell(mindex);
                            FAssertLF(cell.getMortonIndex() == mindex);
                            CellClass* child[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

                            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                                while(iterChildCells != endChildCells && (*iterChildCells)->getEndingIndex() <= ((mindex<<3)+idxChild) ){
                                    ++iterChildCells;
                                }
                                if( iterChildCells == endChildCells ){
                                    break;
                                }
                                if((*iterChildCells)->exists((mindex<<3)+idxChild)){
                                    childData[idxChild] = (*iterChildCells)->getUpCell((mindex<<3)+idxChild);
                                    child[idxChild] = &childData[idxChild];
                                }
                                FAssertLF(child[idxChild] == nullptr || child[idxChild]->getMortonIndex() == ((mindex<<3)+idxChild));
                            }

                            kernels->M2M(&cell, child, idxLevel);
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
        for(int idxLevel = FAbstractAlgorithm::lowerWorkingLevel-1 ; idxLevel >= FAbstractAlgorithm::upperWorkingLevel ; --idxLevel){
            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

            while(iterCells != endCells){
                std::vector<OutOfBlockInteraction> outsideInteractions;

                { // Can be a task(inout:iterCells, out:outsideInteractions)
                    CellClass interactionsData[343];
                    const CellClass* interactions[343];
                    const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterCells)->getEndingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                        if((*iterCells)->exists(mindex)){
                            CellClass cell = (*iterCells)->getDownCell(mindex);
                            FAssertLF(cell.getMortonIndex() == mindex);
                            MortonIndex interactionsIndexes[189];
                            int interactionsPosition[189];
                            const FTreeCoordinate coord(cell.getCoordinate());
                            int counter = coord.getInteractionNeighbors(idxLevel,interactionsIndexes,interactionsPosition);

                            memset(interactions, 0, 343*sizeof(CellClass*));
                            int counterExistingCell = 0;

                            for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                                if( blockStartIdx <= interactionsIndexes[idxInter] && interactionsIndexes[idxInter] < blockEndIdx ){
                                    if((*iterCells)->exists(interactionsIndexes[idxInter])){
                                        CellClass interCell = (*iterCells)->getUpCell(interactionsIndexes[idxInter]);
                                        FAssertLF(interCell.getMortonIndex() == interactionsIndexes[idxInter]);
                                        FAssertLF(interactions[interactionsPosition[idxInter]] == nullptr);
                                        interactionsData[interactionsPosition[idxInter]] = interCell;
                                        interactions[interactionsPosition[idxInter]] = &interactionsData[interactionsPosition[idxInter]];
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

                            kernels->M2L( &cell , interactions, counterExistingCell, idxLevel);
                        }
                    }
                }


                // Manage outofblock interaction
                FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                typename OctreeClass::CellGroupIterator iterLeftCells = tree->cellsBegin(idxLevel);
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
                        const CellClass* interactions[343];
                        memset(interactions, 0, 343*sizeof(CellClass*));

                        for(int outInterIdx = currentOutInteraction ; outInterIdx < lastOutInteraction ; ++outInterIdx){
                            if((*iterLeftCells)->exists(outsideInteractions[outInterIdx].outIndex)){
                                CellClass interCell = (*iterLeftCells)->getCompleteCell(outsideInteractions[outInterIdx].outIndex);
                                FAssertLF(interCell.getMortonIndex() == outsideInteractions[outInterIdx].outIndex);
                                CellClass cell = (*iterCells)->getCompleteCell(outsideInteractions[outInterIdx].insideIndex);
                                FAssertLF(cell.getMortonIndex() == outsideInteractions[outInterIdx].insideIndex);

                                interactions[outsideInteractions[outInterIdx].outPosition] = &interCell;
                                const int counter = 1;
                                kernels->M2L( &cell , interactions, counter, idxLevel);

                                interactions[outsideInteractions[outInterIdx].outPosition] = nullptr;
                                interactions[getOppositeInterIndex(outsideInteractions[outInterIdx].outPosition)] = &cell;
                                kernels->M2L( &interCell , interactions, counter, idxLevel);
                                interactions[getOppositeInterIndex(outsideInteractions[outInterIdx].outPosition)] = nullptr;
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
        for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < FAbstractAlgorithm::lowerWorkingLevel - 1 ; ++idxLevel){
            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(idxLevel);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(idxLevel);

            typename OctreeClass::CellGroupIterator iterChildCells = tree->cellsBegin(idxLevel+1);
            const typename OctreeClass::CellGroupIterator endChildCells = tree->cellsEnd(idxLevel+1);

            while(iterCells != endCells && iterChildCells != endChildCells){
                { // Can be a task(in:iterParticles, inout:iterChildCells ...)
                    const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterCells)->getEndingIndex();
                    CellClass childData[8];

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx && iterChildCells != endChildCells; ++mindex){
                        if((*iterCells)->exists(mindex)){
                            CellClass cell = (*iterCells)->getDownCell(mindex);
                            FAssertLF(cell.getMortonIndex() == mindex);
                            CellClass* child[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

                            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                                while(iterChildCells != endChildCells && (*iterChildCells)->getEndingIndex() <= ((mindex<<3)+idxChild) ){
                                    ++iterChildCells;
                                }
                                if( iterChildCells == endChildCells ){
                                    break;
                                }
                                if((*iterChildCells)->exists((mindex<<3)+idxChild)){
                                    childData[idxChild] = (*iterChildCells)->getCompleteCell((mindex<<3)+idxChild);
                                    child[idxChild] = &childData[idxChild];
                                }
                                FAssertLF(child[idxChild] == nullptr || child[idxChild]->getMortonIndex() == ((mindex<<3)+idxChild));
                            }

                            kernels->L2L(&cell, child, idxLevel);
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
            typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
            const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

            typename OctreeClass::CellGroupIterator iterCells = tree->cellsBegin(tree->getHeight()-1);
            const typename OctreeClass::CellGroupIterator endCells = tree->cellsEnd(tree->getHeight()-1);

            while(iterParticles != endParticles && iterCells != endCells){
                { // Can be a task(in:iterCells, inout:iterParticles)
                    const MortonIndex blockStartIdx = (*iterCells)->getStartingIndex();
                    const MortonIndex blockEndIdx = (*iterCells)->getEndingIndex();

                    for(MortonIndex mindex = blockStartIdx ; mindex < blockEndIdx ; ++mindex){
                        if((*iterCells)->exists(mindex)){
                            CellClass cell = (*iterCells)->getDownCell(mindex);
                            ParticleContainerClass particles = (*iterParticles)->template getLeaf<ParticleContainerClass>(mindex);
                            FAssertLF(particles.isAttachedToSomething());
                            kernels->L2P(&cell, &particles);
                        }
                    }
                }

                ++iterParticles;
                ++iterCells;
            }

            FAssertLF(iterParticles == endParticles && iterCells == endCells);
        }
        {
            typename OctreeClass::ParticleGroupIterator iterParticles = tree->leavesBegin();
            const typename OctreeClass::ParticleGroupIterator endParticles = tree->leavesEnd();

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

                            kernels->P2P( coord, &particles, &particles , interactions, counterExistingCell);
                        }
                    }
                }


                // Manage outofblock interaction
                FQuickSort<OutOfBlockInteraction, int>::QsSequential(outsideInteractions.data(),int(outsideInteractions.size()));

                typename OctreeClass::ParticleGroupIterator iterLeftParticles = tree->leavesBegin();
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
                                kernels->P2PRemote( FTreeCoordinate(outsideInteractions[outInterIdx].insideIndex, tree->getHeight()-1), &particles, &particles , interactions, counter);

                                interactions[outsideInteractions[outInterIdx].outPosition] = nullptr;
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
        return 27-index-1;
    }

    int getOppositeInterIndex(const int index) const {
        // ((( (xdiff+3) * 7) + (ydiff+3))) * 7 + zdiff + 3
        return 343-index-1;
    }
};


#endif // FGROUPSEQALGORITHM_HPP
