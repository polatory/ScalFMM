// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FFMMALGORITHMSTARPUGROUP_HPP
#define FFMMALGORITHMSTARPUGROUP_HPP

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"
#include "../Utils/FMemUtils.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FBoolArray.hpp"


#include "../Extensions/FExtendCoordinate.hpp"
#include "../Extensions/FExtendMortonIndex.hpp"

/*
TODO:
scinder multipole/local
*/


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmStarpuGroup
* @brief
* Please read the license
*/
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmStarpuGroup : protected FAssertable{
    /////////////////////////////////////////////////////////////
    // Utils classes
    /////////////////////////////////////////////////////////////


    struct MortonContainer : public FExtendMortonIndex, public FExtendCoordinate {
        ContainerClass container;
    };

    /** This structure holds the data properties needed
      * by a cell/leaf to finish its computation
      */
    struct TransferProperties{
        explicit TransferProperties(const int inIndex = 0, const int inPosition = 0, const int inDataPos = 0)
            : indexWhoNeedsData(inIndex), positionInComputationArray(inPosition), positionInDataArray(inDataPos) {
        }
        // In the group destination, who need the data?
        int indexWhoNeedsData;
        // where to put the data in the array
        int positionInComputationArray;
        // Where to read the data from?
        int positionInDataArray;
    };

    /** The transfer buffer holds many properties
      * it has enough information to create a copy task and
      * a process task
      */
    struct TransferBuffer {
        TransferBuffer() : groupDestination(0), transferBufferLeaf(0) {
        }
        ~TransferBuffer() {
            delete[] transferBufferLeaf;
        }

        // the group who need the data
        int groupDestination;
        // position in the original group
        FVector<int> originalIndexPosition;
        // transfer properties
        FVector<TransferProperties> compuationProperties;
        // where data will be copied
        int indexToStarCopying;
        // memory to copy before compute remotly
        ContainerClass* FRestrict transferBufferLeaf;
    };

    /** A group contains several cells
      * and some properties
      */
    struct Group {
        Group() : cellArray(0), needOther(0), leavesArray(0), transferBufferCell(0), nbCellToReceive(0) {
        }
        ~Group(){
            delete[] cellArray;
            delete[] needOther;
            delete[] leavesArray;
            for(int idx = 0 ; idx < dataToSend.getSize() ; ++idx){
                delete dataToSend[idx];
            }
            delete[] transferBufferCell;
        }

        // Morton index the group start at
        MortonIndex beginIndex;
        // Morton index the group end at
        MortonIndex endIndex;
        // Number of elements in the group, usually GroupSize
        int nbElements;
        // The data of the group
        CellClass* cellArray;
        bool* needOther;
        // Or the leaves data
        MortonContainer* FRestrict leavesArray;
        // Information needed to compute parent child operations
        int indexOfStartInLowerGroups;
        FVector<Group*> lowerGroups;
        // Information needed in case of transfering data needed
        FVector<TransferBuffer*> dataToSend;

        // memory to copy before compute remotly
        CellClass* FRestrict transferBufferCell;
        int nbCellToReceive;
    };

    /////////////////////////////////////////////////////////////
    // Attributes
    /////////////////////////////////////////////////////////////

    OctreeClass* const tree;      //< The octree to work on
    const int OctreeHeight;       //< Height of the tree
    const int BlockSize;          //< Size of the block
    Group**const blockedTree;     //< Current block tree
    int*const blockedPerLevel;    //< Number of block per level
    KernelClass* const kernel;    //< The kernel

public:
    /** The constructor need the octree and the kernel used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernel to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmStarpuGroup(OctreeClass* const inTree, KernelClass* const inKernel,
                             const int inBlockedSize = 25)
        : tree(inTree), OctreeHeight(tree->getHeight()),
          BlockSize(inBlockedSize),
          blockedTree(new Group*[OctreeHeight + 1]) ,
          blockedPerLevel(new int[OctreeHeight + 1]),
          kernel(inKernel){

        fassert(tree, "tree cannot be null", __LINE__, __FILE__);
        fassert(kernel, "kernel cannot be null", __LINE__, __FILE__);

        memset(blockedTree, 0, sizeof(Group*) * (OctreeHeight + 1));
        memset(blockedPerLevel, 0, (OctreeHeight + 1) * sizeof(int));

        FDEBUG(FDebug::Controller << "FFmmAlgorithmStarpuGroup (block size = " << BlockSize <<")\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmStarpuGroup(){
        delete[] blockedTree;
        delete[] blockedPerLevel;
    }

    /////////////////////////////////////////////////////////////
    // Tree to group functions
    /////////////////////////////////////////////////////////////

    /**
      */
    void buildGroups(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

        // Count leaf to allocate and big array
        typename OctreeClass::Iterator* iterArray = 0;
        {
            int leafsNumber = 0;
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            do{
                ++leafsNumber;
            } while(octreeIterator.moveRight());
            iterArray = new typename OctreeClass::Iterator[leafsNumber];
            fassert(iterArray, "iterArray bad alloc", __LINE__, __FILE__);
        }
        FDEBUG( FDebug::Controller << "\tCopy the tree\n"; );
        // Then we start creating the block
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            typename OctreeClass::Iterator avoidGotLeftIterator(octreeIterator);
            for(int idxLevel = 1; idxLevel < OctreeHeight; ++idxLevel){
                // put every thing in the array
                int counterAtLevel = 0;
                do{
                    iterArray[counterAtLevel++] = octreeIterator;
                } while(octreeIterator.moveRight());
                avoidGotLeftIterator.moveDown();
                octreeIterator = avoidGotLeftIterator;
                // find the number of groups
                const int NbGroups = (counterAtLevel + BlockSize - 1) / BlockSize;

                blockedPerLevel[idxLevel] = NbGroups;
                blockedTree[idxLevel] = new Group[NbGroups];

                // copy data to group
                int copyIndex = 0;
                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    const int cellsInThisGroup = FMath::Min(BlockSize, counterAtLevel-copyIndex);

                    blockedTree[idxLevel][idxGroup].nbElements = cellsInThisGroup;
                    blockedTree[idxLevel][idxGroup].cellArray = new CellClass[cellsInThisGroup];
                    blockedTree[idxLevel][idxGroup].needOther = new bool[cellsInThisGroup];

                    for(int idxCell = 0 ; idxCell < cellsInThisGroup ; ++idxCell, ++copyIndex){
                        blockedTree[idxLevel][idxGroup].cellArray[idxCell].setMortonIndex( iterArray[copyIndex].getCurrentGlobalIndex() );
                        blockedTree[idxLevel][idxGroup].cellArray[idxCell].setCoordinate( iterArray[copyIndex].getCurrentGlobalCoordinate() );
                        blockedTree[idxLevel][idxGroup].needOther[idxCell] = false;
                        blockedTree[idxLevel][idxGroup].cellArray[idxCell].intialCopy( iterArray[copyIndex].getCurrentCell() );
                    }

                    blockedTree[idxLevel][idxGroup].beginIndex = blockedTree[idxLevel][idxGroup].cellArray[0].getMortonIndex();
                    blockedTree[idxLevel][idxGroup].endIndex = blockedTree[idxLevel][idxGroup].cellArray[cellsInThisGroup-1].getMortonIndex();
                }
            }
            // leaf level will have the same groups has head cell level
            const int NbGroups = blockedPerLevel[OctreeHeight-1];
            blockedPerLevel[OctreeHeight] = NbGroups;
            blockedTree[OctreeHeight] = new Group[NbGroups];
            int copyIndex = 0;
            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                blockedTree[OctreeHeight][idxGroup].nbElements = blockedTree[OctreeHeight-1][idxGroup].nbElements;
                blockedTree[OctreeHeight][idxGroup].beginIndex = blockedTree[OctreeHeight-1][idxGroup].beginIndex;
                blockedTree[OctreeHeight][idxGroup].endIndex = blockedTree[OctreeHeight-1][idxGroup].endIndex;

                const int NbLeaves = blockedTree[OctreeHeight][idxGroup].nbElements;
                blockedTree[OctreeHeight][idxGroup].leavesArray = new MortonContainer[NbLeaves];
                for(int idxLeaf = 0 ; idxLeaf < NbLeaves ; ++idxLeaf, ++copyIndex){
                    blockedTree[OctreeHeight][idxGroup].leavesArray[idxLeaf].container = *iterArray[copyIndex].getCurrentListSrc();
                    blockedTree[OctreeHeight][idxGroup].leavesArray[idxLeaf].setMortonIndex(iterArray[copyIndex].getCurrentGlobalIndex());
                    blockedTree[OctreeHeight][idxGroup].leavesArray[idxLeaf].setCoordinate(iterArray[copyIndex].getCurrentGlobalCoordinate());
                }
            }
        }
        delete[] iterArray;
        iterArray = 0;
        FDEBUG( FDebug::Controller << "\tPrepare child parent relations\n"; );
        // All block has been created, now find the parent-child relations
        {
            for(int idxLevel = 1; idxLevel < OctreeHeight - 1; ++idxLevel){
                int currentLowerGroup = 0;
                // find the number of groups
                const int NbGroups = blockedPerLevel[idxLevel];
                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    // copy current group
                    blockedTree[idxLevel][idxGroup].lowerGroups.push( &blockedTree[idxLevel+1][currentLowerGroup] );
                    {
                        int startIndex = 0;
                        while( (blockedTree[idxLevel+1][currentLowerGroup].cellArray[startIndex].getMortonIndex()>>3) != blockedTree[idxLevel][idxGroup].beginIndex){
                            ++startIndex;
                        }
                        blockedTree[idxLevel][idxGroup].indexOfStartInLowerGroups = startIndex;
                    }
                    if((blockedTree[idxLevel+1][currentLowerGroup].endIndex>>3) <= blockedTree[idxLevel][idxGroup].endIndex){
                        ++currentLowerGroup;
                    }
                    // copy until too much on the right
                    while(currentLowerGroup < blockedPerLevel[idxLevel+1] && (blockedTree[idxLevel+1][currentLowerGroup].endIndex>>3) <= blockedTree[idxLevel][idxGroup].endIndex) {
                        blockedTree[idxLevel][idxGroup].lowerGroups.push( &blockedTree[idxLevel+1][currentLowerGroup] );
                        ++currentLowerGroup;
                    }
                    // do we have to move backward for the next group
                    if(currentLowerGroup < blockedPerLevel[idxLevel+1] &&
                            (blockedTree[idxLevel+1][currentLowerGroup].beginIndex>>3) <= blockedTree[idxLevel][idxGroup].endIndex){
                        blockedTree[idxLevel][idxGroup].lowerGroups.push( &blockedTree[idxLevel+1][currentLowerGroup] );
                    }
                }
            }
        }
        FDEBUG( FDebug::Controller << "\tPrepare M2L\n"; );
        // We have the parent-child relation, prepare the Multipole copy
        {
            MortonIndex potentialInteraction[189];
            int potentialPosition[189];
            for(int idxLevel = 1; idxLevel < OctreeHeight; ++idxLevel){
                // for all groups
                const int NbGroups = blockedPerLevel[idxLevel];

                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    const int NbCells = blockedTree[idxLevel][idxGroup].nbElements;
                    // for all cells
                    for( int idxCell = 0 ; idxCell < NbCells ; ++idxCell ){
                        // get interaction
                        const int nbInteraction = getInteractionsFromPosition( blockedTree[idxLevel][idxGroup].cellArray[idxCell].getCoordinate(), idxLevel, potentialInteraction, potentialPosition);
                        // test each interaction
                        for( int idxInteraction = 0 ; idxInteraction < nbInteraction ; ++idxInteraction ){
                            // is there an interaction out of the block?
                            if(potentialInteraction[idxInteraction] < blockedTree[idxLevel][idxGroup].beginIndex
                                    || blockedTree[idxLevel][idxGroup].endIndex < potentialInteraction[idxInteraction]){
                                // find the destinatary
                                const int UNFOUND_GROUP = -1;
                                int receiver = idxGroup;
                                if(potentialInteraction[idxInteraction] < blockedTree[idxLevel][idxGroup].beginIndex){
                                    do {
                                        --receiver;
                                    } while(receiver != -1 && potentialInteraction[idxInteraction] < blockedTree[idxLevel][receiver].beginIndex);
                                    // Does someone hold this index?
                                    if( receiver == -1 || blockedTree[idxLevel][receiver].endIndex < potentialInteraction[idxInteraction]){
                                        receiver = UNFOUND_GROUP;
                                    }
                                }
                                else{
                                    do {
                                        ++receiver;
                                    } while(receiver != NbGroups && blockedTree[idxLevel][receiver].endIndex < potentialInteraction[idxInteraction]);
                                    // Does someone hold this index?
                                    if( receiver == NbGroups || potentialInteraction[idxInteraction] < blockedTree[idxLevel][receiver].beginIndex){
                                        receiver = UNFOUND_GROUP;
                                    }
                                }
                                // we found the group
                                if( receiver != UNFOUND_GROUP){
                                    // does he really need it? TODO move to binary search
                                    const int idxInReceiver = findNeigh(blockedTree[idxLevel][receiver].cellArray,
                                                                        blockedTree[idxLevel][receiver].nbElements,
                                                                        potentialInteraction[idxInteraction]);

                                    const bool targetExistInOtherGroup = (idxInReceiver != -1);
                                    if( targetExistInOtherGroup ){
                                        // then prepare for sending data
                                        // find the right transferBuffer if exists
                                        int transferBufferToSend = 0;
                                        while( transferBufferToSend != blockedTree[idxLevel][idxGroup].dataToSend.getSize()
                                               && receiver != blockedTree[idxLevel][idxGroup].dataToSend[transferBufferToSend]->groupDestination ){
                                            ++transferBufferToSend;
                                        }
                                        // if not data are needed from now
                                        if( transferBufferToSend == blockedTree[idxLevel][idxGroup].dataToSend.getSize()){
                                            blockedTree[idxLevel][idxGroup].dataToSend.push( new TransferBuffer );
                                            blockedTree[idxLevel][idxGroup].dataToSend[transferBufferToSend]->groupDestination = receiver;
                                            blockedTree[idxLevel][idxGroup].dataToSend[transferBufferToSend]->originalIndexPosition.push(idxCell);
                                        }
                                        else if(blockedTree[idxLevel][idxGroup].dataToSend[transferBufferToSend]->originalIndexPosition.getSize() &&
                                                blockedTree[idxLevel][idxGroup].dataToSend[transferBufferToSend]->originalIndexPosition.head() != idxCell) {
                                            blockedTree[idxLevel][idxGroup].dataToSend[transferBufferToSend]->originalIndexPosition.push(idxCell);
                                        }
                                        // add to the sending transferBuffer
                                        blockedTree[idxLevel][idxGroup].dataToSend[transferBufferToSend]->compuationProperties.push(TransferProperties(idxInReceiver,342-potentialPosition[idxInteraction],
                                                                                                                  blockedTree[idxLevel][idxGroup].dataToSend[transferBufferToSend]->originalIndexPosition.getSize()-1));
                                        // need other data
                                        blockedTree[idxLevel][idxGroup].needOther[idxCell] = true;
                                    }
                                }
                            }
                        }
                    }
                }
                // pre-compte buffer to copy data
                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    // Copy transferBuffer to send
                    const int nbRemoteInteractions = blockedTree[idxLevel][idxGroup].dataToSend.getSize();
                    for(int idxRemote = 0 ; idxRemote < nbRemoteInteractions ; ++idxRemote){
                        const int cellToCopy = blockedTree[idxLevel][idxGroup].dataToSend[idxRemote]->originalIndexPosition.getSize();
                        const int receiver = blockedTree[idxLevel][idxGroup].dataToSend[idxRemote]->groupDestination;
                        blockedTree[idxLevel][idxGroup].dataToSend[idxRemote]->indexToStarCopying = blockedTree[idxLevel][receiver].nbCellToReceive;
                        // increase index
                        blockedTree[idxLevel][receiver].nbCellToReceive += cellToCopy;
                    }
                }
                // allocate
                FDEBUG( FReal totalNeeded = 0);
                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    FDEBUG( totalNeeded += blockedTree[idxLevel][idxGroup].nbCellToReceive/FReal(NbGroups); );
                    blockedTree[idxLevel][idxGroup].transferBufferCell = new CellClass[blockedTree[idxLevel][idxGroup].nbCellToReceive];
                }
                FDEBUG( FDebug::Controller << "\t\tAverage cells needed by each group at level " << idxLevel << " = " << totalNeeded << "\n"; );
                // init with right index
                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    // Copy transferBuffer to send
                    const int nbRemoteInteractions = blockedTree[idxLevel][idxGroup].dataToSend.getSize();
                    for(int idxRemote = 0 ; idxRemote < nbRemoteInteractions ; ++idxRemote){
                        const int cellToCopy = blockedTree[idxLevel][idxGroup].dataToSend[idxRemote]->originalIndexPosition.getSize();
                        const int receiver = blockedTree[idxLevel][idxGroup].dataToSend[idxRemote]->groupDestination;
                        const int offset = blockedTree[idxLevel][idxGroup].dataToSend[idxRemote]->indexToStarCopying;
                        for(int idxCell = 0 ; idxCell < cellToCopy ; ++idxCell){
                            const int cellPosition = blockedTree[idxLevel][idxGroup].dataToSend[idxRemote]->originalIndexPosition[idxCell];
                            blockedTree[idxLevel][receiver].transferBufferCell[idxCell + offset].setMortonIndex( blockedTree[idxLevel][idxGroup].cellArray[cellPosition].getMortonIndex() );
                            blockedTree[idxLevel][receiver].transferBufferCell[idxCell + offset].setCoordinate( blockedTree[idxLevel][idxGroup].cellArray[cellPosition].getCoordinate() );
                        }
                    }
                }
            }
        }
        FDEBUG( FDebug::Controller << "\tPrepare P2P\n"; );
        // Prepare sending receiving for the P2P
        {
            MortonIndex potentialInteraction[26];
            int potentialPosition[26];

            const int NbGroups = blockedPerLevel[OctreeHeight];

            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                const int NbCells = blockedTree[OctreeHeight][idxGroup].nbElements;
                // for all cells
                for( int idxCell = 0 ; idxCell < NbCells ; ++idxCell ){
                    // get interaction
                    const int nbInteraction = getNeighborsFromPosition( blockedTree[OctreeHeight][idxGroup].leavesArray[idxCell].getCoordinate(), OctreeHeight, potentialInteraction, potentialPosition);
                    // test each interaction
                    for( int idxInteraction = 0 ; idxInteraction < nbInteraction ; ++idxInteraction ){
                        // is there an interaction out of the block?
                        if(potentialInteraction[idxInteraction] < blockedTree[OctreeHeight][idxGroup].beginIndex
                                || blockedTree[OctreeHeight][idxGroup].endIndex < potentialInteraction[idxInteraction]){
                            // find the destinatary
                            const int UNFOUND_GROUP = -1;
                            int receiver = idxGroup;
                            if(potentialInteraction[idxInteraction] < blockedTree[OctreeHeight][idxGroup].beginIndex){
                                do {
                                    --receiver;
                                } while(receiver != -1 && potentialInteraction[idxInteraction] < blockedTree[OctreeHeight][receiver].beginIndex);
                                // Does someone hold this index?
                                if( receiver == -1 || blockedTree[OctreeHeight][receiver].endIndex < potentialInteraction[idxInteraction]){
                                    receiver = UNFOUND_GROUP;
                                }
                            }
                            else{
                                do {
                                    ++receiver;
                                } while(receiver != NbGroups && blockedTree[OctreeHeight][receiver].endIndex < potentialInteraction[idxInteraction]);
                                // Does someone hold this index?
                                if( receiver == NbGroups || potentialInteraction[idxInteraction] < blockedTree[OctreeHeight][receiver].beginIndex){
                                    receiver = UNFOUND_GROUP;
                                }
                            }
                            // we found the group
                            if( receiver != UNFOUND_GROUP ){
                                // does he really need it? TODO move to binary search
                                const int idxInReceiver = findNeigh( blockedTree[OctreeHeight][receiver].leavesArray,
                                                               blockedTree[OctreeHeight][receiver].nbElements,
                                                               potentialInteraction[idxInteraction]);

                                const bool targetExistInOtherGroup = (idxInReceiver != -1);
                                if( targetExistInOtherGroup ){
                                    // then prepare for sending data
                                    // find the right transferBuffer if exists
                                    int transferBufferToSend = 0;
                                    while( transferBufferToSend != blockedTree[OctreeHeight][idxGroup].dataToSend.getSize()
                                           && receiver != blockedTree[OctreeHeight][idxGroup].dataToSend[transferBufferToSend]->groupDestination ){
                                        ++transferBufferToSend;
                                    }
                                    // if not data are needed from now
                                    if( transferBufferToSend == blockedTree[OctreeHeight][idxGroup].dataToSend.getSize()){
                                        blockedTree[OctreeHeight][idxGroup].dataToSend.push( new TransferBuffer );
                                        blockedTree[OctreeHeight][idxGroup].dataToSend[transferBufferToSend]->groupDestination = receiver;
                                        blockedTree[OctreeHeight][idxGroup].dataToSend[transferBufferToSend]->originalIndexPosition.push(idxCell);
                                    }
                                    else if( blockedTree[OctreeHeight][idxGroup].dataToSend[transferBufferToSend]->originalIndexPosition.getSize() &&
                                             blockedTree[OctreeHeight][idxGroup].dataToSend[transferBufferToSend]->originalIndexPosition.head() != idxCell ){
                                        blockedTree[OctreeHeight][idxGroup].dataToSend[transferBufferToSend]->originalIndexPosition.push(idxCell);
                                    }
                                    // add to the sending transferBuffer
                                    blockedTree[OctreeHeight][idxGroup].dataToSend[transferBufferToSend]->compuationProperties.push(TransferProperties(idxInReceiver,26-potentialPosition[idxInteraction],
                                                                                                                                                       blockedTree[OctreeHeight][idxGroup].dataToSend[transferBufferToSend]->originalIndexPosition.getSize()-1));
                                }
                            }
                        }
                    }
                }
                // Copy transferBuffer to send
                const int nbRemoteInteractions = blockedTree[OctreeHeight][idxGroup].dataToSend.getSize();
                for(int idxRemote = 0 ; idxRemote < nbRemoteInteractions ; ++idxRemote){
                    const int leafToCopy = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->originalIndexPosition.getSize();
                    blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->transferBufferLeaf = new ContainerClass[leafToCopy];
                    for(int idxLeaf = 0 ; idxLeaf < leafToCopy ; ++idxLeaf){
                        const int leafPosition = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->originalIndexPosition[idxLeaf];
                        blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->transferBufferLeaf[idxLeaf] = blockedTree[OctreeHeight][idxGroup].leavesArray[leafPosition].container;
                    }
                }
            }
#ifdef SCALFMM_USE_DEBUG
            FReal totalNeeded = 0;
            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                totalNeeded += blockedTree[OctreeHeight][idxGroup].dataToSend.getSize() / FReal(NbGroups);
            }
            FDebug::Controller << "\t\tEach leaves group are related in average to = " << totalNeeded << " other groups\n";
#endif
        }
    }

    void releaseGroups(){
        // copy from group to octree
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            typename OctreeClass::Iterator avoidGotLeftIterator(octreeIterator);
            for(int idxLevel = 1; idxLevel < OctreeHeight - 1; ++idxLevel){
                const int NbGroups = blockedPerLevel[idxLevel];
                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    const int NbCells = blockedTree[idxLevel][idxGroup].nbElements;
                    for( int idxCell = 0 ; idxCell < NbCells ; ++idxCell ){
                        octreeIterator.getCurrentCell()->restoreCopy( &blockedTree[idxLevel][idxGroup].cellArray[idxCell] );
                        octreeIterator.moveRight();
                    }
                }

                avoidGotLeftIterator.moveDown();
                octreeIterator = avoidGotLeftIterator;
            }
            // restore leaf
            const int NbGroups = blockedPerLevel[OctreeHeight-1];
            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                const int NbCells = blockedTree[OctreeHeight-1][idxGroup].nbElements;
                for( int idxCell = 0 ; idxCell < NbCells ; ++idxCell ){
                    // restore cell
                    octreeIterator.getCurrentCell()->restoreCopy( &blockedTree[OctreeHeight-1][idxGroup].cellArray[idxCell] );
                    // restore particles
                    *octreeIterator.getCurrentListTargets() = blockedTree[OctreeHeight][idxGroup].leavesArray[idxCell].container;
                    // goto right
                    octreeIterator.moveRight();
                }
            }
        }
        // release memory
        for(int idxLevel = 0; idxLevel <= OctreeHeight; ++idxLevel){
            delete[] blockedTree[idxLevel];
        }
        memset(blockedTree, 0, sizeof(Group*) * (OctreeHeight + 1));
        memset(blockedPerLevel, 0, (OctreeHeight + 1) * sizeof(int));
    }

    /////////////////////////////////////////////////////////////
    // Execute functions
    /////////////////////////////////////////////////////////////

    void execute(){
        P2M();

        M2M_M2L();

        P2P();

        L2L();

        L2P();
    }

    /////////////////////////////////////////////////////////////

    // The P2P
    void P2P(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic counterTimeRemote);

        // P2P inside leaves
        const int NbGroups = blockedPerLevel[OctreeHeight];
        for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
            exec_P2P( blockedTree[OctreeHeight][idxGroup].leavesArray,
                      blockedTree[OctreeHeight][idxGroup].nbElements);
        }

        FDEBUG(counterTimeRemote.tic());

        // P2P with partial leaves
        for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
            const int nbRemoteInteraction = blockedTree[OctreeHeight][idxGroup].dataToSend.getSize();
            for( int idxRemote = 0 ; idxRemote < nbRemoteInteraction ; ++idxRemote ){
                const int targetGroup = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->groupDestination;
                exec_P2P_remote(blockedTree[OctreeHeight][targetGroup].leavesArray,
                               blockedTree[OctreeHeight][targetGroup].nbElements,
                               blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->compuationProperties,
                               blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->transferBufferLeaf);
            }
        }

        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (P2P + Remote P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\tRemote P2P only = "  << counterTimeRemote.tacAndElapsed() << "s)\n" );
    }

    void exec_P2P(MortonContainer leafs[], const int size){
        ContainerClass* neighbors[27];
        memset(neighbors, 0, sizeof(ContainerClass*) * 27);

        MortonIndex potentialInteraction[26];
        int potentialPosition[26];

        const MortonIndex begin = leafs[0].getMortonIndex();
        const MortonIndex end = leafs[size-1].getMortonIndex();

        for( int idxLeaf = 0 ; idxLeaf < size ; ++idxLeaf ){
            const int nbInteraction = getNeighborsFromPosition( leafs[idxLeaf].getCoordinate(), OctreeHeight, potentialInteraction, potentialPosition);
            int counter = 0;
            for(int idxInteraction = 0 ; idxInteraction < nbInteraction ; ++idxInteraction){
                if( begin <= potentialInteraction[idxInteraction] && potentialInteraction[idxInteraction] <= end){
                    int neighPosition = -1;
                    int offset = 0;
                    if( potentialInteraction[idxInteraction] < leafs[idxLeaf].getMortonIndex()){
                        neighPosition = findNeigh(leafs, idxLeaf, potentialInteraction[idxInteraction]);
                    }
                    else {
                        neighPosition = findNeigh(leafs + idxLeaf + 1, size - idxLeaf - 1, potentialInteraction[idxInteraction]);
                        offset = idxLeaf+1;
                    }
                    if( neighPosition != -1 ){
                        neighbors[ potentialPosition[idxInteraction] ] = &leafs[neighPosition + offset].container;
                        ++counter;
                    }
                }
            }

            kernel->P2P(leafs[idxLeaf].getCoordinate(),&leafs[idxLeaf].container,
                        &leafs[idxLeaf].container, neighbors, counter);
            if( counter ){
                memset(neighbors, 0, sizeof(ContainerClass*) * 27);
            }
        }
    }

    void exec_P2P_remote(MortonContainer leafs[], const int size,
                   FVector<TransferProperties> compuationProperties,
                   ContainerClass transferBuffer[]){
        // TODO factorize for the same leaf
        ContainerClass* neighbors[27];
        memset(neighbors, 0, sizeof(ContainerClass*) * 27);

        for(int idxInteraction = 0 ; idxInteraction < compuationProperties.getSize() ; ++idxInteraction){
            TransferProperties& properties = compuationProperties[idxInteraction];

            neighbors[properties.positionInComputationArray] = &transferBuffer[properties.positionInDataArray];

            kernel->P2PRemote(leafs[properties.indexWhoNeedsData].getCoordinate(),
                              &leafs[properties.indexWhoNeedsData].container,
                              &leafs[properties.indexWhoNeedsData].container,
                              neighbors, 1);

            neighbors[properties.positionInComputationArray] = 0;
        }
    }

    /////////////////////////////////////////////////////////////

    void P2M(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart P2M\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);

        const int NbGroups = blockedPerLevel[OctreeHeight-1];
        for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
            exec_P2M(blockedTree[OctreeHeight-1][idxGroup].cellArray,
                     blockedTree[OctreeHeight][idxGroup].leavesArray,
                     blockedTree[OctreeHeight-1][idxGroup].nbElements);
        }

        FDEBUG( FDebug::Controller << "\tFinished (@P2M (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    void exec_P2M(CellClass multipole[], const MortonContainer particles[], const int size){
        for( int idxLeaf = 0 ; idxLeaf < size ; ++idxLeaf ){
            kernel->P2M( &multipole[idxLeaf], &particles[idxLeaf].container);
        }
    }

    /////////////////////////////////////////////////////////////

    void M2M_M2L(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward and Transfer pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic counterTimeM2M);
        FDEBUG(FTic counterTimeM2L);
        FDEBUG(FTic counterTimeM2LRemote);

        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel){
            FDEBUG(counterTimeM2M.tic());
            {
                const int NbGroups = blockedPerLevel[idxLevel];
                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    // Todo copy the memory
                    exec_M2M(blockedTree[idxLevel][idxGroup].cellArray,
                             blockedTree[idxLevel][idxGroup].nbElements,
                             blockedTree[idxLevel][idxGroup].indexOfStartInLowerGroups,
                             blockedTree[idxLevel][idxGroup].lowerGroups,
                             idxLevel);
                }
            }
            FDEBUG(counterTimeM2M.tac());
            FDEBUG(counterTimeM2L.tic());
            {
                const int lowerLevel = idxLevel + 1;
                const int NbGroups = blockedPerLevel[lowerLevel];

                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    exec_M2L(blockedTree[lowerLevel][idxGroup].cellArray,
                             blockedTree[lowerLevel][idxGroup].needOther,
                             blockedTree[lowerLevel][idxGroup].nbElements,
                             lowerLevel, blockedTree[lowerLevel][idxGroup].beginIndex,
                             blockedTree[lowerLevel][idxGroup].endIndex);


                    const int nbRemoteInteraction = blockedTree[lowerLevel][idxGroup].dataToSend.getSize();
                    for( int idxRemote = 0 ; idxRemote < nbRemoteInteraction ; ++idxRemote ){
                        // copy
                        const int receiver = blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->groupDestination;
                        exec_M2L_copy(blockedTree[lowerLevel][receiver].transferBufferCell,
                                      blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->indexToStarCopying,
                                      blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->originalIndexPosition,
                                      blockedTree[lowerLevel][idxGroup].cellArray);
                    }
                }

                FDEBUG(counterTimeM2LRemote.tic());
                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                        // remote M2L
                        exec_M2L_remote(blockedTree[lowerLevel][idxGroup].cellArray,
                                        blockedTree[lowerLevel][idxGroup].needOther,
                                        blockedTree[lowerLevel][idxGroup].nbElements,
                                        blockedTree[lowerLevel][idxGroup].transferBufferCell,
                                        blockedTree[lowerLevel][idxGroup].nbCellToReceive,
                                        lowerLevel,
                                        blockedTree[lowerLevel][idxGroup].beginIndex,
                                        blockedTree[lowerLevel][idxGroup].endIndex);
                }                
                FDEBUG(counterTimeM2LRemote.tac());
            }
            FDEBUG(counterTimeM2L.tac());
        }
        {
            FDEBUG(counterTimeM2L.tic());

            const int lowerLevel = 2;
            const int NbGroups = blockedPerLevel[lowerLevel];

            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                exec_M2L(blockedTree[lowerLevel][idxGroup].cellArray,
                         blockedTree[lowerLevel][idxGroup].needOther,
                         blockedTree[lowerLevel][idxGroup].nbElements,
                         lowerLevel, blockedTree[lowerLevel][idxGroup].beginIndex,
                         blockedTree[lowerLevel][idxGroup].endIndex);


                const int nbRemoteInteraction = blockedTree[lowerLevel][idxGroup].dataToSend.getSize();
                for( int idxRemote = 0 ; idxRemote < nbRemoteInteraction ; ++idxRemote ){
                    // copy
                    const int receiver = blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->groupDestination;
                    exec_M2L_copy(blockedTree[lowerLevel][receiver].transferBufferCell,
                                  blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->indexToStarCopying,
                                  blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->originalIndexPosition,
                                  blockedTree[lowerLevel][idxGroup].cellArray);
                }
            }

            FDEBUG(counterTimeM2LRemote.tic());
            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    // remote M2L
                    exec_M2L_remote(blockedTree[lowerLevel][idxGroup].cellArray,
                                    blockedTree[lowerLevel][idxGroup].needOther,
                                    blockedTree[lowerLevel][idxGroup].nbElements,
                                    blockedTree[lowerLevel][idxGroup].transferBufferCell,
                                    blockedTree[lowerLevel][idxGroup].nbCellToReceive,
                                    lowerLevel,
                                    blockedTree[lowerLevel][idxGroup].beginIndex,
                                    blockedTree[lowerLevel][idxGroup].endIndex);
            }
            FDEBUG(counterTimeM2LRemote.tac());

            FDEBUG(counterTimeM2L.tac());
        }

        FDEBUG( FDebug::Controller << "\tFinished (@Upward and transfer (M2M + M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\tM2M only = "  << counterTimeM2M.cumulated() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\tM2L all = "  << counterTimeM2L.cumulated() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\tM2L Remote only = "  << counterTimeM2LRemote.cumulated() << "s)\n" );
    }

    void exec_M2L_remote(CellClass multipole_local[],
                         const bool needOther[],
                   const int size,
                   const CellClass transferBuffer[],
                   const int sizeBuffer,
                   const int level,
                   const MortonIndex begin, const MortonIndex end){
        const CellClass* neighbors[343];
        memset(neighbors, 0, sizeof(CellClass*) * 343);

        MortonIndex potentialInteraction[189];
        int potentialPosition[189];

        for( int idxCell = 0 ; idxCell < size ; ++idxCell ){
            if( needOther[idxCell] ){
                const int nbInteraction = getInteractionsFromPosition( multipole_local[idxCell].getCoordinate(), OctreeHeight, potentialInteraction, potentialPosition);
                int counter = 0;

                for(int idxInteraction = 0 ; idxInteraction < nbInteraction ; ++idxInteraction){
                    if( begin <= potentialInteraction[idxInteraction] && potentialInteraction[idxInteraction] <= end){
                        int neighPosition = -1;
                        int offset = 0;
                        if( potentialInteraction[idxInteraction] < multipole_local[idxCell].getMortonIndex()){
                            neighPosition = findNeigh(multipole_local, idxCell, potentialInteraction[idxInteraction]);
                        }
                        else {
                            neighPosition = findNeigh(multipole_local + idxCell + 1, size - idxCell - 1, potentialInteraction[idxInteraction]);
                            offset = idxCell + 1;
                        }
                        if( neighPosition != -1 ){
                            neighbors[ potentialPosition[idxInteraction] ] = &multipole_local[neighPosition + offset];
                            ++counter;
                        }
                    }
                    else {
                        const int neighPosition = findNeigh(transferBuffer, sizeBuffer, potentialInteraction[idxInteraction]);
                        if( neighPosition != -1 ){
                            neighbors[ potentialPosition[idxInteraction] ] = &transferBuffer[neighPosition];
                            ++counter;
                        }
                    }
                }

                if(counter){
                    kernel->M2L( &multipole_local[idxCell] , neighbors, counter, level);
                    memset(neighbors, 0, sizeof(CellClass*) * 343);
                }
            }
        }
    }

    void exec_M2L_copy(CellClass transferBuffer[], const int indexOfStart,
                       const FVector<int>& originalIndexPosition, const CellClass multipole_local[]){
        for(int idxLeaf = 0 ; idxLeaf < originalIndexPosition.getSize() ; ++idxLeaf){
            const int leafPosition = originalIndexPosition[idxLeaf];
            transferBuffer[idxLeaf + indexOfStart].copyUp(&multipole_local[leafPosition]);
        }
    }

    void exec_M2L(CellClass multipole_local[],
                  const bool needOther[],
                  const int size, const int level,
                  const MortonIndex begin, const MortonIndex end){
        const CellClass* neighbors[343];
        memset(neighbors, 0, sizeof(CellClass*) * 343);

        MortonIndex potentialInteraction[189];
        int potentialPosition[189];

        for( int idxCell = 0 ; idxCell < size ; ++idxCell ){
            if( !needOther[idxCell] ){
                const int nbInteraction = getInteractionsFromPosition( multipole_local[idxCell].getCoordinate(), OctreeHeight, potentialInteraction, potentialPosition);
                int counter = 0;

                for(int idxInteraction = 0 ; idxInteraction < nbInteraction ; ++idxInteraction){
                    if( begin <= potentialInteraction[idxInteraction] && potentialInteraction[idxInteraction] <= end){
                        int neighPosition = -1;
                        int offset = 0;
                        if( potentialInteraction[idxInteraction] < multipole_local[idxCell].getMortonIndex()){
                            neighPosition = findNeigh(multipole_local, idxCell, potentialInteraction[idxInteraction]);
                        }
                        else {
                            neighPosition = findNeigh(multipole_local + idxCell + 1, size - idxCell - 1, potentialInteraction[idxInteraction]);
                            offset = idxCell + 1;
                        }
                        if( neighPosition != -1 ){
                            neighbors[ potentialPosition[idxInteraction] ] = &multipole_local[neighPosition + offset];
                            ++counter;
                        }
                    }
                }

                if(counter){
                    kernel->M2L( &multipole_local[idxCell] , neighbors, counter, level);
                    memset(neighbors, 0, sizeof(CellClass*) * 343);
                }
            }
        }
    }

    void exec_M2M(CellClass multipole[], const int size, const int indexOfStart,
                  const FVector<Group*>& childrenGroup, const int level){
        int childIndex = indexOfStart;
        int childGroup = 0;

        const CellClass* children[8];
        for( int idxCell = 0 ; idxCell < size ; ++idxCell ){
            const MortonIndex currentIndex = multipole[idxCell].getMortonIndex();
            memset(children , 0, 8 * sizeof(CellClass*));

            while( currentIndex == (childrenGroup[childGroup]->cellArray[childIndex].getMortonIndex() >> 3) ){
                children[childrenGroup[childGroup]->cellArray[childIndex].getMortonIndex() & 7] = &childrenGroup[childGroup]->cellArray[childIndex];
                ++childIndex;
                if(childIndex == childrenGroup[childGroup]->nbElements){
                    childIndex = 0;
                    ++childGroup;
                    if( childGroup == childrenGroup.getSize() ){
                        break;
                    }
                }
            }
            kernel->M2M(&multipole[idxCell], children, level);
        }
    }

    /////////////////////////////////////////////////////////////

    void L2L(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Downard pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);

        for(int idxLevel = 2 ; idxLevel < OctreeHeight - 1 ; ++idxLevel){
            const int NbGroups = blockedPerLevel[idxLevel];
            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                // Todo copy the memory
                exec_L2L(blockedTree[idxLevel][idxGroup].cellArray,
                         blockedTree[idxLevel][idxGroup].nbElements,
                         blockedTree[idxLevel][idxGroup].indexOfStartInLowerGroups,
                         blockedTree[idxLevel][idxGroup].lowerGroups,
                         idxLevel);
            }
        }

        FDEBUG( FDebug::Controller << "\tFinished (@Downard (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    void exec_L2L(CellClass multipole[], const int size, const int indexOfStart,
                  FVector<Group*>& childrenGroup, const int level){
        int childIndex = indexOfStart;
        int childGroup = 0;

        CellClass* children[8];
        for( int idxCell = 0 ; idxCell < size ; ++idxCell ){
            const MortonIndex currentIndex = multipole[idxCell].getMortonIndex();
            memset( children, 0, 8 * sizeof(CellClass*));

            while(currentIndex == (childrenGroup[childGroup]->cellArray[childIndex].getMortonIndex() >> 3)){
                children[childrenGroup[childGroup]->cellArray[childIndex].getMortonIndex() & 7] = &childrenGroup[childGroup]->cellArray[childIndex];
                ++childIndex;
                if(childIndex == childrenGroup[childGroup]->nbElements){
                    childIndex = 0;
                    ++childGroup;
                    if( childGroup == childrenGroup.getSize() ){
                        break;
                    }
                }
            }

            kernel->L2L(&multipole[idxCell], children, level);
        }
    }

    /////////////////////////////////////////////////////////////

    void L2P(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tL2P\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);

        const int NbGroups = blockedPerLevel[OctreeHeight-1];
        for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
            // Todo copy the memory
            exec_L2P(blockedTree[OctreeHeight-1][idxGroup].cellArray,
                     blockedTree[OctreeHeight][idxGroup].leavesArray,
                     blockedTree[OctreeHeight][idxGroup].nbElements);
        }

        FDEBUG( FDebug::Controller << "\tFinished (@L2P (L2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    void exec_L2P(const CellClass local[], MortonContainer particles[], const int size){
        for( int idxLeaf = 0 ; idxLeaf < size ; ++idxLeaf ){
            kernel->L2P(&local[idxLeaf], &particles[idxLeaf].container);
        }
    }

    /////////////////////////////////////////////////////////////
    // Utils functions
    /////////////////////////////////////////////////////////////

    template <class Object>
    int findNeigh(const Object leafs[], const int size, const MortonIndex indexToFound) const {

        int idxLeft = 0;
        int idxRight = size - 1;

        while( idxLeft <= idxRight ){
            const int idxMiddle = (idxRight + idxLeft) / 2;
            const MortonIndex mortonMiddle = leafs[idxMiddle].getMortonIndex();

            if( mortonMiddle == indexToFound){
                return idxMiddle;
            }

            if( indexToFound < mortonMiddle ){
                idxRight = idxMiddle - 1;
            }
            else{
                idxLeft = idxMiddle + 1;
            }
        }
        return -1;
    }

    int getInteractionsFromPosition(const FTreeCoordinate& workingCell,const int inLevel, MortonIndex inNeighbors[189], int inNeighborsPosition[189]) const{

        // Then take each child of the parent's neighbors if not in directNeighbors
        // Father coordinate
        const FTreeCoordinate parentCell(workingCell.getX()>>1,workingCell.getY()>>1,workingCell.getZ()>>1);

        // Limite at parent level number of box (split by 2 by level)
        const int limite = FMath::pow2(inLevel-1);

        int idxNeighbors = 0;
        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(parentCell.getX() + idxX,0,limite)) continue;

            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(parentCell.getY() + idxY,0,limite)) continue;

                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(parentCell.getZ() + idxZ,0,limite)) continue;

                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){
                        const FTreeCoordinate otherParent(parentCell.getX() + idxX,parentCell.getY() + idxY,parentCell.getZ() + idxZ);
                        const MortonIndex mortonOther = otherParent.getMortonIndex(inLevel-1);

                        // For each child
                        for(int idxCousin = 0 ; idxCousin < 8 ; ++idxCousin){
                            const int xdiff  = ((otherParent.getX()<<1) | ( (idxCousin>>2) & 1)) - workingCell.getX();
                            const int ydiff  = ((otherParent.getY()<<1) | ( (idxCousin>>1) & 1)) - workingCell.getY();
                            const int zdiff  = ((otherParent.getZ()<<1) | (idxCousin&1)) - workingCell.getZ();

                            // Test if it is a direct neighbor
                            if(FMath::Abs(xdiff) > 1 || FMath::Abs(ydiff) > 1 || FMath::Abs(zdiff) > 1){
                                // add to neighbors
                                inNeighborsPosition[idxNeighbors] = ((( (xdiff+3) * 7) + (ydiff+3))) * 7 + zdiff + 3;
                                inNeighbors[idxNeighbors++] = (mortonOther << 3) | idxCousin;
                            }
                        }
                    }
                }
            }
        }

        return idxNeighbors;
    }

    int getNeighborsFromPosition(const FTreeCoordinate& center, const int inLevel, MortonIndex indexes[26], int indexInArray[26]) const{

        const int limite = FMath::pow2(inLevel-1);
        int idxNeig = 0;
        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(center.getX() + idxX,0, limite)) continue;

            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(center.getY() + idxY,0, limite)) continue;

                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(center.getZ() + idxZ,0, limite)) continue;

                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){
                        const FTreeCoordinate other(center.getX() + idxX,center.getY() + idxY,center.getZ() + idxZ);
                        indexes[ idxNeig ] = other.getMortonIndex(this->OctreeHeight - 1);
                        indexInArray[ idxNeig ] = ((idxX+1)*3 + (idxY+1)) * 3 + (idxZ+1);
                        ++idxNeig;
                    }
                }
            }
        }
        return idxNeig;
    }
};

#endif // FFMMALGORITHMSTARPUGROUP_HPP
