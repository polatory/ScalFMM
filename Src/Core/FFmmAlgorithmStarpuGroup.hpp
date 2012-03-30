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

#include <starpu.h>

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
        TransferBuffer() : groupDestination(0) {
        }
        // the group who need the data
        int groupDestination;
        // position in the original group
        FVector<int> originalIndexPosition;
        // transfer properties
        FVector<TransferProperties> compuationProperties;
        // where data will be copied
        int indexToStarCopying;
    };

    /** A group contains several cells
      * and some properties
      */
    struct Group {
        Group() : cellArray(0), needOther(0), leavesArray(0), transferBufferCell(0),
            nbCellToReceive(0), transferBufferLeaf(0), nbLeafToReceive(0) {
            handleCellArray = 0;
            handleLeafArray = 0;
            handleLeafArrayRead = 0;
            handleTransferCell = 0;
            handleTransferLeaf = 0;
        }
        ~Group(){
            delete[] cellArray;
            delete[] needOther;
            delete[] leavesArray;
            for(int idx = 0 ; idx < dataToSend.getSize() ; ++idx){
                delete dataToSend[idx];
            }
            delete[] transferBufferCell;
            delete[] transferBufferLeaf;

            if( handleCellArray != starpu_data_handle_t(0)) starpu_data_unregister(handleCellArray);
            if( handleLeafArray != starpu_data_handle_t(0)) starpu_data_unregister(handleLeafArray);
            if( handleLeafArrayRead != starpu_data_handle_t(0)) starpu_data_unregister(handleLeafArrayRead);
            if( handleTransferCell != starpu_data_handle_t(0)) starpu_data_unregister(handleTransferCell);
            if( handleTransferLeaf != starpu_data_handle_t(0)) starpu_data_unregister(handleTransferLeaf);
        }

        // Morton index the group start at
        MortonIndex beginIndex;
        // Morton index the group end at
        MortonIndex endIndex;
        // Number of elements in the group, usually GroupSize
        int nbElements;
        // The data of the group
        CellClass* FRestrict cellArray;
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

        // memory to copy before compute remotly
        MortonContainer* FRestrict transferBufferLeaf;
        int nbLeafToReceive;

        // Starpu data
        starpu_data_handle_t handleCellArray;
        starpu_data_handle_t handleLeafArray;
        starpu_data_handle_t handleLeafArrayRead;
        starpu_data_handle_t handleTransferCell;
        starpu_data_handle_t handleTransferLeaf;
    };

    //////////////////////////////////////////////////////////////////
    // Init Kernels
    //////////////////////////////////////////////////////////////////

    // Init the fmm kernel (1 per thread)
    void initKernels(){
        globalKernels = new KernelClass*[starpu_worker_get_count()];
        memset(globalKernels, 0, sizeof(KernelClass*) * starpu_worker_get_count());

        for(unsigned int workerid = 0; workerid < starpu_worker_get_count(); ++workerid){
            if( starpu_worker_get_type(workerid) == STARPU_CPU_WORKER ){
                globalKernels[workerid] = new KernelClass(*kernel);
            }
        }
    }

    // Delete kernels
    void releaseKernels(){
        for(unsigned int workerid = 0; workerid < starpu_worker_get_count(); ++workerid){
            delete globalKernels[workerid];
        }
        delete[] globalKernels;
    }

    /////////////////////////////////////////////////////////////
    // Attributes
    /////////////////////////////////////////////////////////////

    OctreeClass* const tree;      //< The octree to work on
    const int OctreeHeight;       //< Height of the tree
    const int BlockSize;          //< Size of the block
    Group**const blockedTree;     //< Current block tree
    int*const blockedPerLevel;    //< Number of block per level
    KernelClass* const kernel;    //< The kernel
    const bool useStarpuPerfModel;//< to know if perf model has to be used

    static const int MaxChild = 9;

    starpu_codelet p2m_cl;
    starpu_codelet p2p_cl;
    starpu_codelet p2p_restore_cl;
    starpu_codelet m2m_cl[MaxChild];
    starpu_codelet m2l_cl;
    starpu_codelet m2l_other_cl;
    starpu_codelet m2l_copy_cl;
    starpu_codelet l2l_cl[MaxChild];
    starpu_codelet l2p_cl;

    starpu_perfmodel p2p_model;
    starpu_perfmodel p2p_restore_model;
    starpu_perfmodel p2m_model;
    starpu_perfmodel m2m_model;
    starpu_perfmodel m2l_model;
    starpu_perfmodel m2l_other_model;
    starpu_perfmodel m2l_copy_model;
    starpu_perfmodel l2l_model;
    starpu_perfmodel l2p_model;

    void initCodelet(){
        // init perf model
        memset(&p2p_model, 0, sizeof(p2p_model));
        p2p_model.type = STARPU_HISTORY_BASED;
        p2p_model.symbol = "P2P";
        memset(&p2p_restore_model, 0, sizeof(p2p_restore_model));
        p2p_restore_model.type = STARPU_HISTORY_BASED;
        p2p_restore_model.symbol = "P2P Restore";
        memset(&p2m_model, 0, sizeof(p2m_model));
        p2m_model.type = STARPU_HISTORY_BASED;
        p2m_model.symbol = "P2M";
        memset(&m2l_model, 0, sizeof(m2l_model));
        m2l_model.type = STARPU_HISTORY_BASED;
        m2l_model.symbol = "M2L";
        memset(&m2l_other_model, 0, sizeof(m2l_other_model));
        m2l_other_model.type = STARPU_HISTORY_BASED;
        m2l_other_model.symbol = "M2L Other";
        memset(&m2l_copy_model, 0, sizeof(m2l_model));
        m2l_copy_model.type = STARPU_HISTORY_BASED;
        m2l_copy_model.symbol = "M2L Copy";
        memset(&l2p_model, 0, sizeof(l2p_model));
        l2p_model.type = STARPU_HISTORY_BASED;
        l2p_model.symbol = "L2P";
        memset(&l2l_model, 0, sizeof(l2l_model));
        l2l_model.type = STARPU_HISTORY_BASED;
        l2l_model.symbol = "L2L";
        memset(&m2m_model, 0, sizeof(m2m_model));
        m2m_model.type = STARPU_HISTORY_BASED;
        m2m_model.symbol = "M2M";

        // P2M
        memset(&p2m_cl, 0, sizeof(p2m_cl));
        p2m_cl.where = STARPU_CPU;
        p2m_cl.cpu_funcs[0] = p2m_cpu;
        p2m_cl.nbuffers = 2;
        p2m_cl.modes[0] = STARPU_W;
        p2m_cl.modes[1] = STARPU_R;
        if(useStarpuPerfModel) p2m_cl.model = &p2m_model;
        // P2P
        memset(&p2p_cl, 0, sizeof(starpu_codelet) );
        p2p_cl.where = STARPU_CPU;
        p2p_cl.cpu_funcs[0] = p2p_cpu;
        p2p_cl.nbuffers = 2;
        p2p_cl.modes[0] = STARPU_RW;
        p2p_cl.modes[1] = STARPU_RW;
        if( useStarpuPerfModel ) p2p_cl.model = &p2p_model;
        // P2P restore
        memset(&p2p_restore_cl, 0, sizeof(starpu_codelet) );
        p2p_restore_cl.where = STARPU_CPU;
        p2p_restore_cl.cpu_funcs[0] = p2p_restore_cpu;
        p2p_restore_cl.nbuffers = 2;
        p2p_restore_cl.modes[0] = STARPU_RW;
        p2p_restore_cl.modes[1] = STARPU_R;
        if( useStarpuPerfModel ) p2p_restore_cl.model = &p2p_restore_model;

        // L2P
        memset(&l2p_cl, 0, sizeof(l2p_cl));
        l2p_cl.where = STARPU_CPU;
        l2p_cl.cpu_funcs[0] = l2p_cpu;
        l2p_cl.nbuffers = 2;
        l2p_cl.modes[0] = STARPU_R;
        l2p_cl.modes[1] = STARPU_RW;
        if(useStarpuPerfModel)  l2p_cl.model = &l2p_model;
        // M2L
        memset(&m2l_cl, 0, sizeof(starpu_codelet) );
        m2l_cl.where = STARPU_CPU;
        m2l_cl.cpu_funcs[0] = m2l_cpu;
        m2l_cl.nbuffers = 1;
        m2l_cl.modes[0] = STARPU_RW;
        if( useStarpuPerfModel ) m2l_cl.model = &m2l_model;
        // M2L other
        memset(&m2l_other_cl, 0, sizeof(starpu_codelet) );
        m2l_other_cl.where = STARPU_CPU;
        m2l_other_cl.cpu_funcs[0] = m2l_other_cpu;
        m2l_other_cl.nbuffers = 2;
        m2l_other_cl.modes[0] = STARPU_RW;
        m2l_other_cl.modes[1] = STARPU_R;
        if( useStarpuPerfModel ) m2l_other_cl.model = &m2l_other_model;
        // M2L copy
        memset(&m2l_copy_cl, 0, sizeof(starpu_codelet) );
        m2l_copy_cl.where = STARPU_CPU;
        m2l_copy_cl.cpu_funcs[0] = m2l_copy_cpu;
        m2l_copy_cl.nbuffers = 2;
        m2l_copy_cl.modes[0] = STARPU_RW;
        m2l_copy_cl.modes[1] = STARPU_R;
        if( useStarpuPerfModel ) m2l_copy_cl.model = &m2l_copy_model;
        // M2M & L2L
        memset(m2m_cl, 0, sizeof(starpu_codelet) * MaxChild);
        memset(l2l_cl, 0, sizeof(starpu_codelet) * MaxChild);
        for( int idxChild = 0 ; idxChild < MaxChild ; ++idxChild){
            m2m_cl[idxChild].where = STARPU_CPU;
            m2m_cl[idxChild].cpu_funcs[0] = m2m_cpu;
            m2m_cl[idxChild].nbuffers = idxChild + 2;
            m2m_cl[idxChild].modes[0] = STARPU_W;
            if( useStarpuPerfModel) m2m_cl[idxChild].model = &m2m_model;

            l2l_cl[idxChild].where = STARPU_CPU;
            l2l_cl[idxChild].cpu_funcs[0] = l2l_cpu;
            l2l_cl[idxChild].nbuffers = idxChild + 2;
            l2l_cl[idxChild].modes[0] = STARPU_R;
            if( useStarpuPerfModel) l2l_cl[idxChild].model = &l2l_model;

            for( int idxMode = 0 ; idxMode <= idxChild ; ++idxMode){
                m2m_cl[idxChild].modes[idxMode+1] = STARPU_R;
                l2l_cl[idxChild].modes[idxMode+1] = STARPU_RW;
            }
        }
    }

public:
    /** The constructor need the octree and the kernel used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernel to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmStarpuGroup(OctreeClass* const inTree, KernelClass* const inKernel,
                             const int inBlockedSize = 25, const bool inUseStarpuPerfModel = false)
        : tree(inTree), OctreeHeight(tree->getHeight()),
          BlockSize(inBlockedSize),
          blockedTree(new Group*[OctreeHeight + 1]) ,
          blockedPerLevel(new int[OctreeHeight + 1]),
          kernel(inKernel),
          useStarpuPerfModel(inUseStarpuPerfModel) {

        fassert(tree, "tree cannot be null", __LINE__, __FILE__);
        fassert(kernel, "kernel cannot be null", __LINE__, __FILE__);

        memset(blockedTree, 0, sizeof(Group*) * (OctreeHeight + 1));
        memset(blockedPerLevel, 0, (OctreeHeight + 1) * sizeof(int));

        FDEBUG(FDebug::Controller << "FFmmAlgorithmStarpuGroup (Block size = " << BlockSize <<")\n");
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

        // star starpu
        starpu_init(NULL);
        FDEBUG(FDebug::Controller << "Start starpu runtime, Nb Workers = " << starpu_worker_get_count() << "\n");

        // create codelet
        initCodelet();

        // create kernel for all thread
        initKernels();

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
                FDEBUG( FDebug::Controller << "\t\tAt level " << idxLevel << " there are " << NbGroups << " groups\n"; );
                blockedPerLevel[idxLevel] = NbGroups;
                blockedTree[idxLevel] = new Group[NbGroups];

                // copy data to group
                int copyIndex = 0;
                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    const int cellsInThisGroup = FMath::Min(BlockSize, counterAtLevel-copyIndex);

                    blockedTree[idxLevel][idxGroup].nbElements = cellsInThisGroup;
                    blockedTree[idxLevel][idxGroup].cellArray = new CellClass[cellsInThisGroup];
                    blockedTree[idxLevel][idxGroup].needOther = new bool[cellsInThisGroup];

                    // starpu
                    starpu_vector_data_register(&blockedTree[idxLevel][idxGroup].handleCellArray, 0,
                                                (uintptr_t)blockedTree[idxLevel][idxGroup].cellArray,
                                                blockedTree[idxLevel][idxGroup].nbElements, sizeof(CellClass));

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

                // starpu
                starpu_vector_data_register(&blockedTree[OctreeHeight][idxGroup].handleLeafArray, 0,
                                            (uintptr_t)blockedTree[OctreeHeight][idxGroup].leavesArray,
                                            NbLeaves, sizeof(MortonContainer));
                starpu_vector_data_register(&blockedTree[OctreeHeight][idxGroup].handleLeafArrayRead, 0,
                                            (uintptr_t)blockedTree[OctreeHeight][idxGroup].leavesArray,
                                            NbLeaves, sizeof(MortonContainer));

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
                FDEBUG( FReal totalDependencies = 0 );
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
                    FDEBUG( totalDependencies += blockedTree[idxLevel][idxGroup].lowerGroups.getSize()/FReal(NbGroups) );
                }
                FDEBUG( FDebug::Controller << "\t\tAt level " << idxLevel << " average parent-child dependencies " << totalDependencies << "\n"; );
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
                    // starpu
                    starpu_vector_data_register(&blockedTree[idxLevel][idxGroup].handleTransferCell, 0,
                                                (uintptr_t)blockedTree[idxLevel][idxGroup].transferBufferCell,
                                                blockedTree[idxLevel][idxGroup].nbCellToReceive, sizeof(CellClass));
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
            }
            // pre-compte buffer to copy data
            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                // Copy transferBuffer to send
                const int nbRemoteInteractions = blockedTree[OctreeHeight][idxGroup].dataToSend.getSize();
                for(int idxRemote = 0 ; idxRemote < nbRemoteInteractions ; ++idxRemote){
                    const int leafToCopy = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->originalIndexPosition.getSize();
                    const int receiver = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->groupDestination;
                    blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->indexToStarCopying = blockedTree[OctreeHeight][receiver].nbLeafToReceive;
                    // increase index
                    blockedTree[OctreeHeight][receiver].nbLeafToReceive += leafToCopy;
                }
            }
            // allocate needed memory
            FDEBUG( FReal totalNeeded = 0 );
            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                FDEBUG( totalNeeded += blockedTree[OctreeHeight][idxGroup].nbLeafToReceive/FReal(NbGroups); );
                blockedTree[OctreeHeight][idxGroup].transferBufferLeaf = new MortonContainer[blockedTree[OctreeHeight][idxGroup].nbLeafToReceive];
                // starpu
                starpu_vector_data_register(&blockedTree[OctreeHeight][idxGroup].handleTransferLeaf, 0,
                                            (uintptr_t)blockedTree[OctreeHeight][idxGroup].transferBufferLeaf,
                                            blockedTree[OctreeHeight][idxGroup].nbLeafToReceive, sizeof(MortonContainer));
            }
            FDEBUG( FDebug::Controller << "\t\tAverage leaves needed by each group = " << totalNeeded << "\n"; );
            // copy data
            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                // Copy transferBuffer to send
                const int nbRemoteInteractions = blockedTree[OctreeHeight][idxGroup].dataToSend.getSize();
                for(int idxRemote = 0 ; idxRemote < nbRemoteInteractions ; ++idxRemote){
                    const int leafToCopy = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->originalIndexPosition.getSize();
                    const int receiver = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->groupDestination;
                    const int offset = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->indexToStarCopying;
                    for(int idxLeaf = 0 ; idxLeaf < leafToCopy ; ++idxLeaf){
                        const int leafPosition = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->originalIndexPosition[idxLeaf];
                        blockedTree[OctreeHeight][receiver].transferBufferLeaf[idxLeaf + offset].setMortonIndex( blockedTree[OctreeHeight][idxGroup].leavesArray[leafPosition].getMortonIndex() );
                        blockedTree[OctreeHeight][receiver].transferBufferLeaf[idxLeaf + offset].setCoordinate( blockedTree[OctreeHeight][idxGroup].leavesArray[leafPosition].getCoordinate() );
                        blockedTree[OctreeHeight][receiver].transferBufferLeaf[idxLeaf + offset].container = blockedTree[OctreeHeight][idxGroup].leavesArray[leafPosition].container;
                    }
                }
            }
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

        // release kernels
        releaseKernels();

        // shutdown starpu
        FDEBUG( FDebug::Controller << "Shutdown starpu\n" );
        starpu_shutdown();
    }

    /////////////////////////////////////////////////////////////
    // Execute functions
    /////////////////////////////////////////////////////////////

    void execute(){
        P2P_P2M();

        M2M_M2L();

        L2L();

        L2P();

        FDEBUG( FDebug::Controller << "Wait task to be finished...\n" );
        starpu_task_wait_for_all();
    }

    /////////////////////////////////////////////////////////////

    // The P2P
    void P2P_P2M(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass (P2P + P2M)\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTimeP2PP2M);
        FDEBUG(FTic counterTime);
        FDEBUG(FTic counterTimeCopy);

        // P2P inside leaves
        const int NbGroups = blockedPerLevel[OctreeHeight];
        for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
            /*exec_P2M(blockedTree[OctreeHeight-1][idxGroup].cellArray,
                     blockedTree[OctreeHeight][idxGroup].leavesArray,
                     blockedTree[OctreeHeight-1][idxGroup].nbElements);*/
            starpu_insert_task( &p2m_cl, STARPU_W, blockedTree[OctreeHeight-1][idxGroup].handleCellArray,
                                STARPU_R, blockedTree[OctreeHeight][idxGroup].handleLeafArrayRead, 0);

            /*exec_P2P( blockedTree[OctreeHeight][idxGroup].leavesArray,
                      blockedTree[OctreeHeight][idxGroup].nbElements,
                      blockedTree[OctreeHeight][idxGroup].transferBufferLeaf,
                      blockedTree[OctreeHeight][idxGroup].nbLeafToReceive);*/

            starpu_insert_task( &p2p_cl, STARPU_VALUE, &OctreeHeight, sizeof(OctreeHeight),
                                STARPU_RW, blockedTree[OctreeHeight][idxGroup].handleLeafArray,
                                STARPU_RW, blockedTree[OctreeHeight][idxGroup].handleTransferLeaf, 0);
        }

        FDEBUG(counterTimeP2PP2M.tac());
        FDEBUG(counterTimeCopy.tic());

        // P2P restore
        for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
            const int nbRemoteInteraction = blockedTree[OctreeHeight][idxGroup].dataToSend.getSize();
            for( int idxRemote = 0 ; idxRemote < nbRemoteInteraction ; ++idxRemote ){
                const int receiver = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->groupDestination;
                const int offset = blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->indexToStarCopying;

                /*exec_P2P_restore(blockedTree[OctreeHeight][idxGroup].leavesArray,
                               blockedTree[OctreeHeight][idxGroup].nbElements,
                               blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->originalIndexPosition,
                                 &blockedTree[OctreeHeight][receiver].transferBufferLeaf[offset]);*/
                FVector<int>* indexPosition = &blockedTree[OctreeHeight][idxGroup].dataToSend[idxRemote]->originalIndexPosition;
                starpu_insert_task( &p2p_restore_cl,
                                    STARPU_VALUE, &indexPosition, sizeof(indexPosition),
                                    STARPU_VALUE, &offset, sizeof(offset),
                                    STARPU_RW, blockedTree[OctreeHeight][idxGroup].handleLeafArray,
                                    STARPU_R, blockedTree[OctreeHeight][receiver].handleTransferLeaf, 0);
            }
        }

        FDEBUG( FDebug::Controller << "\tFinished (@Bottom Pass (P2M + P2P + Copy P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\tP2P + P2M = "  << counterTimeP2PP2M.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\tRestore P2P leaves only = "  << counterTimeCopy.tacAndElapsed() << "s)\n" );
    }

    void exec_P2P(MortonContainer leafs[], const int size, MortonContainer otherleafs[], const int othersize){
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
                else if( othersize ){
                    const int neighPosition = findNeigh(otherleafs, othersize, potentialInteraction[idxInteraction]);
                    if( neighPosition != -1 ){
                        neighbors[ potentialPosition[idxInteraction] ] = &otherleafs[neighPosition].container;
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

    void exec_P2P_restore(MortonContainer leafs[], const int size,
                   FVector<int> originalPosition,
                   const MortonContainer transferBuffer[]){

        for(int idxLeaf = 0 ; idxLeaf < originalPosition.getSize() ; ++idxLeaf){
            leafs[originalPosition[idxLeaf]].container.reduce( &transferBuffer[idxLeaf].container );
        }
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
                    /*exec_M2M(blockedTree[idxLevel][idxGroup].cellArray,
                             blockedTree[idxLevel][idxGroup].nbElements,
                             blockedTree[idxLevel][idxGroup].indexOfStartInLowerGroups,
                             blockedTree[idxLevel][idxGroup].lowerGroups,
                             idxLevel);*/
                    struct starpu_task* const task = starpu_task_create();
                    // buffer 0 is current leaf
                    task->handles[0] = blockedTree[idxLevel][idxGroup].handleCellArray;

                    const int nbLowerGroups = blockedTree[idxLevel][idxGroup].lowerGroups.getSize();
                    for(int idxLower = 0 ; idxLower < nbLowerGroups ; ++idxLower){
                        task->handles[idxLower + 1] = blockedTree[idxLevel][idxGroup].lowerGroups[idxLower]->handleCellArray;
                    }
                    // put the right codelet
                    task->cl = &m2m_cl[nbLowerGroups-1];
                    // put args values
                    char *arg_buffer;
                    size_t arg_buffer_size;
                    starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
                            STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                            STARPU_VALUE, &blockedTree[idxLevel][idxGroup].indexOfStartInLowerGroups, sizeof(blockedTree[idxLevel][idxGroup].indexOfStartInLowerGroups),
                            STARPU_VALUE, &nbLowerGroups, sizeof(nbLowerGroups),
                            0);
                    task->cl_arg = arg_buffer;
                    task->cl_arg_size = arg_buffer_size;

                    // submit task
                    starpu_task_submit(task);
                }
            }
            FDEBUG(counterTimeM2M.tac());

            FDEBUG(counterTimeM2L.tic());
            {
                const int lowerLevel = idxLevel + 1;
                const int NbGroups = blockedPerLevel[lowerLevel];

                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    /*exec_M2L(blockedTree[lowerLevel][idxGroup].cellArray,
                             blockedTree[lowerLevel][idxGroup].needOther,
                             blockedTree[lowerLevel][idxGroup].nbElements,
                             lowerLevel, blockedTree[lowerLevel][idxGroup].beginIndex,
                             blockedTree[lowerLevel][idxGroup].endIndex);*/
                    const MortonIndex begin = blockedTree[lowerLevel][idxGroup].beginIndex;
                    const MortonIndex end = blockedTree[lowerLevel][idxGroup].endIndex;
                    bool*const needOther = blockedTree[lowerLevel][idxGroup].needOther;
                    starpu_insert_task( &m2l_cl,
                                        STARPU_VALUE, &needOther, sizeof(needOther),
                                        STARPU_VALUE, &lowerLevel, sizeof(lowerLevel),
                                        STARPU_VALUE, &begin, sizeof(begin),
                                        STARPU_VALUE, &end, sizeof(end),
                                        STARPU_RW, blockedTree[lowerLevel][idxGroup].handleCellArray, 0);


                    const int nbRemoteInteraction = blockedTree[lowerLevel][idxGroup].dataToSend.getSize();
                    for( int idxRemote = 0 ; idxRemote < nbRemoteInteraction ; ++idxRemote ){
                        // copy
                        const int receiver = blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->groupDestination;
                        const int indexStart = blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->indexToStarCopying;
                        FVector<int>* originalPosition = &blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->originalIndexPosition;
                        /*exec_M2L_copy(blockedTree[lowerLevel][receiver].transferBufferCell,
                                      blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->indexToStarCopying,
                                      blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->originalIndexPosition,
                                      blockedTree[lowerLevel][idxGroup].cellArray);*/
                        starpu_insert_task( &m2l_copy_cl,
                                            STARPU_VALUE, &indexStart, sizeof(indexStart),
                                            STARPU_VALUE, &originalPosition, sizeof(originalPosition),
                                            STARPU_RW, blockedTree[lowerLevel][receiver].handleTransferCell,
                                            STARPU_R, blockedTree[lowerLevel][idxGroup].handleCellArray, 0);
                    }
                }

                FDEBUG(counterTimeM2LRemote.tic());
                for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                        // remote M2L
                        /*exec_M2L_remote(blockedTree[lowerLevel][idxGroup].cellArray,
                                        blockedTree[lowerLevel][idxGroup].needOther,
                                        blockedTree[lowerLevel][idxGroup].nbElements,
                                        blockedTree[lowerLevel][idxGroup].transferBufferCell,
                                        blockedTree[lowerLevel][idxGroup].nbCellToReceive,
                                        lowerLevel,
                                        blockedTree[lowerLevel][idxGroup].beginIndex,
                                        blockedTree[lowerLevel][idxGroup].endIndex);*/
                    const MortonIndex begin = blockedTree[lowerLevel][idxGroup].beginIndex;
                    const MortonIndex end = blockedTree[lowerLevel][idxGroup].endIndex;
                    bool*const needOther = blockedTree[lowerLevel][idxGroup].needOther;
                    starpu_insert_task( &m2l_other_cl,
                                        STARPU_VALUE, &needOther, sizeof(needOther),
                                        STARPU_VALUE, &lowerLevel, sizeof(lowerLevel),
                                        STARPU_VALUE, &begin, sizeof(begin),
                                        STARPU_VALUE, &end, sizeof(end),
                                        STARPU_RW, blockedTree[lowerLevel][idxGroup].handleCellArray,
                                        STARPU_R, blockedTree[lowerLevel][idxGroup].handleTransferCell, 0);
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
                /*exec_M2L(blockedTree[lowerLevel][idxGroup].cellArray,
                         blockedTree[lowerLevel][idxGroup].needOther,
                         blockedTree[lowerLevel][idxGroup].nbElements,
                         lowerLevel, blockedTree[lowerLevel][idxGroup].beginIndex,
                         blockedTree[lowerLevel][idxGroup].endIndex);*/

                const MortonIndex begin = blockedTree[lowerLevel][idxGroup].beginIndex;
                const MortonIndex end = blockedTree[lowerLevel][idxGroup].endIndex;
                bool*const needOther = blockedTree[lowerLevel][idxGroup].needOther;
                starpu_insert_task( &m2l_cl,
                                    STARPU_VALUE, &needOther, sizeof(needOther),
                                    STARPU_VALUE, &lowerLevel, sizeof(lowerLevel),
                                    STARPU_VALUE, &begin, sizeof(begin),
                                    STARPU_VALUE, &end, sizeof(end),
                                    STARPU_RW, blockedTree[lowerLevel][idxGroup].handleCellArray, 0);


                const int nbRemoteInteraction = blockedTree[lowerLevel][idxGroup].dataToSend.getSize();
                for( int idxRemote = 0 ; idxRemote < nbRemoteInteraction ; ++idxRemote ){
                    // copy
                    const int receiver = blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->groupDestination;
                    const int indexStart = blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->indexToStarCopying;
                    FVector<int>* originalPosition = &blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->originalIndexPosition;
                    /*exec_M2L_copy(blockedTree[lowerLevel][receiver].transferBufferCell,
                                  blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->indexToStarCopying,
                                  blockedTree[lowerLevel][idxGroup].dataToSend[idxRemote]->originalIndexPosition,
                                  blockedTree[lowerLevel][idxGroup].cellArray);*/
                    starpu_insert_task( &m2l_copy_cl,
                                        STARPU_VALUE, &indexStart, sizeof(indexStart),
                                        STARPU_VALUE, &originalPosition, sizeof(originalPosition),
                                        STARPU_RW, blockedTree[lowerLevel][receiver].handleTransferCell,
                                        STARPU_R, blockedTree[lowerLevel][idxGroup].handleCellArray, 0);
                }

            }

            FDEBUG(counterTimeM2LRemote.tic());
            for( int idxGroup = 0 ; idxGroup < NbGroups ; ++idxGroup ){
                    // remote M2L
                    /*exec_M2L_remote(blockedTree[lowerLevel][idxGroup].cellArray,
                                    blockedTree[lowerLevel][idxGroup].needOther,
                                    blockedTree[lowerLevel][idxGroup].nbElements,
                                    blockedTree[lowerLevel][idxGroup].transferBufferCell,
                                    blockedTree[lowerLevel][idxGroup].nbCellToReceive,
                                    lowerLevel,
                                    blockedTree[lowerLevel][idxGroup].beginIndex,
                                    blockedTree[lowerLevel][idxGroup].endIndex);*/

                    const MortonIndex begin = blockedTree[lowerLevel][idxGroup].beginIndex;
                    const MortonIndex end = blockedTree[lowerLevel][idxGroup].endIndex;
                    bool*const needOther = blockedTree[lowerLevel][idxGroup].needOther;
                    starpu_insert_task( &m2l_other_cl,
                                        STARPU_VALUE, &needOther, sizeof(needOther),
                                        STARPU_VALUE, &lowerLevel, sizeof(lowerLevel),
                                        STARPU_VALUE, &begin, sizeof(begin),
                                        STARPU_VALUE, &end, sizeof(end),
                                        STARPU_RW, blockedTree[lowerLevel][idxGroup].handleCellArray,
                                        STARPU_R, blockedTree[lowerLevel][idxGroup].handleTransferCell, 0);
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
                const int nbInteraction = getInteractionsFromPosition( multipole_local[idxCell].getCoordinate(), level, potentialInteraction, potentialPosition);
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
                const int nbInteraction = getInteractionsFromPosition( multipole_local[idxCell].getCoordinate(), level, potentialInteraction, potentialPosition);
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
                /*exec_L2L(blockedTree[idxLevel][idxGroup].cellArray,
                         blockedTree[idxLevel][idxGroup].nbElements,
                         blockedTree[idxLevel][idxGroup].indexOfStartInLowerGroups,
                         blockedTree[idxLevel][idxGroup].lowerGroups,
                         idxLevel);*/
                struct starpu_task* const task = starpu_task_create();
                // buffer 0 is current leaf
                task->handles[0] = blockedTree[idxLevel][idxGroup].handleCellArray;

                const int nbLowerGroups = blockedTree[idxLevel][idxGroup].lowerGroups.getSize();
                for(int idxLower = 0 ; idxLower < nbLowerGroups ; ++idxLower){
                    task->handles[idxLower + 1] = blockedTree[idxLevel][idxGroup].lowerGroups[idxLower]->handleCellArray;
                }
                // put the right codelet
                task->cl = &l2l_cl[nbLowerGroups-1];

                // put args values
                char *arg_buffer;
                size_t arg_buffer_size;
                const int indexOfStart = blockedTree[idxLevel][idxGroup].indexOfStartInLowerGroups;
                starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
                        STARPU_VALUE, &idxLevel, sizeof(idxLevel),
                        STARPU_VALUE, &indexOfStart, sizeof(indexOfStart),
                        STARPU_VALUE, &nbLowerGroups, sizeof(nbLowerGroups),
                        0);
                task->cl_arg = arg_buffer;
                task->cl_arg_size = arg_buffer_size;

                // submit task
                starpu_task_submit(task);
            }

        }

        FDEBUG( FDebug::Controller << "\tFinished (@Downard (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    void exec_L2L(const CellClass local[], const int size, const int indexOfStart,
                  FVector<Group*>& childrenGroup, const int level){
        int childIndex = indexOfStart;
        int childGroup = 0;

        CellClass* children[8];
        for( int idxCell = 0 ; idxCell < size ; ++idxCell ){
            const MortonIndex currentIndex = local[idxCell].getMortonIndex();
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

            kernel->L2L(&local[idxCell], children, level);
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
            /*exec_L2P(blockedTree[OctreeHeight-1][idxGroup].cellArray,
                     blockedTree[OctreeHeight][idxGroup].leavesArray,
                     blockedTree[OctreeHeight][idxGroup].nbElements);*/
            starpu_insert_task( &l2p_cl, STARPU_R, blockedTree[OctreeHeight-1][idxGroup].handleCellArray,
                                STARPU_RW, blockedTree[OctreeHeight][idxGroup].handleLeafArray, 0);
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
    static int findNeigh(const Object leafs[], const int size, const MortonIndex indexToFound) {

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

    static int getInteractionsFromPosition(const FTreeCoordinate& workingCell,const int inLevel, MortonIndex inNeighbors[189], int inNeighborsPosition[189]) {

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

    static int getNeighborsFromPosition(const FTreeCoordinate& center, const int inLevel, MortonIndex indexes[26], int indexInArray[26]) {

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
                        indexes[ idxNeig ] = other.getMortonIndex(inLevel-1);
                        indexInArray[ idxNeig ] = ((idxX+1)*3 + (idxY+1)) * 3 + (idxZ+1);
                        ++idxNeig;
                    }
                }
            }
        }
        return idxNeig;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Callback
    /////////////////////////////////////////////////////////////////////////////
    // The current solution uses a global array of kernel
    static KernelClass** globalKernels;

    // P2M
    static void p2m_cpu(void *descr[], void *){
        CellClass*const multipole = (CellClass*)STARPU_VECTOR_GET_PTR(descr[0]);
        const MortonContainer*const particles = (const MortonContainer*)STARPU_VECTOR_GET_PTR(descr[1]);
        const int size = STARPU_VECTOR_GET_NX(descr[0]);

        for( int idxLeaf = 0 ; idxLeaf < size ; ++idxLeaf ){
            globalKernels[starpu_worker_get_id()]->P2M( &multipole[idxLeaf], &particles[idxLeaf].container);
        }
    }

    // M2M
    static void m2m_cpu(void *descr[], void *cl_arg){
        int level;
        int indexOfStart;
        int nbLowerGroup;
        starpu_codelet_unpack_args(cl_arg, &level, &indexOfStart, &nbLowerGroup);

        CellClass*const multipole = (CellClass*)STARPU_VECTOR_GET_PTR(descr[0]);
        const int size = STARPU_VECTOR_GET_NX(descr[0]);

        const CellClass* child_multipole = (CellClass*)STARPU_VECTOR_GET_PTR(descr[1]);
        int child_size = STARPU_VECTOR_GET_NX(descr[1]);

        int childIndex = indexOfStart;
        int childGroup = 0;

        const CellClass* children[8];
        for( int idxCell = 0 ; idxCell < size ; ++idxCell ){
            const MortonIndex currentIndex = multipole[idxCell].getMortonIndex();
            memset(children , 0, 8 * sizeof(CellClass*));

            while( currentIndex == (child_multipole[childIndex].getMortonIndex() >> 3) ){
                children[child_multipole[childIndex].getMortonIndex() & 7] = &child_multipole[childIndex];
                ++childIndex;
                if(childIndex == child_size){
                    childIndex = 0;
                    ++childGroup;
                    if( childGroup == nbLowerGroup ){
                        break;
                    }
                    child_multipole = (const CellClass*)STARPU_VECTOR_GET_PTR(descr[childGroup+1]);
                    child_size = STARPU_VECTOR_GET_NX(descr[childGroup+1]);
                }
            }
            globalKernels[starpu_worker_get_id()]->M2M(&multipole[idxCell], children, level);
        }
    }

    // M2L
    static void m2l_cpu(void *descr[], void *cl_arg)
    {
        bool* needOther;
        int level;
        MortonIndex begin;
        MortonIndex end;
        starpu_codelet_unpack_args(cl_arg, &needOther, &level, &begin, &end);

        CellClass*const multipole_local = (CellClass*)STARPU_VECTOR_GET_PTR(descr[0]);
        const int size = STARPU_VECTOR_GET_NX(descr[0]);


        const CellClass* neighbors[343];
        memset(neighbors, 0, sizeof(CellClass*) * 343);

        MortonIndex potentialInteraction[189];
        int potentialPosition[189];

        for( int idxCell = 0 ; idxCell < size ; ++idxCell ){
            if( !needOther[idxCell] ){
                const int nbInteraction = getInteractionsFromPosition( multipole_local[idxCell].getCoordinate(), level, potentialInteraction, potentialPosition);
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
                    globalKernels[starpu_worker_get_id()]->M2L( &multipole_local[idxCell] , neighbors, counter, level);
                    memset(neighbors, 0, sizeof(CellClass*) * 343);
                }
            }
        }
    }

    // M2L other
    static void m2l_other_cpu(void *descr[], void *cl_arg)
    {

        bool* needOther;
        int level;
        MortonIndex begin;
        MortonIndex end;
        starpu_codelet_unpack_args(cl_arg, &needOther, &level, &begin, &end);

        CellClass*const multipole_local = (CellClass*)STARPU_VECTOR_GET_PTR(descr[0]);
        const int size = STARPU_VECTOR_GET_NX(descr[0]);

        const CellClass*const transferBuffer = (CellClass*)STARPU_VECTOR_GET_PTR(descr[1]);
        const int sizeBuffer = STARPU_VECTOR_GET_NX(descr[1]);

        const CellClass* neighbors[343];
        memset(neighbors, 0, sizeof(CellClass*) * 343);

        MortonIndex potentialInteraction[189];
        int potentialPosition[189];

        for( int idxCell = 0 ; idxCell < size ; ++idxCell ){
            if( needOther[idxCell] ){
                const int nbInteraction = getInteractionsFromPosition( multipole_local[idxCell].getCoordinate(), level, potentialInteraction, potentialPosition);
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
                    globalKernels[starpu_worker_get_id()]->M2L( &multipole_local[idxCell] , neighbors, counter, level);
                    memset(neighbors, 0, sizeof(CellClass*) * 343);
                }
            }
        }
    }

    // M2L copy
    static void m2l_copy_cpu(void *descr[], void *cl_arg)
    {
        int indexOfStart;
        FVector<int>* originalIndexPosition;
        starpu_codelet_unpack_args(cl_arg, &indexOfStart, &originalIndexPosition);

        CellClass*const transferBuffer = (CellClass*)STARPU_VECTOR_GET_PTR(descr[0]);
        const CellClass*const multipole_local = ((const CellClass*)STARPU_VECTOR_GET_PTR(descr[1]));

        for(int idxLeaf = 0 ; idxLeaf < (*originalIndexPosition).getSize() ; ++idxLeaf){
            const int leafPosition = (*originalIndexPosition)[idxLeaf];
            transferBuffer[idxLeaf + indexOfStart].copyUp(&multipole_local[leafPosition]);
        }
    }

    // L2L
    static void l2l_cpu(void *descr[], void * cl_arg) {
        int level;
        int indexOfStart;
        int nbLowerGroup;
        starpu_codelet_unpack_args(cl_arg, &level, &indexOfStart, &nbLowerGroup);

        const CellClass*const local = (CellClass*)STARPU_VECTOR_GET_PTR(descr[0]);
        const int size = STARPU_VECTOR_GET_NX(descr[0]);

        CellClass* child_local = (CellClass*)STARPU_VECTOR_GET_PTR(descr[1]);
        int child_size = STARPU_VECTOR_GET_NX(descr[1]);

        int childIndex = indexOfStart;
        int childGroup = 0;

        CellClass* children[8];
        for( int idxCell = 0 ; idxCell < size ; ++idxCell ){
            const MortonIndex currentIndex = local[idxCell].getMortonIndex();
            memset(children , 0, 8 * sizeof(CellClass*));

            while( currentIndex == (child_local[childIndex].getMortonIndex() >> 3) ){
                children[child_local[childIndex].getMortonIndex() & 7] = &child_local[childIndex];
                ++childIndex;
                if(childIndex == child_size){
                    childIndex = 0;
                    ++childGroup;
                    if( childGroup == nbLowerGroup ){
                        break;
                    }
                    child_local = (CellClass*)STARPU_VECTOR_GET_PTR(descr[childGroup+1]);
                    child_size = STARPU_VECTOR_GET_NX(descr[childGroup+1]);
                }
            }
            globalKernels[starpu_worker_get_id()]->L2L(&local[idxCell], children, level);
        }

    }

    // L2L
    static void l2p_cpu(void *descr[], void *) {
        const CellClass*const local = (const CellClass*)STARPU_VECTOR_GET_PTR(descr[0]);
        MortonContainer*const particles = (MortonContainer*)STARPU_VECTOR_GET_PTR(descr[1]);
        const int size = STARPU_VECTOR_GET_NX(descr[0]);

        for( int idxLeaf = 0 ; idxLeaf < size ; ++idxLeaf ){
            globalKernels[starpu_worker_get_id()]->L2P( &local[idxLeaf], &particles[idxLeaf].container);
        }
    }

    // P2P restore
    static void p2p_restore_cpu(void *descr[], void* cl_arg){
        int offset;
        FVector<int>* originalPosition;
        starpu_codelet_unpack_args(cl_arg, &originalPosition, &offset);

        MortonContainer*const leafs = (MortonContainer*)STARPU_VECTOR_GET_PTR(descr[0]);
        const MortonContainer*const transferBuffer = ((const MortonContainer*)STARPU_VECTOR_GET_PTR(descr[1])) + offset;

        for(int idxLeaf = 0 ; idxLeaf < originalPosition->getSize() ; ++idxLeaf){
           leafs[(*originalPosition)[idxLeaf]].container.reduce( &transferBuffer[idxLeaf].container );
        }
    }

    // P2P
    static void p2p_cpu(void *descr[], void* cl_arg) {
        MortonContainer*const leafs = (MortonContainer*)STARPU_VECTOR_GET_PTR(descr[0]);
        const int size = STARPU_VECTOR_GET_NX(descr[0]);
        MortonContainer*const otherleafs = (MortonContainer*)STARPU_VECTOR_GET_PTR(descr[1]);
        const int othersize = STARPU_VECTOR_GET_NX(descr[1]);

        int OctreeHeight;
        starpu_codelet_unpack_args(cl_arg, &OctreeHeight);

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
                else if( othersize ){
                    const int neighPosition = findNeigh(otherleafs, othersize, potentialInteraction[idxInteraction]);
                    if( neighPosition != -1 ){
                        neighbors[ potentialPosition[idxInteraction] ] = &otherleafs[neighPosition].container;
                        ++counter;
                    }
                }
            }

            globalKernels[starpu_worker_get_id()]->P2P(leafs[idxLeaf].getCoordinate(),&leafs[idxLeaf].container,
                        &leafs[idxLeaf].container, neighbors, counter);
            if( counter ){
                memset(neighbors, 0, sizeof(ContainerClass*) * 27);
            }
        }
    }
};
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
KernelClass** FFmmAlgorithmStarpuGroup<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass>::globalKernels = 0;

#endif // FFMMALGORITHMSTARPUGROUP_HPP
