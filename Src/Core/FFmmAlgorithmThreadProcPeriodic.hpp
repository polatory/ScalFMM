// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FFMMALGORITHMTHREADPROCPPERIODIC_HPP
#define FFMMALGORITHMTHREADPROCPPERIODIC_HPP


#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"
#include "../Utils/FMemUtils.hpp"

#include "../Containers/FBoolArray.hpp"
#include "../Containers/FOctree.hpp"
#include "../Containers/FLightOctree.hpp"

#include "../Containers/FBufferWriter.hpp"
#include "../Containers/FBufferReader.hpp"

#include "../Utils/FMpi.hpp"

#include <omp.h>


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmThreadProcPeriodic
* @brief
* Please read the license
*
* This class is a threaded FMM algorithm with mpi.
* It just iterates on a tree and call the kernels with good arguments.
* It used the inspector-executor model :
* iterates on the tree and builds an array to work in parallel on this array
*
* Of course this class does not deallocate pointer given in arguements.
*
* Threaded & based on the inspector-executor model
* schedule(runtime) export OMP_NUM_THREADS=2
* export OMPI_CXX=`which g++-4.4`
* mpirun -np 2 valgrind --suppressions=/usr/share/openmpi/openmpi-valgrind.supp
* --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes
* ./Tests/testFmmAlgorithmProc ../Data/testLoaderSmall.fma.tmp
*/
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmThreadProcPeriodic : protected FAssertable {

    static const int MaxSizePerCell = 2048;

    OctreeClass* const tree;                 //< The octree to work on
    KernelClass** kernels;                   //< The kernels    

    const FMpi::FComm& comm;                 //< MPI comm

    CellClass rootCellFromProc;     //< root of tree needed by the periodicity
    const int nbLevelsAboveRoot;    //< The nb of level the user ask to go above the tree (>= -1)
    const int offsetRealTree;       //< nbLevelsAboveRoot GetFackLevel
    const int periodicDirections;

    typename OctreeClass::Iterator* iterArray;
    int numberOfLeafs;                          //< To store the size at the previous level

    const int MaxThreads;               //< the max number of thread allowed by openmp

    const int nbProcess;                //< Number of process
    const int idProcess;                //< Id of current process

    const int OctreeHeight;


    struct Interval{
        MortonIndex min;
        MortonIndex max;
    };
    Interval*const intervals;
    Interval*const workingIntervalsPerLevel;


    Interval& getWorkingInterval(const int level, const int proc){
        return workingIntervalsPerLevel[OctreeHeight * proc + level];
    }

    static int GetFackLevel(const int inLevelAboveRequiered){
        if( inLevelAboveRequiered == -1 ) return 1;
        if( inLevelAboveRequiered == 0  ) return 2;
        return inLevelAboveRequiered + 3;
    }


public:

    void setKernel(KernelClass*const inKernels){
        this->kernels = new KernelClass*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass(*inKernels);
        }
    }

    Interval& getWorkingInterval(const int level){
        return getWorkingInterval(level, idProcess);
    }

    bool hasWorkAtLevel(const int level){
        return idProcess == 0 || getWorkingInterval(level, idProcess - 1).max < getWorkingInterval(level, idProcess).max;
    }

    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThreadProcPeriodic(const FMpi::FComm& inComm, OctreeClass* const inTree,
                                    const int inUpperLevel = 2, const int inPeriodicDirections = AllDirs)
        : tree(inTree) , kernels(0), comm(inComm), nbLevelsAboveRoot(inUpperLevel), offsetRealTree(GetFackLevel(inUpperLevel)),
          periodicDirections(inPeriodicDirections), numberOfLeafs(0),
          MaxThreads(omp_get_max_threads()), nbProcess(inComm.processCount()), idProcess(inComm.processId()),
          OctreeHeight(tree->getHeight()),intervals(new Interval[inComm.processCount()]),
          workingIntervalsPerLevel(new Interval[inComm.processCount() * tree->getHeight()]) {

        fassert(tree, "tree cannot be null", __LINE__, __FILE__);
        fassert(-1 <= inUpperLevel, "inUpperLevel cannot be < -1", __LINE__, __FILE__);

        FDEBUG(FDebug::Controller << "FFmmAlgorithmThreadProcPeriodic\n");
        FDEBUG(FDebug::Controller << "Max threads = "  << MaxThreads << ", Procs = " << nbProcess << ", I am " << idProcess << ".\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmThreadProcPeriodic(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;

        delete [] intervals;
        delete [] workingIntervalsPerLevel;
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::FFunction functionTrace( __FUNCTION__, "Fmm" , __FILE__ , __LINE__ ) );

        // Count leaf
        this->numberOfLeafs = 0;
        {
            FTRACE( FTrace::FRegion regionTrace( "Preprocess" , __FUNCTION__ , __FILE__ , __LINE__) );

            Interval myLastInterval;
            {
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoBottomLeft();
                myLastInterval.min = octreeIterator.getCurrentGlobalIndex();
                do{
                    ++this->numberOfLeafs;
                } while(octreeIterator.moveRight());
                myLastInterval.max = octreeIterator.getCurrentGlobalIndex();
            }
            iterArray = new typename OctreeClass::Iterator[numberOfLeafs];
            fassert(iterArray, "iterArray bad alloc", __LINE__, __FILE__);

            // We get the min/max indexes from each procs
            FMpi::MpiAssert( MPI_Allgather( &myLastInterval, sizeof(Interval), MPI_BYTE, intervals, sizeof(Interval), MPI_BYTE, comm.getComm()),  __LINE__ );

            Interval*const myIntervals = new Interval[OctreeHeight];
            myIntervals[OctreeHeight - 1] = myLastInterval;
            for(int idxLevel = OctreeHeight - 2 ; idxLevel >= 0 ; --idxLevel){
                myIntervals[idxLevel].min = myIntervals[idxLevel+1].min >> 3;
                myIntervals[idxLevel].max = myIntervals[idxLevel+1].max >> 3;
            }
            if(idProcess != 0){
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoBottomLeft();
                octreeIterator.moveUp();

                MortonIndex currentLimit = intervals[idProcess-1].max >> 3;

                for(int idxLevel = OctreeHeight - 2 ; idxLevel >= 1 ; --idxLevel){
                    while(octreeIterator.getCurrentGlobalIndex() <= currentLimit){
                        if( !octreeIterator.moveRight() ) break;
                    }
                    myIntervals[idxLevel].min = octreeIterator.getCurrentGlobalIndex();
                    octreeIterator.moveUp();
                    currentLimit >>= 3;
                }
            }

            // We get the min/max indexes from each procs
            FMpi::MpiAssert( MPI_Allgather( myIntervals, int(sizeof(Interval)) * OctreeHeight, MPI_BYTE,
                                            workingIntervalsPerLevel, int(sizeof(Interval)) * OctreeHeight, MPI_BYTE, comm.getComm()),  __LINE__ );
            delete[] myIntervals;
        }

        // run;
        bottomPass();

        upwardPass();

        transferPass();

        downardPass();

        directPass();

        // delete array
        delete [] iterArray;
        iterArray = 0;
    }

private:

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M Bottom Pass */
    void bottomPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        int leafs = 0;
        do{
            iterArray[leafs++] = octreeIterator;
        } while(octreeIterator.moveRight());

        FDEBUG(FTic computationCounter);
        #pragma omp parallel
        {
            KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for nowait
            for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentListSrc());
            }
        }
        FDEBUG(computationCounter.tac());


        FDEBUG( FDebug::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );

    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic computationCounter);
        FDEBUG(FTic prepareCounter);
        FDEBUG(FTic waitCounter);

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // This variable is the proc responsible
        // of the shared cells
        int sendToProc = idProcess;

        // There are a maximum of 8-1 sends and 8-1 receptions
        MPI_Request requests[14];
        MPI_Status status[14];

        // Maximum data per message is:
        FBufferWriter sendBuffer;
        const int recvBufferOffset = 8 * MaxSizePerCell + 1;
        FBufferReader recvBuffer(nbProcess * recvBufferOffset);
        CellClass recvBufferCells[8];

        int firstProcThatSend = idProcess + 1;

        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 0 ; --idxLevel ){
            const int fackLevel = idxLevel + offsetRealTree;
            // No more work for me
            if(idProcess != 0
                    && getWorkingInterval((idxLevel+1), idProcess).max <= getWorkingInterval((idxLevel+1), idProcess - 1).max){
                break;
            }

            // copy cells to work with
            int numberOfCells = 0;
            // for each cells
            do{
                iterArray[numberOfCells++] = octreeIterator;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;

            // We may need to send something
            int iterRequests = 0;
            int cellsToSend = -1;

            while(iterArray[cellsToSend+1].getCurrentGlobalIndex() < getWorkingInterval(idxLevel, idProcess).min){
                ++cellsToSend;
            }

            FTRACE( FTrace::FRegion regionTrace( "Preprocess" , __FUNCTION__ , __FILE__ , __LINE__) );

            FDEBUG(prepareCounter.tic());
            if(idProcess != 0
                    && (getWorkingInterval((idxLevel+1), idProcess).min >>3) <= (getWorkingInterval((idxLevel+1), idProcess - 1).max >>3)){

                char state = 0;
                sendBuffer.write(char(0));

                const CellClass* const* const child = iterArray[cellsToSend].getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    if( child[idxChild] && getWorkingInterval((idxLevel+1), idProcess).min <= child[idxChild]->getMortonIndex() ){
                        child[idxChild]->serializeUp(sendBuffer);

                        state = char(state | (0x1 << idxChild));
                    }
                }
                sendBuffer.writeAt(0,state);

                while( sendToProc && iterArray[cellsToSend].getCurrentGlobalIndex() == getWorkingInterval(idxLevel , sendToProc - 1).max){
                    --sendToProc;
                }

                MPI_Isend(sendBuffer.data(), sendBuffer.getSize(), MPI_BYTE, sendToProc, FMpi::TagFmmM2M, comm.getComm(), &requests[iterRequests++]);
            }

            // We may need to receive something
            bool hasToReceive = false;
            int endProcThatSend = firstProcThatSend;

            if(idProcess != nbProcess - 1){
                while(firstProcThatSend < nbProcess
                      && getWorkingInterval((idxLevel+1), firstProcThatSend).max < getWorkingInterval((idxLevel+1), idProcess).max){
                    ++firstProcThatSend;
                }

                if(firstProcThatSend < nbProcess &&
                        (getWorkingInterval((idxLevel+1), firstProcThatSend).min >>3) <= (getWorkingInterval((idxLevel+1) , idProcess).max>>3) ){

                    endProcThatSend = firstProcThatSend;

                    while( endProcThatSend < nbProcess &&
                           (getWorkingInterval((idxLevel+1) ,endProcThatSend).min >>3) <= (getWorkingInterval((idxLevel+1) , idProcess).max>>3)){
                        ++endProcThatSend;
                    }


                    if(firstProcThatSend != endProcThatSend){
                        hasToReceive = true;

                        for(int idxProc = firstProcThatSend ; idxProc < endProcThatSend ; ++idxProc ){
                            MPI_Irecv(&recvBuffer.data()[idxProc * recvBufferOffset], recvBufferOffset, MPI_BYTE,
                                      idxProc, FMpi::TagFmmM2M, comm.getComm(), &requests[iterRequests++]);
                        }
                    }
                }
            }
            FDEBUG(prepareCounter.tac());
            FTRACE( regionTrace.end() );

            // Compute
            const int endIndex = (hasToReceive?numberOfCells-1:numberOfCells);
            FDEBUG(computationCounter.tic());
            #pragma omp parallel
            {
                KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
                #pragma omp for nowait
                for( int idxCell = cellsToSend + 1 ; idxCell < endIndex ; ++idxCell){
                    myThreadkernels.M2M( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), fackLevel);
                }
            }
            FDEBUG(computationCounter.tac());

            // Are we sending or waiting anything?
            if(iterRequests){
                FDEBUG(waitCounter.tic());
                MPI_Waitall( iterRequests, requests, status);
                FDEBUG(waitCounter.tac());

                // we were receiving data
                if( hasToReceive ){
                    CellClass* currentChild[8];
                    memcpy(currentChild, iterArray[numberOfCells - 1].getCurrentChild(), 8 * sizeof(CellClass*));

                    // retreive data and merge my child and the child from others
                    for(int idxProc = firstProcThatSend ; idxProc < endProcThatSend ; ++idxProc){
                        recvBuffer.seek(idxProc * recvBufferOffset);
                        int state = int(recvBuffer.getValue<char>());

                        int position = 0;
                        while( state && position < 8){
                            while(!(state & 0x1)){
                                state >>= 1;
                                ++position;
                            }

                            fassert(!currentChild[position], "Already has a cell here", __LINE__, __FILE__);

                            recvBufferCells[position].deserializeUp(recvBuffer);
                            currentChild[position] = (CellClass*) &recvBufferCells[position];

                            state >>= 1;
                            ++position;
                        }
                    }

                    // Finally compute
                    FDEBUG(computationCounter.tic());
                    (*kernels[0]).M2M( iterArray[numberOfCells - 1].getCurrentCell() , currentChild, fackLevel);
                    FDEBUG(computationCounter.tac());


                    firstProcThatSend = endProcThatSend - 1;
                }
            }

            sendBuffer.reset();
            recvBuffer.seek(0);
        }


        FDEBUG( FDebug::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Prepare : " << prepareCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Wait : " << waitCounter.cumulated() << " s\n" );

        //////////////////////////////////////////////////////////////////
        //Periodicity
        //////////////////////////////////////////////////////////////////

        octreeIterator = typename OctreeClass::Iterator(tree);

        if( idProcess == 0){
            int iterRequests = 0;

            CellClass* currentChild[8];
            memcpy(currentChild, octreeIterator.getCurrentBox(), 8 * sizeof(CellClass*));

            for(int idxProc = 1 ; idxProc < nbProcess ; ++idxProc ){
                if( getWorkingInterval(1, idxProc - 1).max < getWorkingInterval(1, idxProc).max ){
                    MPI_Irecv(&recvBuffer.data()[idxProc * recvBufferOffset], recvBufferOffset, MPI_BYTE, idxProc,
                              FMpi::TagFmmM2M, comm.getComm(), &requests[iterRequests++]);
                }
            }

            MPI_Waitall( iterRequests, requests, MPI_STATUSES_IGNORE);

            // retreive data and merge my child and the child from others
            for(int idxProc = 1 ; idxProc < nbProcess ; ++idxProc){
                if( getWorkingInterval(1, idxProc - 1).max < getWorkingInterval(1, idxProc).max ){
                    recvBuffer.seek(idxProc * recvBufferOffset);
                    int state = int(recvBuffer.getValue<char>());

                    int position = 0;
                    while( state && position < 8){
                        while(!(state & 0x1)){
                            state >>= 1;
                            ++position;
                        }
                        fassert(!currentChild[position], "Already has a cell here", __LINE__, __FILE__);

                        recvBufferCells[position].deserializeUp(recvBuffer);

                        currentChild[position] = (CellClass*) &recvBufferCells[position];

                        state >>= 1;
                        ++position;
                    }
                }
            }

            (*kernels[0]).M2M( &rootCellFromProc , currentChild, offsetRealTree);

            processPeriodicLevels();
        }
        else {
            if( hasWorkAtLevel(1) ){
                const int firstChild = getWorkingInterval(1, idProcess).min & 7;
                const int lastChild = getWorkingInterval(1, idProcess).max & 7;

                CellClass** child = octreeIterator.getCurrentBox();

                char state = 0;
                sendBuffer.write(state);

                for(int idxChild = firstChild ; idxChild <= lastChild ; ++idxChild){
                    if( child[idxChild] ){
                        child[idxChild]->serializeUp(sendBuffer);
                        state = char( state | (0x1 << idxChild));

                    }
                }
                sendBuffer.writeAt(0,state);

                MPI_Send(sendBuffer.data(), sendBuffer.getSize(), MPI_BYTE, 0, FMpi::TagFmmM2M, comm.getComm());
            }
        }

    }

    /////////////////////////////////////////////////////////////////////////////
    // Downard
    /////////////////////////////////////////////////////////////////////////////

    /** M2L  */
    void transferPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic computationCounter);
        FDEBUG(FTic sendCounter);
        FDEBUG(FTic receiveCounter);
        FDEBUG(FTic prepareCounter);
        FDEBUG(FTic gatherCounter);

        //////////////////////////////////////////////////////////////////
        // First know what to send to who
        //////////////////////////////////////////////////////////////////

        // pointer to send
        FVector<typename OctreeClass::Iterator> toSend[nbProcess * OctreeHeight];
        // index
        int*const indexToSend = new int[nbProcess * OctreeHeight];
        memset(indexToSend, 0, sizeof(int) * nbProcess * OctreeHeight);
        // To know which one has need someone
        FBoolArray** const leafsNeedOther = new FBoolArray*[OctreeHeight];
        memset(leafsNeedOther, 0, sizeof(FBoolArray*) * OctreeHeight);

        {
            FTRACE( FTrace::FRegion regionTrace( "Preprocess" , __FUNCTION__ , __FILE__ , __LINE__) );
            FDEBUG(prepareCounter.tic());

            // To know if a leaf has been already sent to a proc
            bool*const alreadySent = new bool[nbProcess];
            memset(alreadySent, 0, sizeof(bool) * nbProcess);

            typename OctreeClass::Iterator octreeIterator(tree);
            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
            // for each levels
            for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
                if(idProcess != 0
                        && getWorkingInterval(idxLevel, idProcess).max <= getWorkingInterval(idxLevel, idProcess - 1).max){
                    avoidGotoLeftIterator.moveDown();
                    octreeIterator = avoidGotoLeftIterator;

                    continue;
                }

                int numberOfCells = 0;

                while(octreeIterator.getCurrentGlobalIndex() <  getWorkingInterval(idxLevel , idProcess).min){
                    octreeIterator.moveRight();
                }

                // for each cells
                do{
                    iterArray[numberOfCells] = octreeIterator;
                    ++numberOfCells;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                leafsNeedOther[idxLevel] = new FBoolArray(numberOfCells);


                // Which cell potentialy needs other data and in the same time
                // are potentialy needed by other
                int neighborsPosition[189];
                MortonIndex neighborsIndexes[189];
                for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                    // Find the M2L neigbors of a cell
                    const int counter = getPeriodicInteractionNeighbors(iterArray[idxCell].getCurrentGlobalCoordinate(),idxLevel,
                                                                        neighborsIndexes, neighborsPosition, periodicDirections);

                    memset(alreadySent, false, sizeof(bool) * nbProcess);
                    bool needOther = false;
                    // Test each negibors to know which one do not belong to us
                    for(int idxNeigh = 0 ; idxNeigh < counter ; ++idxNeigh){
                        if(neighborsIndexes[idxNeigh] < getWorkingInterval(idxLevel , idProcess).min
                                || getWorkingInterval(idxLevel , idProcess).max < neighborsIndexes[idxNeigh]){
                            int procToReceive = idProcess;
                            while( 0 != procToReceive && neighborsIndexes[idxNeigh] < getWorkingInterval(idxLevel , procToReceive).min ){
                                --procToReceive;
                            }
                            while( procToReceive != nbProcess -1 && getWorkingInterval(idxLevel , procToReceive).max < neighborsIndexes[idxNeigh]){
                                ++procToReceive;
                            }

                            // Maybe already sent to that proc?
                            if( !alreadySent[procToReceive]
                                    && getWorkingInterval(idxLevel , procToReceive).min <= neighborsIndexes[idxNeigh]
                                    && neighborsIndexes[idxNeigh] <= getWorkingInterval(idxLevel , procToReceive).max){

                                alreadySent[procToReceive] = true;

                                needOther = true;

                                toSend[idxLevel * nbProcess + procToReceive].push(iterArray[idxCell]);
                                ++indexToSend[idxLevel * nbProcess + procToReceive];
                            }
                        }
                    }
                    if(needOther){
                        leafsNeedOther[idxLevel]->set(idxCell,true);
                    }

                }

            }
            FDEBUG(prepareCounter.tac());

            delete[] alreadySent;
        }

        //////////////////////////////////////////////////////////////////
        // Gather this information
        //////////////////////////////////////////////////////////////////

        FDEBUG(gatherCounter.tic());
        // All process say to each others
        // what the will send to who
        int*const globalReceiveMap = new int[nbProcess * nbProcess * OctreeHeight];
        memset(globalReceiveMap, 0, sizeof(int) * nbProcess * nbProcess * OctreeHeight);
        FMpi::MpiAssert( MPI_Allgather( indexToSend, nbProcess * OctreeHeight, MPI_INT, globalReceiveMap, nbProcess * OctreeHeight, MPI_INT, comm.getComm()),  __LINE__ );
        FDEBUG(gatherCounter.tac());

        //////////////////////////////////////////////////////////////////
        // Send and receive for real
        //////////////////////////////////////////////////////////////////

        FDEBUG(sendCounter.tic());
        // Then they can send and receive (because they know what they will receive)
        // To send in asynchrone way
        MPI_Request*const requests = new MPI_Request[2 * nbProcess * OctreeHeight];
        MPI_Status*const status = new MPI_Status[2 * nbProcess * OctreeHeight];
        int iterRequest = 0;

        const int SizeOfCellToSend = sizeof(MortonIndex) + sizeof(int) + MaxSizePerCell;

        FBufferWriter**const sendBuffer = new FBufferWriter*[nbProcess * OctreeHeight];
        memset(sendBuffer, 0, sizeof(FBufferWriter*) * nbProcess * OctreeHeight);

        FBufferReader**const recvBuffer = new FBufferReader*[nbProcess * OctreeHeight];
        memset(recvBuffer, 0, sizeof(FBufferReader*) * nbProcess * OctreeHeight);


        for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
            for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                const int toSendAtProcAtLevel = indexToSend[idxLevel * nbProcess + idxProc];
                if(toSendAtProcAtLevel != 0){
                    sendBuffer[idxLevel * nbProcess + idxProc] = new FBufferWriter(toSendAtProcAtLevel * SizeOfCellToSend);

                    for(int idxLeaf = 0 ; idxLeaf < toSendAtProcAtLevel; ++idxLeaf){
                        const MortonIndex cellIndex = toSend[idxLevel * nbProcess + idxProc][idxLeaf].getCurrentGlobalIndex();
                        sendBuffer[idxLevel * nbProcess + idxProc]->write(cellIndex);
                        toSend[idxLevel * nbProcess + idxProc][idxLeaf].getCurrentCell()->serializeUp(*sendBuffer[idxLevel * nbProcess + idxProc]);
                    }

                    FMpi::MpiAssert( MPI_Isend( sendBuffer[idxLevel * nbProcess + idxProc]->data(), sendBuffer[idxLevel * nbProcess + idxProc]->getSize()
                                                , MPI_BYTE , idxProc, FMpi::TagLast + idxLevel, comm.getComm(), &requests[iterRequest++]) , __LINE__ );

                }

                const int toReceiveFromProcAtLevel = globalReceiveMap[(idxProc * nbProcess * OctreeHeight) + idxLevel * nbProcess + idProcess];
                if(toReceiveFromProcAtLevel){
                    recvBuffer[idxLevel * nbProcess + idxProc] = new FBufferReader(toReceiveFromProcAtLevel * SizeOfCellToSend);

                    FMpi::MpiAssert( MPI_Irecv(recvBuffer[idxLevel * nbProcess + idxProc]->data(), recvBuffer[idxLevel * nbProcess + idxProc]->getSize(), MPI_BYTE,
                                               idxProc, FMpi::TagLast + idxLevel, comm.getComm(), &requests[iterRequest++]) , __LINE__ );
                }
            }
        }
        FDEBUG(sendCounter.tac());

        //////////////////////////////////////////////////////////////////
        // Do M2L
        //////////////////////////////////////////////////////////////////

        {
            FTRACE( FTrace::FRegion regionTrace("Compute", __FUNCTION__ , __FILE__ , __LINE__) );
            typename OctreeClass::Iterator octreeIterator(tree);
            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
            // Now we can compute all the data
            // for each levels
            for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){                
                const int fackLevel = idxLevel + offsetRealTree;
                if(idProcess != 0
                        && getWorkingInterval(idxLevel, idProcess).max <= getWorkingInterval(idxLevel, idProcess - 1).max){

                    avoidGotoLeftIterator.moveDown();
                    octreeIterator = avoidGotoLeftIterator;

                    continue;
                }

                int numberOfCells = 0;
                while(octreeIterator.getCurrentGlobalIndex() <  getWorkingInterval(idxLevel , idProcess).min){
                    octreeIterator.moveRight();
                }
                // for each cells
                do{
                    iterArray[numberOfCells] = octreeIterator;
                    ++numberOfCells;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                FDEBUG(computationCounter.tic());
                #pragma omp parallel
                {
                    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                    const CellClass* neighbors[343];

                    #pragma omp for  schedule(dynamic) nowait
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        const int counter = tree->getPeriodicInteractionNeighbors(neighbors, iterArray[idxCell].getCurrentGlobalCoordinate(),idxLevel, periodicDirections);
                        if(counter) myThreadkernels->M2L( iterArray[idxCell].getCurrentCell() , neighbors, counter, fackLevel);
                    }

                    myThreadkernels->finishedLevelM2L(fackLevel);
                }
                FDEBUG(computationCounter.tac());
            }
        }

        //////////////////////////////////////////////////////////////////
        // Wait received data and compute
        //////////////////////////////////////////////////////////////////

        // Wait to receive every things (and send every things)
        MPI_Waitall(iterRequest, requests, status);

        {
            FTRACE( FTrace::FRegion regionTrace("Compute Received data", __FUNCTION__ , __FILE__ , __LINE__) );
            FDEBUG(receiveCounter.tic());
            typename OctreeClass::Iterator octreeIterator(tree);
            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
            // compute the second time
            // for each levels
            for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
                const int fackLevel = idxLevel + offsetRealTree;
                if(idProcess != 0
                        && getWorkingInterval(idxLevel, idProcess).max <= getWorkingInterval(idxLevel, idProcess - 1).max){

                    avoidGotoLeftIterator.moveDown();
                    octreeIterator = avoidGotoLeftIterator;

                    continue;
                }

                // put the received data into a temporary tree
                FLightOctree<CellClass> tempTree;
                for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                    const int toReceiveFromProcAtLevel = globalReceiveMap[(idxProc * nbProcess * OctreeHeight) + idxLevel * nbProcess + idProcess];

                    for(int idxCell = 0 ; idxCell < toReceiveFromProcAtLevel ; ++idxCell){
                        const MortonIndex cellIndex = recvBuffer[idxLevel * nbProcess + idxProc]->FBufferReader::getValue<MortonIndex>();

                        CellClass* const newCell = new CellClass;
                        newCell->setMortonIndex(cellIndex);
                        newCell->deserializeUp(*recvBuffer[idxLevel * nbProcess + idxProc]);

                        tempTree.insertCell(cellIndex, idxLevel, newCell);
                    }
                }


                // take cells from our octree only if they are
                // linked to received data
                int numberOfCells = 0;
                int realCellId = 0;

                while(octreeIterator.getCurrentGlobalIndex() <  getWorkingInterval(idxLevel , idProcess).min){
                    octreeIterator.moveRight();
                }
                // for each cells
                do{
                    // copy cells that need data from others
                    if(leafsNeedOther[idxLevel]->get(realCellId++)){
                        iterArray[numberOfCells++] = octreeIterator;
                    }
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                delete leafsNeedOther[idxLevel];
                leafsNeedOther[idxLevel] = 0;

                // Compute this cells
                FDEBUG(computationCounter.tic());
                #pragma omp parallel
                {
                    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                    MortonIndex neighborsIndex[189];
                    int neighborsPosition[189];
                    const CellClass* neighbors[343];

                    #pragma omp for schedule(dynamic) nowait
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        // compute indexes
                        memset(neighbors, 0, 343 * sizeof(CellClass*));
                        const int counterNeighbors = getPeriodicInteractionNeighbors(iterArray[idxCell].getCurrentGlobalCoordinate(), idxLevel,
                                                                         neighborsIndex, neighborsPosition, periodicDirections);
                        int counter = 0;
                        // does we receive this index from someone?
                        for(int idxNeig = 0 ;idxNeig < counterNeighbors ; ++idxNeig){

                            if(neighborsIndex[idxNeig] < getWorkingInterval(idxLevel , idProcess).min
                                    || getWorkingInterval(idxLevel , idProcess).max < neighborsIndex[idxNeig]){

                                CellClass*const otherCell = tempTree.getCell(neighborsIndex[idxNeig], idxLevel);

                                if(otherCell){
                                    //otherCell->setMortonIndex(neighborsIndex[idxNeig]);
                                    neighbors[ neighborsPosition[idxNeig] ] = otherCell;
                                    ++counter;
                                }
                            }
                        }

                        // need to compute
                        if(counter){
                            myThreadkernels->M2L( iterArray[idxCell].getCurrentCell() , neighbors, counter, fackLevel);
                        }
                    }

                    myThreadkernels->finishedLevelM2L(fackLevel);
                }
                FDEBUG(computationCounter.tac());
            }
            FDEBUG(receiveCounter.tac());
        }

        for(int idxComm = 0 ; idxComm < nbProcess * OctreeHeight; ++idxComm){
            delete sendBuffer[idxComm];
            delete recvBuffer[idxComm];
        }
        for(int idxComm = 0 ; idxComm < OctreeHeight; ++idxComm){
            delete leafsNeedOther[idxComm];
        }
        delete[] sendBuffer;
        delete[] recvBuffer;
        delete[] indexToSend;
        delete[] leafsNeedOther;
        delete[] globalReceiveMap;
        delete[] requests;
        delete[] status;

        FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Gather : " << gatherCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Prepare : " << prepareCounter.cumulated() << " s\n" );
    }

    //////////////////////////////////////////////////////////////////
    // ---------------- L2L ---------------
    //////////////////////////////////////////////////////////////////

    void downardPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic computationCounter);
        FDEBUG(FTic prepareCounter);
        FDEBUG(FTic waitCounter);

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(tree);
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        MPI_Request*const requests = new MPI_Request[nbProcess];
        MPI_Status*const status = new MPI_Status[nbProcess];

        const int heightMinusOne = OctreeHeight - 1;

        FBufferWriter sendBuffer;
        FBufferReader recvBuffer(MaxSizePerCell);

        // Periodic
        if( idProcess == 0){
            rootCellFromProc.serializeDown(sendBuffer);
            int sizeOfSerialization = sendBuffer.getSize();
            FMpi::MpiAssert( MPI_Bcast( &sizeOfSerialization, 1, MPI_INT, 0, comm.getComm() ), __LINE__ );
            FMpi::MpiAssert( MPI_Bcast( sendBuffer.data(), sendBuffer.getSize(), MPI_BYTE, 0, comm.getComm() ), __LINE__ );
            sendBuffer.reset();
        }
        else{
            int sizeOfSerialization = -1;
            FMpi::MpiAssert( MPI_Bcast( &sizeOfSerialization, 1, MPI_INT, 0, comm.getComm() ), __LINE__ );
            FMpi::MpiAssert( MPI_Bcast( recvBuffer.data(), sizeOfSerialization, MPI_BYTE, 0, comm.getComm() ), __LINE__ );
            rootCellFromProc.deserializeDown(recvBuffer);
            recvBuffer.seek(0);
        }

        kernels[0]->L2L(&rootCellFromProc, octreeIterator.getCurrentBox(), offsetRealTree);

        // for each levels exepted leaf level
        for(int idxLevel = 1 ; idxLevel < heightMinusOne ; ++idxLevel ){
            const int fackLevel = idxLevel + offsetRealTree;
            if(idProcess != 0
                    && getWorkingInterval((idxLevel+1) , idProcess).max <= getWorkingInterval((idxLevel+1) , idProcess - 1).max){

                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                continue;
            }

            // copy cells to work with
            int numberOfCells = 0;
            // for each cells
            do{
                iterArray[numberOfCells++] = octreeIterator;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;

            int firstCellWork = -1;
            while(iterArray[firstCellWork+1].getCurrentGlobalIndex() < getWorkingInterval(idxLevel , idProcess).min){
                ++firstCellWork;
            }

            bool needToRecv = false;
            int iterRequests = 0;

            FDEBUG(prepareCounter.tic());

            // do we need to receive one or zeros cell
            if(idProcess != 0
                    && (getWorkingInterval((idxLevel + 1) , idProcess).min >> 3 ) <= (getWorkingInterval((idxLevel+1) , idProcess - 1).max >> 3 ) ){
                needToRecv = true;

                MPI_Irecv( recvBuffer.data(), recvBuffer.getSize(), MPI_BYTE, MPI_ANY_SOURCE,
                           FMpi::TagFmmL2L, comm.getComm(), &requests[iterRequests++]);
            }


            if(idProcess != nbProcess - 1){
                int firstProcThatRecv = idProcess + 1;
                while( firstProcThatRecv < nbProcess &&
                       getWorkingInterval((idxLevel + 1) , firstProcThatRecv).max <= getWorkingInterval((idxLevel+1) , idProcess).max){
                    ++firstProcThatRecv;
                }

                int endProcThatRecv = firstProcThatRecv;
                while( endProcThatRecv < nbProcess &&
                       (getWorkingInterval((idxLevel + 1) , endProcThatRecv).min >> 3) <= (getWorkingInterval((idxLevel+1) , idProcess).max >> 3) ){
                    ++endProcThatRecv;
                }

                if(firstProcThatRecv != endProcThatRecv){
                    iterArray[numberOfCells - 1].getCurrentCell()->serializeDown(sendBuffer);

                    for(int idxProc = firstProcThatRecv ; idxProc < endProcThatRecv ; ++idxProc ){
                        MPI_Isend(sendBuffer.data(), sendBuffer.getSize(), MPI_BYTE, idxProc,
                                  FMpi::TagFmmL2L, comm.getComm(), &requests[iterRequests++]);
                    }

                }
            }
            FDEBUG(prepareCounter.tac());

            FDEBUG(computationCounter.tic());
            #pragma omp parallel
            {
                KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
                #pragma omp for nowait
                for(int idxCell = firstCellWork + 1 ; idxCell < numberOfCells ; ++idxCell){
                    myThreadkernels.L2L( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), fackLevel);
                }
            }
            FDEBUG(computationCounter.tac());

            // are we sending or receiving?
            if(iterRequests){
                // process
                FDEBUG(waitCounter.tic());
                MPI_Waitall( iterRequests, requests, status);
                FDEBUG(waitCounter.tac());

                if(needToRecv){
                    // Need to compute
                    FDEBUG(computationCounter.tic());
                    iterArray[firstCellWork].getCurrentCell()->deserializeDown(recvBuffer);

                    kernels[0]->L2L( iterArray[firstCellWork].getCurrentCell() , iterArray[firstCellWork].getCurrentChild(), fackLevel);
                    FDEBUG(computationCounter.tac());
                }
            }

            sendBuffer.reset();
            recvBuffer.seek(0);
        }

        delete[] requests;
        delete[] status;

        FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Prepare : " << prepareCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Wait : " << waitCounter.cumulated() << " s\n" );
    }


    /////////////////////////////////////////////////////////////////////////////
    // Direct
    /////////////////////////////////////////////////////////////////////////////
    struct LeafData{
        MortonIndex index;
        CellClass* cell;
        ContainerClass* targets;
        ContainerClass* sources;
    };
    /** P2P */
    void directPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG( FTic counterTime);
        FDEBUG( FTic prepareCounter);
        FDEBUG( FTic gatherCounter);
        FDEBUG( FTic waitCounter);

        ///////////////////////////////////////////////////
        // Prepare data to send receive
        ///////////////////////////////////////////////////
        FDEBUG(prepareCounter.tic());

        // To send in asynchrone way
        MPI_Request requests[2 * nbProcess];
        MPI_Status status[2 * nbProcess];
        int iterRequest = 0;
        int nbMessagesToRecv = 0;

        FBufferWriter**const sendBuffer = new FBufferWriter*[nbProcess];
        memset(sendBuffer, 0, sizeof(FBufferWriter*) * nbProcess);

        FBufferReader**const recvBuffer = new FBufferReader*[nbProcess];
        memset(recvBuffer, 0, sizeof(FBufferReader*) * nbProcess);

        int*const globalReceiveMap = new int[nbProcess * nbProcess];
        memset(globalReceiveMap, 0, sizeof(int) * nbProcess * nbProcess);

        FBoolArray leafsNeedOther(this->numberOfLeafs);
        int countNeedOther = 0;

        {
            FTRACE( FTrace::FRegion regionTrace( "Preprocess" , __FUNCTION__ , __FILE__ , __LINE__) );
            // Copy leafs
            {
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoBottomLeft();
                int idxLeaf = 0;
                do{
                    this->iterArray[idxLeaf++] = octreeIterator;
                } while(octreeIterator.moveRight());
            }

            // Box limite
            const int limite = 1 << (this->OctreeHeight - 1);
            // pointer to send
            FVector<typename OctreeClass::Iterator>*const toSend = new FVector<typename OctreeClass::Iterator>[nbProcess];

            // index
            int partsToSend[nbProcess];
            memset(partsToSend, 0, sizeof(int) * nbProcess);

            // To know if a leaf has been already sent to a proc
            int alreadySent[nbProcess];

            MortonIndex indexesNeighbors[26];
            int uselessIndexArray[26];

            for(int idxLeaf = 0 ; idxLeaf < this->numberOfLeafs ; ++idxLeaf){
                FTreeCoordinate center;
                center.setPositionFromMorton(iterArray[idxLeaf].getCurrentGlobalIndex(), OctreeHeight - 1);

                memset(alreadySent, 0, sizeof(int) * nbProcess);
                bool needOther = false;

                const int neighCount = getNeighborsIndexesPeriodic(iterArray[idxLeaf].getCurrentGlobalCoordinate(),
                                                                   limite, indexesNeighbors, uselessIndexArray, AllDirs);

                for(int idxNeigh = 0 ; idxNeigh < neighCount ; ++idxNeigh){
                    if(indexesNeighbors[idxNeigh] < intervals[idProcess].min || intervals[idProcess].max < indexesNeighbors[idxNeigh]){
                        needOther = true;

                        // find the proc that need this information
                        int procToReceive = idProcess;
                        while( procToReceive != 0 && indexesNeighbors[idxNeigh] < intervals[procToReceive].min){
                            --procToReceive;
                        }

                        while( procToReceive != nbProcess - 1 && intervals[procToReceive].max < indexesNeighbors[idxNeigh]){
                            ++procToReceive;
                        }

                        if( !alreadySent[procToReceive] && intervals[procToReceive].min <= indexesNeighbors[idxNeigh] && indexesNeighbors[idxNeigh] <= intervals[procToReceive].max){

                            alreadySent[procToReceive] = 1;
                            toSend[procToReceive].push( iterArray[idxLeaf] );
                            partsToSend[procToReceive] += iterArray[idxLeaf].getCurrentListSrc()->getSavedSize();
                            partsToSend[procToReceive] += int(sizeof(MortonIndex));
                        }
                    }
                }

                if(needOther){
                    leafsNeedOther.set(idxLeaf,true);
                    ++countNeedOther;
                }
            }

            for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                if(partsToSend[idxProc]){
                    partsToSend[idxProc] += int(sizeof(int));
                }
            }

            FDEBUG(gatherCounter.tic());
            FMpi::MpiAssert( MPI_Allgather( partsToSend, nbProcess, MPI_INT, globalReceiveMap, nbProcess, MPI_INT, comm.getComm()),  __LINE__ );
            FDEBUG(gatherCounter.tac());


            // Prepare receive
            for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                if(globalReceiveMap[idxProc * nbProcess + idProcess]){
                    recvBuffer[idxProc] = new FBufferReader(globalReceiveMap[idxProc * nbProcess + idProcess]);
                    FMpi::MpiAssert( MPI_Irecv(recvBuffer[idxProc]->data(), recvBuffer[idxProc]->getSize(), MPI_BYTE,
                                               idxProc, FMpi::TagFmmP2P, comm.getComm(), &requests[iterRequest++]) , __LINE__ );
                }
            }

            nbMessagesToRecv = iterRequest;
            // Prepare send
            for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
                if(toSend[idxProc].getSize() != 0){
                    sendBuffer[idxProc] = new FBufferWriter(partsToSend[idxProc]);

                    (*sendBuffer[idxProc]) << toSend[idxProc].getSize();

                    for(int idxLeaf = 0 ; idxLeaf < toSend[idxProc].getSize() ; ++idxLeaf){
                        (*sendBuffer[idxProc]) << toSend[idxProc][idxLeaf].getCurrentGlobalIndex();
                        toSend[idxProc][idxLeaf].getCurrentListSrc()->save(*sendBuffer[idxProc]);
                    }

                    FMpi::MpiAssert( MPI_Isend( sendBuffer[idxProc]->data(), sendBuffer[idxProc]->getSize() , MPI_BYTE ,
                                                idxProc, FMpi::TagFmmP2P, comm.getComm(), &requests[iterRequest++]) , __LINE__ );

                }
            }

            delete[] toSend;
        }
        FDEBUG(prepareCounter.tac());

        ///////////////////////////////////////////////////
        // Prepare data for thread P2P
        ///////////////////////////////////////////////////

        // init
        const int LeafIndex = OctreeHeight - 1;
        const int SizeShape = 3*3*3;

        int shapeLeaf[SizeShape];
        memset(shapeLeaf,0,SizeShape*sizeof(int));

        LeafData* const leafsDataArray = new LeafData[this->numberOfLeafs];

        FVector<LeafData> leafsNeedOtherData(countNeedOther);

        // split data
        {
            FTRACE( FTrace::FRegion regionTrace( "Split" , __FUNCTION__ , __FILE__ , __LINE__) );

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();

            // to store which shape for each leaf
            typename OctreeClass::Iterator* const myLeafs = new typename OctreeClass::Iterator[this->numberOfLeafs];
            int*const shapeType = new int[this->numberOfLeafs];

            for(int idxLeaf = 0 ; idxLeaf < this->numberOfLeafs ; ++idxLeaf){
                myLeafs[idxLeaf] = octreeIterator;

                const FTreeCoordinate& coord = octreeIterator.getCurrentCell()->getCoordinate();
                const int shape = (coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3);
                shapeType[idxLeaf] = shape;

                ++shapeLeaf[shape];

                octreeIterator.moveRight();
            }

            int startPosAtShape[SizeShape];
            startPosAtShape[0] = 0;
            for(int idxShape = 1 ; idxShape < SizeShape ; ++idxShape){
                startPosAtShape[idxShape] = startPosAtShape[idxShape-1] + shapeLeaf[idxShape-1];
            }

            int idxInArray = 0;
            for(int idxLeaf = 0 ; idxLeaf < this->numberOfLeafs ; ++idxLeaf, ++idxInArray){
                const int shapePosition = shapeType[idxInArray];

                leafsDataArray[startPosAtShape[shapePosition]].index = myLeafs[idxInArray].getCurrentGlobalIndex();
                leafsDataArray[startPosAtShape[shapePosition]].cell = myLeafs[idxInArray].getCurrentCell();
                leafsDataArray[startPosAtShape[shapePosition]].targets = myLeafs[idxInArray].getCurrentListTargets();
                leafsDataArray[startPosAtShape[shapePosition]].sources = myLeafs[idxInArray].getCurrentListSrc();
                if( leafsNeedOther.get(idxLeaf) ) leafsNeedOtherData.push(leafsDataArray[startPosAtShape[shapePosition]]);

                ++startPosAtShape[shapePosition];
            }

            delete[] shapeType;
            delete[] myLeafs;
        }


        //////////////////////////////////////////////////////////
        // Computation P2P that DO NOT need others data
        //////////////////////////////////////////////////////////
        FTRACE( FTrace::FRegion regionP2PTrace("Compute P2P", __FUNCTION__ , __FILE__ , __LINE__) );

        FDEBUG(FTic computationCounter);

        #pragma omp parallel
        {
            KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
            // There is a maximum of 26 neighbors
            ContainerClass* neighbors[27];
            FTreeCoordinate offsets[27];
            const FReal boxWidth = tree->getBoxWidth();
            bool hasPeriodicLeaves;
            int previous = 0;

            for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
                const int endAtThisShape = shapeLeaf[idxShape] + previous;

                #pragma omp for
                for(int idxLeafs = previous ; idxLeafs < endAtThisShape ; ++idxLeafs){
                    LeafData& currentIter = leafsDataArray[idxLeafs];
                    myThreadkernels.L2P(currentIter.cell, currentIter.targets);

                    // need the current particles and neighbors particles
                    const int counter = tree->getPeriodicLeafsNeighbors(neighbors, offsets, &hasPeriodicLeaves,
                                           currentIter.cell->getCoordinate(), LeafIndex, periodicDirections);
                    int periodicNeighborsCounter = 0;

                    if(hasPeriodicLeaves){
                        ContainerClass* periodicNeighbors[27];
                        memset(periodicNeighbors, 0, 27 * sizeof(ContainerClass*));

                        for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
                            if( neighbors[idxNeig] && !offsets[idxNeig].equals(0,0,0) ){
                                // Put periodic neighbors into other array
                                periodicNeighbors[idxNeig] = neighbors[idxNeig];
                                neighbors[idxNeig] = 0;
                                ++periodicNeighborsCounter;

                                FReal*const positionsX = periodicNeighbors[idxNeig]->getWPositions()[0];
                                FReal*const positionsY = periodicNeighbors[idxNeig]->getWPositions()[1];
                                FReal*const positionsZ = periodicNeighbors[idxNeig]->getWPositions()[2];

                                for(int idxPart = 0; idxPart < periodicNeighbors[idxNeig]->getNbParticles() ; ++idxPart){
                                    positionsX[idxPart] += boxWidth * FReal(offsets[idxNeig].getX());
                                    positionsY[idxPart] += boxWidth * FReal(offsets[idxNeig].getY());
                                    positionsZ[idxPart] += boxWidth * FReal(offsets[idxNeig].getZ());
                                }
                            }
                        }

                        myThreadkernels.P2PRemote(currentIter.cell->getCoordinate(),currentIter.targets,
                                     currentIter.sources, periodicNeighbors, periodicNeighborsCounter);

                        for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
                            if( periodicNeighbors[idxNeig] ){                                
                                FReal*const positionsX = periodicNeighbors[idxNeig]->getWPositions()[0];
                                FReal*const positionsY = periodicNeighbors[idxNeig]->getWPositions()[1];
                                FReal*const positionsZ = periodicNeighbors[idxNeig]->getWPositions()[2];

                                for(int idxPart = 0; idxPart < periodicNeighbors[idxNeig]->getNbParticles() ; ++idxPart){
                                    positionsX[idxPart] -= boxWidth * FReal(offsets[idxNeig].getX());
                                    positionsY[idxPart] -= boxWidth * FReal(offsets[idxNeig].getY());
                                    positionsZ[idxPart] -= boxWidth * FReal(offsets[idxNeig].getZ());
                                }
                            }
                        }
                    }

                    myThreadkernels.P2P( currentIter.cell->getCoordinate(), currentIter.targets,
                                         currentIter.sources, neighbors, counter - periodicNeighborsCounter);
                }

                previous = endAtThisShape;
            }
        }
        FDEBUG(computationCounter.tac());
        FTRACE( regionP2PTrace.end() );

        //////////////////////////////////////////////////////////
        // Wait send receive
        //////////////////////////////////////////////////////////

        FDEBUG(FTic computation2Counter);

        // Create an octree with leaves from others
        OctreeClass otherP2Ptree( tree->getHeight(), tree->getSubHeight(), tree->getBoxWidth(), tree->getBoxCenter() );
        int complete = 0;
        int*const indexMessage = new int[nbProcess * 2];
        while( complete != iterRequest){
            memset(indexMessage, 0, sizeof(int) * nbProcess * 2);
            int countMessages = 0;
            // Wait data
            FDEBUG(waitCounter.tic());
            MPI_Waitsome(iterRequest, requests, &countMessages, indexMessage, status);
            FDEBUG(waitCounter.tac());
            complete += countMessages;


            for(int idxRcv = 0 ; idxRcv < countMessages ; ++idxRcv){
                if( indexMessage[idxRcv] < nbMessagesToRecv ){
                    const int idxProc = status[idxRcv].MPI_SOURCE;
                    int nbLeaves;
                    (*recvBuffer[idxProc]) >> nbLeaves;
                    for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
                        MortonIndex leafIndex;
                        (*recvBuffer[idxProc]) >> leafIndex;
                        otherP2Ptree.createLeaf(leafIndex)->getSrc()->restore((*recvBuffer[idxProc]));
                    }
                    delete recvBuffer[idxProc];
                    recvBuffer[idxProc] = 0;
                }
            }
        }
        delete[] indexMessage;

        //////////////////////////////////////////////////////////
        // Computation P2P that need others data
        //////////////////////////////////////////////////////////

        FTRACE( FTrace::FRegion regionOtherTrace("Compute P2P Other", __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( computation2Counter.tic() );

        #pragma omp parallel
        {
            KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
            // There is a maximum of 26 neighbors
            ContainerClass* neighbors[27];
            MortonIndex indexesNeighbors[27];
            int indexArray[26];
            // Box limite
            const int limite = 1 << (this->OctreeHeight - 1);
            const int nbLeafToProceed = leafsNeedOtherData.getSize();

            #pragma omp for
            for(int idxLeafs = 0 ; idxLeafs < nbLeafToProceed ; ++idxLeafs){
                LeafData currentIter = leafsNeedOtherData[idxLeafs];

                // need the current particles and neighbors particles
                int counter = 0;
                memset( neighbors, 0, sizeof(ContainerClass*) * 27);

                // Take possible data
                const int nbNeigh = getNeighborsIndexesPeriodic(currentIter.cell->getCoordinate(), limite, indexesNeighbors, indexArray, periodicDirections);

                for(int idxNeigh = 0 ; idxNeigh < nbNeigh ; ++idxNeigh){
                    if(indexesNeighbors[idxNeigh] < intervals[idProcess].min || intervals[idProcess].max < indexesNeighbors[idxNeigh]){
                        ContainerClass*const hypotheticNeighbor = otherP2Ptree.getLeafSrc(indexesNeighbors[idxNeigh]);
                        if(hypotheticNeighbor){
                            neighbors[ indexArray[idxNeigh] ] = hypotheticNeighbor;
                            ++counter;
                        }
                    }
                }

                myThreadkernels.P2PRemote( currentIter.cell->getCoordinate(), currentIter.targets,
                                     currentIter.sources, neighbors, counter);
            }

        }

        for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
            delete sendBuffer[idxProc];
            delete recvBuffer[idxProc];
        }
        delete[] globalReceiveMap;
        delete[] leafsDataArray;

        FDEBUG(computation2Counter.tac());


        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation L2P + P2P : " << computationCounter.elapsed() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation P2P 2 : " << computation2Counter.elapsed() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Prepare P2P : " << prepareCounter.elapsed() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Gather P2P : " << gatherCounter.elapsed() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Wait : " << waitCounter.elapsed() << " s\n" );

    }



    int getNeighborsIndexesPeriodic(const FTreeCoordinate& center, const int boxLimite, MortonIndex indexes[26], int indexInArray[26], const int inDirection) const{
        if(center.getX() != 0 && center.getY() != 0 && center.getZ() != 0 &&
                center.getX() != boxLimite -1 && center.getY() != boxLimite -1  && center.getZ() != boxLimite -1 ){
            return getNeighborsIndexes(center, indexes, indexInArray);
        }

        int idxNeig = 0;

        const int startX = (TestPeriodicCondition(inDirection, DirMinusX) || center.getX() != 0 ?-1:0);
        const int endX = (TestPeriodicCondition(inDirection, DirPlusX) || center.getX() != boxLimite - 1 ?1:0);
        const int startY = (TestPeriodicCondition(inDirection, DirMinusY) || center.getY() != 0 ?-1:0);
        const int endY = (TestPeriodicCondition(inDirection, DirPlusY) || center.getY() != boxLimite - 1 ?1:0);
        const int startZ = (TestPeriodicCondition(inDirection, DirMinusZ) || center.getZ() != 0 ?-1:0);
        const int endZ = (TestPeriodicCondition(inDirection, DirPlusZ) || center.getZ() != boxLimite - 1 ?1:0);

        // We test all cells around
        for(int idxX = startX ; idxX <= endX ; ++idxX){
            for(int idxY = startY ; idxY <= endY ; ++idxY){
                for(int idxZ = startZ ; idxZ <= endZ ; ++idxZ){
                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){
                        FTreeCoordinate other(center.getX() + idxX,center.getY() + idxY,center.getZ() + idxZ);

                        if( other.getX() < 0 ) other.setX( other.getX() + boxLimite );
                        else if( boxLimite <= other.getX() ) other.setX( other.getX() - boxLimite );
                        if( other.getY() < 0 ) other.setY( other.getY() + boxLimite );
                        else if( boxLimite <= other.getY() ) other.setY( other.getY() - boxLimite );
                        if( other.getZ() < 0 ) other.setZ( other.getZ() + boxLimite );
                        else if( boxLimite <= other.getZ() ) other.setZ( other.getZ() - boxLimite );

                        indexes[idxNeig] = other.getMortonIndex(this->OctreeHeight - 1);
                        indexInArray[idxNeig] = ((idxX+1)*3  + (idxY+1))*3 + (idxZ+1);
                        ++idxNeig;
                    }
                }
            }
        }
        return idxNeig;
    }

    int getNeighborsIndexes(const FTreeCoordinate& center, MortonIndex indexes[26], int indexInArray[26]) const{
        int idxNeig = 0;

        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){
                        const FTreeCoordinate other(center.getX() + idxX,center.getY() + idxY,center.getZ() + idxZ);
                        indexes[idxNeig] = other.getMortonIndex(this->OctreeHeight - 1);
                        indexInArray[idxNeig] = ((idxX+1)*3  + (idxY+1))*3 + (idxZ+1);
                        ++idxNeig;
                    }
                }
            }
        }
        return idxNeig;
    }


    int getPeriodicInteractionNeighbors(const FTreeCoordinate& workingCell,const int inLevel, MortonIndex inNeighbors[189], int inNeighborsPosition[189], const int inDirection) const{

        // Then take each child of the parent's neighbors if not in directNeighbors
        // Father coordinate
        const FTreeCoordinate parentCell(workingCell.getX()>>1,workingCell.getY()>>1,workingCell.getZ()>>1);

        // Limite at parent level number of box (split by 2 by level)
        const int boxLimite = FMath::pow2(inLevel-1);

        // This is not on a border we can use normal interaction list method
        if( !(parentCell.getX() == 0 || parentCell.getY() == 0 || parentCell.getZ() == 0 ||
              parentCell.getX() == boxLimite - 1 || parentCell.getY() == boxLimite - 1 || parentCell.getZ() == boxLimite - 1 ) ) {
            return getInteractionNeighbors( workingCell, inLevel, inNeighbors, inNeighborsPosition);
        }

        const int startX =  (TestPeriodicCondition(inDirection, DirMinusX) || parentCell.getX() != 0 ?-1:0);
        const int endX =    (TestPeriodicCondition(inDirection, DirPlusX)  || parentCell.getX() != boxLimite - 1 ?1:0);
        const int startY =  (TestPeriodicCondition(inDirection, DirMinusY) || parentCell.getY() != 0 ?-1:0);
        const int endY =    (TestPeriodicCondition(inDirection, DirPlusY)  || parentCell.getY() != boxLimite - 1 ?1:0);
        const int startZ =  (TestPeriodicCondition(inDirection, DirMinusZ) || parentCell.getZ() != 0 ?-1:0);
        const int endZ =    (TestPeriodicCondition(inDirection, DirPlusZ)  || parentCell.getZ() != boxLimite - 1 ?1:0);

        int idxNeighbors = 0;
        // We test all cells around
        for(int idxX = startX ; idxX <= endX ; ++idxX){
            for(int idxY = startY ; idxY <= endY ; ++idxY){
                for(int idxZ = startZ ; idxZ <= endZ ; ++idxZ){
                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){

                        const FTreeCoordinate otherParent(parentCell.getX() + idxX,parentCell.getY() + idxY,parentCell.getZ() + idxZ);
                        FTreeCoordinate otherParentInBox(otherParent);
                        // periodic
                        if( otherParentInBox.getX() < 0 ){
                            otherParentInBox.setX( otherParentInBox.getX() + boxLimite );
                        }
                        else if( boxLimite <= otherParentInBox.getX() ){
                            otherParentInBox.setX( otherParentInBox.getX() - boxLimite );
                        }

                        if( otherParentInBox.getY() < 0 ){
                            otherParentInBox.setY( otherParentInBox.getY() + boxLimite );
                        }
                        else if( boxLimite <= otherParentInBox.getY() ){
                            otherParentInBox.setY( otherParentInBox.getY() - boxLimite );
                        }

                        if( otherParentInBox.getZ() < 0 ){
                            otherParentInBox.setZ( otherParentInBox.getZ() + boxLimite );
                        }
                        else if( boxLimite <= otherParentInBox.getZ() ){
                            otherParentInBox.setZ( otherParentInBox.getZ() - boxLimite );
                        }

                        const MortonIndex mortonOtherParent = otherParentInBox.getMortonIndex(inLevel-1);

                        // For each child
                        for(int idxCousin = 0 ; idxCousin < 8 ; ++idxCousin){
                            const int xdiff  = ((otherParent.getX()<<1) | ( (idxCousin>>2) & 1)) - workingCell.getX();
                            const int ydiff  = ((otherParent.getY()<<1) | ( (idxCousin>>1) & 1)) - workingCell.getY();
                            const int zdiff  = ((otherParent.getZ()<<1) | (idxCousin&1))         - workingCell.getZ();

                            // Test if it is a direct neighbor
                            if(FMath::Abs(xdiff) > 1 || FMath::Abs(ydiff) > 1 || FMath::Abs(zdiff) > 1){
                                // add to neighbors
                                inNeighborsPosition[idxNeighbors] = (((xdiff+3) * 7) + (ydiff+3)) * 7 + zdiff + 3;
                                inNeighbors[idxNeighbors++] = (mortonOtherParent << 3) | idxCousin;
                            }
                        }
                    }
                }
            }
        }

        return idxNeighbors;
    }

    int getInteractionNeighbors(const FTreeCoordinate& workingCell,const int inLevel, MortonIndex inNeighbors[189], int inNeighborsPosition[189]) const{

        // Then take each child of the parent's neighbors if not in directNeighbors
        // Father coordinate
        const FTreeCoordinate parentCell(workingCell.getX()>>1,workingCell.getY()>>1,workingCell.getZ()>>1);

        int idxNeighbors = 0;
        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){

                        const FTreeCoordinate otherParent(parentCell.getX() + idxX,parentCell.getY() + idxY,parentCell.getZ() + idxZ);

                        const MortonIndex mortonOtherParent = otherParent.getMortonIndex(inLevel-1);

                        // For each child
                        for(int idxCousin = 0 ; idxCousin < 8 ; ++idxCousin){
                            const int xdiff  = ((otherParent.getX()<<1) | ( (idxCousin>>2) & 1)) - workingCell.getX();
                            const int ydiff  = ((otherParent.getY()<<1) | ( (idxCousin>>1) & 1)) - workingCell.getY();
                            const int zdiff  = ((otherParent.getZ()<<1) | (idxCousin&1)) - workingCell.getZ();

                            // Test if it is a direct neighbor
                            if(FMath::Abs(xdiff) > 1 || FMath::Abs(ydiff) > 1 || FMath::Abs(zdiff) > 1){
                                // add to neighbors
                                inNeighborsPosition[idxNeighbors] = ((( (xdiff+3) * 7) + (ydiff+3))) * 7 + zdiff + 3;
                                inNeighbors[idxNeighbors++] = (mortonOtherParent << 3) | idxCousin;
                            }
                        }
                    }
                }
            }
        }

        return idxNeighbors;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Periodic levels = levels <= 0
    /////////////////////////////////////////////////////////////////////////////
public:

    /** This function process several M2M from level nbLevelsAboveRoot to level 0
      * and give the final result
      * @param result the cell at the last M2M
      * @param root the starting cell
      * @param startX the beginning of the index in x [0;endX]
      * @param endX the end of the index in x [startX;1]
      * @param startY the beginning of the index in y [0;endY]
      * @param endY the end of the index in y [startY;1]
      * @param startZ the beginning of the index in z [0;endZ]
      * @param endZ the end of the index in z [startZ;1]
      */
    void processTopM2MInIntervals( CellClass*const result, const CellClass& root, const int startX,
                              const int endX, const int startY, const int endY, const int startZ,
                              const int endZ){
        // allocate array
        CellClass*const cellsAtLevel = new CellClass[nbLevelsAboveRoot+2];
        // process by using other function
        processM2MInIntervals(cellsAtLevel,root,startX,endX,startY,endY,startZ,endZ);
        // copy result
        *result = cellsAtLevel[0];
        delete[] cellsAtLevel;
    }

    /** This function process several M2M from level nbLevelsAboveRoot to level 0
      * @param cellsAtLevel the intermediate results
      * @param root the starting cell
      * @param startX the beginning of the index in x [0;endX]
      * @param endX the end of the index in x [startX;1]
      * @param startY the beginning of the index in y [0;endY]
      * @param endY the end of the index in y [startY;1]
      * @param startZ the beginning of the index in z [0;endZ]
      * @param endZ the end of the index in z [startZ;1]
      */
    void  processM2MInIntervals( CellClass cellsAtLevel[], const CellClass& root, const int startX,
                              const int endX, const int startY, const int endY, const int startZ,
                              const int endZ){
        // start from the initial cell
        cellsAtLevel[nbLevelsAboveRoot+1] = root;
        // to create virtual children
        CellClass* virtualChild[8];
        // for all levels
        for(int idxLevel = nbLevelsAboveRoot ; idxLevel >= 0  ; --idxLevel){
            // reset children
            memset(virtualChild, 0, sizeof(CellClass*)*8);
            // fill the vector with previous result
            for(int idxX = startX ; idxX <= endX ; ++idxX){
                for(int idxY = startY ; idxY <= endY ; ++idxY){
                    for(int idxZ = startZ ; idxZ <= endZ ; ++idxZ){
                        virtualChild[childIndex(idxX,idxY,idxZ)] = &cellsAtLevel[idxLevel+1];
                    }
                }
            }
            // compute the M2M
            kernels[0]->M2M( &cellsAtLevel[idxLevel], virtualChild, idxLevel + 2);
        }
    }

    /** Fill an interactions neighbors with some intervals
      * @param neighbors the vector to fill
      * @param source the source cell to fill the vector
      * @param startX the beginning of the index in x [-3;0]
      * @param endX the end of the index in x  [0;3]
      * @param startY the beginning of the index in y [-3;0]
      * @param endY the end of the index in y [0;3]
      * @param startZ the beginning of the index in z [-3;0]
      * @param endZ the end of the index in z [0;3]
      * @return the number of position filled
      */
    int  fillM2LVectorFromIntervals(const CellClass* neighbors[343], const CellClass& source,
                     const int startX, const int endX, const int startY, const int endY,
                     const int startZ, const int endZ){
        int counter = 0;
        // for all x in interval
        for(int idxX = startX ; idxX <= endX ; ++idxX){
            // for all y in interval
            for(int idxY = startY ; idxY <= endY ; ++idxY){
                // for all z in interval
                for(int idxZ = startZ ; idxZ <= endZ ; ++idxZ){
                    // do not fill close neigbors
                    if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1 ){
                        neighbors[neighIndex(idxX,idxY,idxZ)] = &source;
                        ++counter;
                    }
                }
            }
        }
        // return the number of position filled
        return counter;
    }

    /** Get the index of a child (for the M2M and the L2L)
      * @param x the x position in the children  (from 0 to +1)
      * @param y the y position in the children  (from 0 to +1)
      * @param z the z position in the children  (from 0 to +1)
      * @return the index (from 0 to 7)
      */
    int childIndex(const int x, const int y, const int z) const {
        return (x<<2) | (y<<1) | z;
    }

    /** Get the index of a interaction neighbors (for M2L)
      * @param x the x position in the interactions (from -3 to +3)
      * @param y the y position in the interactions (from -3 to +3)
      * @param z the z position in the interactions (from -3 to +3)
      * @return the index (from 0 to 342)
      */
    int neighIndex(const int x, const int y, const int z) const {
        return (((x+3)*7) + (y+3))*7 + (z + 3);
    }

    /** To know how many times the box is repeated in each direction
      * -x +x -y +y -z +z
      * If the periodicy is not in all direction this number is unchanged
      * because it contains the theorical periodic box width for the
      * nbLevelsAboveRoot choosen
      */
    int theoricalRepetition() const {
        return nbLevelsAboveRoot == -1 ? 3 : 3 * (1<<(nbLevelsAboveRoot+1)) + 1;
    }

    /** To know the number of box repeated in each direction
      * @param min the number of repetition for -x,-y,-z
      * @param max the number of repetition for x,y,z
      * The mins value are contains between [-(theoricalRepetition-1 / 2); 0]
      * The maxs value are contains between [0;(theoricalRepetition-1 / 2)]
      */
    void repetitionsIntervals(FTreeCoordinate*const min, FTreeCoordinate*const max) const {
        const int halfRepeated = (theoricalRepetition()-1) /2;
        min->setPosition(-ifDir(DirMinusX,halfRepeated,0),-ifDir(DirMinusY,halfRepeated,0),
                         -ifDir(DirMinusZ,halfRepeated,0));
        max->setPosition(ifDir(DirPlusX,halfRepeated,0),ifDir(DirPlusY,halfRepeated,0),
                         ifDir(DirPlusZ,halfRepeated,0));
    }

    /** To get the number of repetition in all direction (x,y,z)
      * @return the number of box in all directions
      * Each value is between [1;theoricalRepetition]
      */
    FTreeCoordinate repetitions() const {
        const int halfRepeated = (theoricalRepetition()-1) /2;
        return FTreeCoordinate(ifDir(DirMinusX,halfRepeated,0) + ifDir(DirPlusX,halfRepeated,0) + 1,
                               ifDir(DirMinusY,halfRepeated,0) + ifDir(DirPlusY,halfRepeated,0) + 1,
                               ifDir(DirMinusZ,halfRepeated,0) + ifDir(DirPlusZ,halfRepeated,0) + 1);
    }

    /** This function has to be used to init the kernel with correct args
      * it return the box seen from a kernel point of view from the periodicity the user ask for
      * this is computed using the originalBoxWidth given in parameter
      * @param originalBoxWidth the real system size
      * @return the size the kernel should use
      */
    FReal extendedBoxWidth() const {
        return tree->getBoxWidth() * FReal(1<<offsetRealTree);
    }

    /** This function has to be used to init the kernel with correct args
      * it return the box cneter seen from a kernel point of view from the periodicity the user ask for
      * this is computed using the originalBoxWidth and originalBoxCenter given in parameter
      * @param originalBoxCenter the real system center
      * @param originalBoxWidth the real system size
      * @return the center the kernel should use
      */
    FPoint extendedBoxCenter() const {
        const FReal originalBoxWidth = tree->getBoxWidth();
        const FPoint originalBoxCenter = tree->getBoxCenter();
        const FReal offset = originalBoxWidth * FReal(1<<(offsetRealTree-1)) - originalBoxWidth/FReal(2.0);
        return FPoint( originalBoxCenter.getX() + offset,
                       originalBoxCenter.getY() + offset,
                       originalBoxCenter.getZ() + offset);
    }

    /** This function has to be used to init the kernel with correct args
      * it return the tree heigh seen from a kernel point of view from the periodicity the user ask for
      * this is computed using the originalTreeHeight given in parameter
      * @param originalTreeHeight the real tree heigh
      * @return the heigh the kernel should use
      */
    int extendedTreeHeight() const {
        // The real height
        return OctreeHeight + offsetRealTree;
    }

    /** To know if a direction is used
      * @param testDir a direction to test
      * @return true if the direction is used else false
      */
    bool usePerDir(const int testDir) const{
        return TestPeriodicCondition(periodicDirections , PeriodicCondition(testDir));
    }

    /** To enable quick test of the direction
      * @param testDir the direction to test
      * @param correctValue the value to return if direction is use
      * @param wrongValue the value to return if direction is not use
      * @return correctValue if testDir is used, else wrongValue
      */
    template <class T>
    int ifDir(const PeriodicCondition testDir, const T& correctValue, const T& wrongValue) const {
        return (periodicDirections & testDir ? correctValue : wrongValue);
    }

    /** Periodicity Core
      * This function is split in several part:
      * 1 - special case managment
      * There is nothing to do if nbLevelsAboveRoot == -1 and only
      * a M2L if nbLevelsAboveRoot == 0
      * 2 - if nbLevelsAboveRoot > 0
      * First we compute M2M and special M2M if needed for the border
      * Then the M2L by taking into account the periodicity directions
      * Then the border by using the precomputed M2M
      * Finally the L2L
      */
    void processPeriodicLevels(){
        /////////////////////////////////////////////////////
        // If nb level == -1 nothing to do
        if( nbLevelsAboveRoot == -1 ){
            return;
        }
        /////////////////////////////////////////////////////
        // if nb level == 0 only M2L at real root level
        if( nbLevelsAboveRoot == 0 ){
            CellClass rootUp;
            // compute the root
            rootUp = rootCellFromProc;

            // build fack M2L vector from -3/+3 x/y/z
            const CellClass* neighbors[343];
            memset(neighbors, 0, sizeof(CellClass*) * 343);
            int counter = 0;
            for(int idxX = ifDir(DirMinusX,-3,0) ; idxX <= ifDir(DirPlusX,3,0) ; ++idxX){
                for(int idxY = ifDir(DirMinusY,-3,0) ; idxY <= ifDir(DirPlusY,3,0) ; ++idxY){
                    for(int idxZ = ifDir(DirMinusZ,-3,0) ; idxZ <= ifDir(DirPlusZ,3,0) ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            neighbors[neighIndex(idxX,idxY,idxZ)] = &rootUp;
                            ++counter;
                        }
                    }
                }
            }
            // compute M2L
            CellClass rootDown;
            kernels[0]->M2L( &rootDown , neighbors, counter, 2);

            // put result in level 1
            rootCellFromProc = rootDown;

            return;
        }
        /////////////////////////////////////////////////////
        // in other situation, we have to compute M2L from 0 to nbLevelsAboveRoot
        // but also at nbLevelsAboveRoot +1 for the rest
        CellClass*const upperCells = new CellClass[nbLevelsAboveRoot+2];

        CellClass*const cellsXAxis = new CellClass[nbLevelsAboveRoot+2];
        CellClass*const cellsYAxis = new CellClass[nbLevelsAboveRoot+2];
        CellClass*const cellsZAxis = new CellClass[nbLevelsAboveRoot+2];
        CellClass*const cellsXYAxis = new CellClass[nbLevelsAboveRoot+2];
        CellClass*const cellsYZAxis = new CellClass[nbLevelsAboveRoot+2];
        CellClass*const cellsXZAxis = new CellClass[nbLevelsAboveRoot+2];


        // First M2M from level 1 to level 0
        upperCells[nbLevelsAboveRoot+1] = rootCellFromProc;

        // Then M2M from level 0 to level -LIMITE
        {
            CellClass* virtualChild[8];
            for(int idxLevel = nbLevelsAboveRoot ; idxLevel > 0  ; --idxLevel){
                FMemUtils::setall(virtualChild,&upperCells[idxLevel+1],8);
                kernels[0]->M2M( &upperCells[idxLevel], virtualChild, idxLevel + 2);
            }

            // Cells on the axis of the center should be computed separatly.

            if( usePerDir(DirMinusX) && usePerDir(DirMinusY) ){
                FMemUtils::copyall(cellsZAxis,upperCells,nbLevelsAboveRoot+2);
            }
            else{
                 processM2MInIntervals(cellsZAxis,upperCells[nbLevelsAboveRoot+1],
                                    ifDir(DirMinusX,0,1),1,ifDir(DirMinusY,0,1),1,0,1);
            }
            if( usePerDir(DirMinusX) && usePerDir(DirMinusZ) ){
                FMemUtils::copyall(cellsYAxis,upperCells,nbLevelsAboveRoot+2);
            }
            else{
                 processM2MInIntervals(cellsYAxis,upperCells[nbLevelsAboveRoot+1],
                                    ifDir(DirMinusX,0,1),1,0,1,ifDir(DirMinusZ,0,1),1);
            }
            if( usePerDir(DirMinusY) && usePerDir(DirMinusZ) ){
                FMemUtils::copyall(cellsXAxis,upperCells,nbLevelsAboveRoot+2);
            }
            else{
                 processM2MInIntervals(cellsXAxis,upperCells[nbLevelsAboveRoot+1],
                                    0,1,ifDir(DirMinusY,0,1),1,ifDir(DirMinusZ,0,1),1);
            }

            // Then cells on the spaces should be computed separatly

            if( !usePerDir(DirMinusX) ){
                 processM2MInIntervals(cellsYZAxis,upperCells[nbLevelsAboveRoot+1],1,1,0,1,0,1);
            }
            else {
                FMemUtils::copyall(cellsYZAxis,upperCells,nbLevelsAboveRoot+2);
            }
            if( !usePerDir(DirMinusY) ){
                 processM2MInIntervals(cellsXZAxis,upperCells[nbLevelsAboveRoot+1],0,1,1,1,0,1);
            }
            else {
                FMemUtils::copyall(cellsXZAxis,upperCells,nbLevelsAboveRoot+2);
            }
            if( !usePerDir(DirMinusZ) ){
                 processM2MInIntervals(cellsXYAxis,upperCells[nbLevelsAboveRoot+1],0,1,0,1,1,1);
            }
            else {
                FMemUtils::copyall(cellsXYAxis,upperCells,nbLevelsAboveRoot+2);
            }

        }

        // Then M2L at all level
        {
            CellClass* positionedCells[343];
            memset(positionedCells, 0, 343 * sizeof(CellClass**));

            for(int idxX = ifDir(DirMinusX,-3,0) ; idxX <= ifDir(DirPlusX,3,0) ; ++idxX){
                for(int idxY = ifDir(DirMinusY,-3,0) ; idxY <= ifDir(DirPlusY,3,0) ; ++idxY){
                    for(int idxZ = ifDir(DirMinusZ,-3,0) ; idxZ <= ifDir(DirPlusZ,3,0) ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            if(idxX == 0 && idxY == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsZAxis;
                            }
                            else if(idxX == 0 && idxZ == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsYAxis;
                            }
                            else if(idxY == 0 && idxZ == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsXAxis;
                            }
                            else if(idxX == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsYZAxis;
                            }
                            else if(idxY == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsXZAxis;
                            }
                            else if(idxZ == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsXYAxis;
                            }
                            else{
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = upperCells;
                            }
                        }
                    }
                }
            }

            // We say that we are in the child index 0
            // So we can compute one time the relative indexes
            const CellClass* neighbors[343];
            memset(neighbors, 0, sizeof(CellClass*) * 343);
            int counter = 0;
            for(int idxX = ifDir(DirMinusX,-3,0) ; idxX <= ifDir(DirPlusX,2,0) ; ++idxX){
                for(int idxY = ifDir(DirMinusY,-3,0) ; idxY <= ifDir(DirPlusY,2,0) ; ++idxY){
                    for(int idxZ = ifDir(DirMinusZ,-3,0) ; idxZ <= ifDir(DirPlusZ,2,0) ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            neighbors[neighIndex(idxX,idxY,idxZ)] = reinterpret_cast<const CellClass*>(~0);
                            ++counter;
                        }
                    }
                }
            }

            for(int idxLevel = nbLevelsAboveRoot + 1 ; idxLevel > 1 ; --idxLevel ){
                for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
                    if(neighbors[idxNeigh]){
                        neighbors[idxNeigh] = &positionedCells[idxNeigh][idxLevel];
                    }
                }
                kernels[0]->M2L( &upperCells[idxLevel] , neighbors, counter , idxLevel + 2);
            }

            memset(neighbors, 0, sizeof(CellClass*) * 343);
            counter = 0;
            for(int idxX = ifDir(DirMinusX,-2,0) ; idxX <= ifDir(DirPlusX,3,0) ; ++idxX){
                for(int idxY = ifDir(DirMinusY,-2,0) ; idxY <= ifDir(DirPlusY,3,0) ; ++idxY){
                    for(int idxZ = ifDir(DirMinusZ,-2,0) ; idxZ <= ifDir(DirPlusZ,3,0) ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            const int index = neighIndex(idxX,idxY,idxZ);
                            neighbors[index] = &positionedCells[index][1];
                            ++counter;
                        }
                    }
                }
            }
            kernels[0]->M2L( &upperCells[1] , neighbors, 189, 3);
        }

        { // compute the border
            if( usePerDir(AllDirs) ){
                CellClass leftborder, bottomborder, frontborder, angleborderlb,
                        angleborderfb, angleborderlf, angleborder;

                const CellClass* neighbors[343];
                memset(neighbors, 0, sizeof(CellClass*) * 343);
                int counter = 0;

                processTopM2MInIntervals( &leftborder, upperCells[nbLevelsAboveRoot+1],     1,1 , 0,1 , 0,1);
                counter +=  fillM2LVectorFromIntervals(neighbors, leftborder,     -2,-2 , -1,1,  -1,1 );
                processTopM2MInIntervals( &bottomborder, upperCells[nbLevelsAboveRoot+1],   0,1 , 0,1 , 1,1);
                counter +=  fillM2LVectorFromIntervals(neighbors, bottomborder,   -1,1  , -1,1,  -2,-2);
                processTopM2MInIntervals( &frontborder, upperCells[nbLevelsAboveRoot+1],    0,1 , 1,1 , 0,1);
                counter +=  fillM2LVectorFromIntervals(neighbors, frontborder,    -1,1  , -2,-2, -1,1 );
                processTopM2MInIntervals( &angleborderlb, upperCells[nbLevelsAboveRoot+1],  1,1 , 0,1 , 1,1);
                counter +=  fillM2LVectorFromIntervals(neighbors, angleborderlb,  -2,-2 , -1,1,  -2,-2);
                processTopM2MInIntervals( &angleborderfb, upperCells[nbLevelsAboveRoot+1],  0,1 , 1,1 , 1,1);
                counter +=  fillM2LVectorFromIntervals(neighbors, angleborderfb,  -1,1 ,  -2,-2, -2,-2);
                processTopM2MInIntervals( &angleborderlf, upperCells[nbLevelsAboveRoot+1],  1,1 , 1,1 , 0,1);
                counter +=  fillM2LVectorFromIntervals(neighbors, angleborderlf,  -2,-2 , -2,-2, -1,1 );
                processTopM2MInIntervals( &angleborder, upperCells[nbLevelsAboveRoot+1],    1,1 , 1,1 , 1,1);
                counter +=  fillM2LVectorFromIntervals(neighbors, angleborder,    -2,-2 , -2,-2, -2,-2);


                kernels[0]->M2L( &upperCells[0] , neighbors, counter, 2);

                CellClass* virtualChild[8];
                memset(virtualChild, 0, sizeof(CellClass*) * 8);
                virtualChild[childIndex(0,0,0)] = &upperCells[1];
                kernels[0]->L2L( &upperCells[0], virtualChild, 2);
            }
            else {
                CellClass*const leftborder = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const bottomborder = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const frontborder = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const angleborderlb = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const angleborderfb = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const angleborderlf = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const angleborder = new CellClass[nbLevelsAboveRoot+2];

                 processM2MInIntervals( leftborder,   upperCells[nbLevelsAboveRoot+1],     1,1 , 0,1 , 0,1);
                 processM2MInIntervals( bottomborder, upperCells[nbLevelsAboveRoot+1],   0,1 , 0,1 , 1,1);
                 processM2MInIntervals( frontborder,  upperCells[nbLevelsAboveRoot+1],    0,1 , 1,1 , 0,1);
                 processM2MInIntervals( angleborderlb,upperCells[nbLevelsAboveRoot+1],  1,1 , 0,1 , 1,1);
                 processM2MInIntervals( angleborderfb,upperCells[nbLevelsAboveRoot+1],  0,1 , 1,1 , 1,1);
                 processM2MInIntervals( angleborderlf,upperCells[nbLevelsAboveRoot+1],  1,1 , 1,1 , 0,1);
                 processM2MInIntervals( angleborder,  upperCells[nbLevelsAboveRoot+1],    1,1 , 1,1 , 1,1);

                const CellClass* neighbors[343];
                memset(neighbors, 0, sizeof(CellClass*) * 343);
                int counter = 0;

                if(usePerDir(DirMinusX) && usePerDir(DirMinusY) && usePerDir(DirMinusZ)){
                    neighbors[neighIndex(-2,-1,-1)] = &leftborder[0];
                    neighbors[neighIndex(-1,-2,-1)] = &frontborder[0];
                    neighbors[neighIndex(-1,-1,-2)] = &bottomborder[0];
                    neighbors[neighIndex(-2,-2,-1)] = &angleborderlf[0];
                    neighbors[neighIndex(-2,-1,-2)] = &angleborderlb[0];
                    neighbors[neighIndex(-1,-2,-2)] = &angleborderfb[0];
                    neighbors[neighIndex(-2,-2,-2)] = &angleborder[0];
                    counter += 7;
                }
                if(usePerDir(DirMinusX) && usePerDir(DirPlusY) && usePerDir(DirMinusZ)){
                    neighbors[neighIndex(-2,1,-1)] = &leftborder[0];
                    neighbors[neighIndex(-1,1,-2)] = &bottomborder[0];
                    neighbors[neighIndex(-2,1,-2)] = &angleborderlb[0];
                    counter += 3;
                }
                if(usePerDir(DirMinusX) && usePerDir(DirMinusY) && usePerDir(DirPlusZ)){
                    neighbors[neighIndex(-2,-1,1)] = &leftborder[0];
                    neighbors[neighIndex(-1,-2,1)] = &frontborder[0];
                    neighbors[neighIndex(-2,-2,1)] = &angleborderlf[0];
                    counter += 3;
                }
                if(usePerDir(DirMinusX) && usePerDir(DirPlusY) && usePerDir(DirPlusZ)){
                    neighbors[neighIndex(-2,1,1)] = &leftborder[0];
                    counter += 1;
                }
                if(usePerDir(DirPlusX) && usePerDir(DirMinusY) && usePerDir(DirMinusZ)){
                    neighbors[neighIndex(1,-2,-1)] = &frontborder[0];
                    neighbors[neighIndex(1,-1,-2)] = &bottomborder[0];
                    neighbors[neighIndex(1,-2,-2)] = &angleborderfb[0];
                    counter += 3;
                }
                if(usePerDir(DirPlusX) && usePerDir(DirMinusY) && usePerDir(DirPlusZ)){
                    neighbors[neighIndex(1,-2,1)] = &frontborder[0];
                    counter += 1;
                }
                if(usePerDir(DirPlusX) && usePerDir(DirPlusY) && usePerDir(DirMinusZ)){
                    neighbors[neighIndex(1,1,-2)] = &bottomborder[0];
                    counter += 1;
                }

                CellClass centerXFace;
                CellClass centerYFace;
                CellClass centerZFace;
                CellClass centerXZFace;
                CellClass centerXYFace;
                CellClass centerYXFace;
                CellClass centerYZFace;
                CellClass centerZXFace;
                CellClass centerZYFace;
                CellClass angleXZFace;
                CellClass angleXYFace;
                CellClass angleYZFace;

                if(usePerDir(DirMinusX)){
                    {
                        CellClass* virtualChild[8];
                        memset(virtualChild, 0, 8*sizeof(CellClass*));
                        if(usePerDir(DirPlusY) && usePerDir(DirPlusZ))  virtualChild[childIndex(1,1,1)] = &leftborder[1];
                        if(usePerDir(DirMinusY) && usePerDir(DirPlusZ)) virtualChild[childIndex(1,0,1)] = &leftborder[1];
                        else if( usePerDir(DirPlusZ))                   virtualChild[childIndex(1,0,1)] = &angleborderlf[1];
                        if(usePerDir(DirPlusY) && usePerDir(DirMinusZ)) virtualChild[childIndex(1,1,0)] = &leftborder[1];
                        else if(usePerDir(DirPlusY) )                   virtualChild[childIndex(1,1,0)] = &angleborderlb[1];

                        if(usePerDir(DirMinusY) && usePerDir(DirMinusZ)) virtualChild[childIndex(1,0,0)] = &leftborder[1];
                        else if(usePerDir(DirMinusZ))                    virtualChild[childIndex(1,0,0)] = &angleborderlf[1];
                        else if(usePerDir(DirMinusY))                    virtualChild[childIndex(1,0,0)] = &angleborderlb[1];
                        else                                             virtualChild[childIndex(1,0,0)] = &angleborder[1];

                        kernels[0]->M2M( &centerXFace, virtualChild, 2);
                        neighbors[neighIndex(-2,0,0)] = &centerXFace;
                        counter += 1;
                    }
                    if(usePerDir(DirMinusZ) || usePerDir(DirPlusZ)){
                        if(usePerDir(DirY)){
                            centerXZFace = leftborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusY)){
                                virtualChild[childIndex(1,1,0)] = &leftborder[1];
                                virtualChild[childIndex(1,1,1)] = &leftborder[1];
                            }
                            if(usePerDir(DirMinusY)){
                                virtualChild[childIndex(1,0,0)] = &leftborder[1];
                                virtualChild[childIndex(1,0,1)] = &leftborder[1];
                            }
                            else{
                                virtualChild[childIndex(1,0,0)] = &angleborderlf[1];
                                virtualChild[childIndex(1,0,1)] = &angleborderlf[1];
                            }
                            kernels[0]->M2M( &centerXZFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusZ) ){
                            neighbors[neighIndex(-2,0,1)] = &centerXZFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusZ) ){
                            neighbors[neighIndex(-2,0,-1)] = &centerXZFace;
                            counter += 1;

                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusY)){
                                virtualChild[childIndex(1,1,1)] = &angleborderlb[1];
                            }
                            if(usePerDir(DirMinusY)){
                                virtualChild[childIndex(1,0,1)] = &angleborderlb[1];
                            }
                            else{
                                virtualChild[childIndex(1,0,1)] = &angleborder[1];
                            }
                            kernels[0]->M2M( &angleXZFace, virtualChild, 2);

                            neighbors[neighIndex(-2,0,-2)] = &angleXZFace;
                            counter += 1;
                        }
                    }
                    if(usePerDir(DirMinusY) || usePerDir(DirPlusY)){
                        if(usePerDir(DirZ)){
                            centerXYFace = leftborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusZ)){
                                virtualChild[childIndex(1,0,1)] = &leftborder[1];
                                virtualChild[childIndex(1,1,1)] = &leftborder[1];
                            }
                            if(usePerDir(DirMinusZ)){
                                virtualChild[childIndex(1,0,0)] = &leftborder[1];
                                virtualChild[childIndex(1,1,0)] = &leftborder[1];
                            }
                            else{
                                virtualChild[childIndex(1,0,0)] = &angleborderlb[1];
                                virtualChild[childIndex(1,1,0)] = &angleborderlb[1];
                            }
                            kernels[0]->M2M( &centerXYFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusY) ){
                            neighbors[neighIndex(-2,1,0)] = &centerXYFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusY) ){
                            neighbors[neighIndex(-2,-1,0)] = &centerXYFace;
                            counter += 1;

                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusZ)){
                                virtualChild[childIndex(1,1,1)] = &angleborderlf[1];
                            }
                            if(usePerDir(DirMinusZ)){
                                virtualChild[childIndex(1,1,0)] = &angleborderlf[1];
                            }
                            else{
                                virtualChild[childIndex(1,1,0)] = &angleborder[1];
                            }
                            kernels[0]->M2M( &angleXYFace, virtualChild, 2);

                            neighbors[neighIndex(-2,-2,0)] = &angleXYFace;
                            counter += 1;
                        }
                    }
                }
                if(usePerDir(DirMinusY)){
                    {
                        CellClass* virtualChild[8];
                        memset(virtualChild, 0, 8*sizeof(CellClass*));
                        if(usePerDir(DirPlusX) && usePerDir(DirPlusZ))  virtualChild[childIndex(1,1,1)] = &frontborder[1];
                        if(usePerDir(DirMinusX) && usePerDir(DirPlusZ)) virtualChild[childIndex(0,1,1)] = &frontborder[1];
                        else if(usePerDir(DirPlusZ))                    virtualChild[childIndex(0,1,1)] = &angleborderlf[1];
                        if(usePerDir(DirPlusX) && usePerDir(DirMinusZ)) virtualChild[childIndex(1,1,0)] = &frontborder[1];
                        else if(usePerDir(DirPlusX))                    virtualChild[childIndex(1,1,0)] = &angleborderfb[1];

                        if(usePerDir(DirMinusX) && usePerDir(DirMinusZ)) virtualChild[childIndex(0,1,0)] = &frontborder[1];
                        else if(usePerDir(DirMinusZ))                    virtualChild[childIndex(0,1,0)] = &angleborderlf[1];
                        else if(usePerDir(DirMinusX))                    virtualChild[childIndex(0,1,0)] = &angleborderfb[1];
                        else                                             virtualChild[childIndex(0,1,0)] = &angleborder[1];

                        kernels[0]->M2M( &centerYFace, virtualChild, 2);
                        neighbors[neighIndex(0,-2,0)] = &centerYFace;
                        counter += 1;
                    }
                    if(usePerDir(DirMinusZ) || usePerDir(DirPlusZ)){
                        if(usePerDir(DirX)){
                            centerYZFace = frontborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusX)){
                                virtualChild[childIndex(1,1,0)] = &frontborder[1];
                                virtualChild[childIndex(1,1,1)] = &frontborder[1];
                            }
                            if(usePerDir(DirMinusX)){
                                virtualChild[childIndex(0,1,0)] = &frontborder[1];
                                virtualChild[childIndex(0,1,1)] = &frontborder[1];
                            }
                            else{
                                virtualChild[childIndex(0,1,0)] = &angleborderlf[1];
                                virtualChild[childIndex(0,1,1)] = &angleborderlf[1];
                            }
                            kernels[0]->M2M( &centerYZFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusZ) ){
                            neighbors[neighIndex(0,-2,1)] = &centerYZFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusZ) ){
                            neighbors[neighIndex(0,-2,-1)] = &centerYZFace;
                            counter += 1;

                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusX)){
                                virtualChild[childIndex(1,1,1)] = &angleborderfb[1];
                            }
                            if(usePerDir(DirMinusX)){
                                virtualChild[childIndex(0,1,1)] = &angleborderfb[1];
                            }
                            else{
                                virtualChild[childIndex(0,1,1)] = &angleborder[1];
                            }
                            kernels[0]->M2M( &angleYZFace, virtualChild, 2);

                            neighbors[neighIndex(0,-2,-2)] = &angleYZFace;
                            counter += 1;
                        }
                    }
                    if(usePerDir(DirMinusX) || usePerDir(DirPlusX)){
                        if(usePerDir(DirZ)){
                            centerYXFace = frontborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusZ)){
                                virtualChild[childIndex(0,1,1)] = &frontborder[1];
                                virtualChild[childIndex(1,1,1)] = &frontborder[1];
                            }
                            if(usePerDir(DirMinusZ)){
                                virtualChild[childIndex(0,1,0)] = &frontborder[1];
                                virtualChild[childIndex(1,1,0)] = &frontborder[1];
                            }
                            else{
                                virtualChild[childIndex(0,1,0)] = &angleborderfb[1];
                                virtualChild[childIndex(1,1,0)] = &angleborderfb[1];
                            }
                            kernels[0]->M2M( &centerYXFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusX) ){
                            neighbors[neighIndex(1,-2,0)] = &centerYXFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusX) ){
                            neighbors[neighIndex(-1,-2,0)] = &centerYXFace;
                            counter += 1;
                        }
                    }
                }
                if(usePerDir(DirMinusZ)){
                    {
                        CellClass* virtualChild[8];
                        memset(virtualChild, 0, 8*sizeof(CellClass*));
                        if(usePerDir(DirPlusX) && usePerDir(DirPlusY))  virtualChild[childIndex(1,1,1)] = &bottomborder[1];
                        if(usePerDir(DirMinusX) && usePerDir(DirPlusY)) virtualChild[childIndex(0,1,1)] = &bottomborder[1];
                        else if( usePerDir(DirPlusY))                   virtualChild[childIndex(0,1,1)] = &angleborderlb[1];
                        if(usePerDir(DirPlusX) && usePerDir(DirMinusY)) virtualChild[childIndex(1,0,1)] = &bottomborder[1];
                        else if(usePerDir(DirPlusX) )                   virtualChild[childIndex(1,0,1)] = &angleborderfb[1];

                        if(usePerDir(DirMinusX) && usePerDir(DirMinusY)) virtualChild[childIndex(0,0,1)] = &bottomborder[1];
                        else if(usePerDir(DirMinusY))                    virtualChild[childIndex(0,0,1)] = &angleborderfb[1];
                        else if(usePerDir(DirMinusX))                    virtualChild[childIndex(0,0,1)] = &angleborderlb[1];
                        else                                             virtualChild[childIndex(0,0,1)] = &angleborder[1];

                        kernels[0]->M2M( &centerZFace, virtualChild, 2);
                        neighbors[neighIndex(0,0,-2)] = &centerZFace;
                        counter += 1;
                    }
                    if(usePerDir(DirMinusY) || usePerDir(DirPlusY)){
                        if(usePerDir(DirX)){
                            centerZYFace = bottomborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusX)){
                                virtualChild[childIndex(1,0,1)] = &bottomborder[1];
                                virtualChild[childIndex(1,1,1)] = &bottomborder[1];
                            }
                            if(usePerDir(DirMinusX)){
                                virtualChild[childIndex(0,0,1)] = &bottomborder[1];
                                virtualChild[childIndex(0,1,1)] = &bottomborder[1];
                            }
                            else{
                                virtualChild[childIndex(0,0,1)] = &angleborderlb[1];
                                virtualChild[childIndex(0,1,1)] = &angleborderlb[1];
                            }
                            kernels[0]->M2M( &centerZYFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusY) ){
                            neighbors[neighIndex(0,1,-2)] = &centerZYFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusY) ){
                            neighbors[neighIndex(0,-1,-2)] = &centerZYFace;
                            counter += 1;
                        }
                    }
                    if(usePerDir(DirMinusX) || usePerDir(DirPlusX)){
                        if(usePerDir(DirY)){
                            centerZXFace = bottomborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusY)){
                                virtualChild[childIndex(0,1,1)] = &bottomborder[1];
                                virtualChild[childIndex(1,1,1)] = &bottomborder[1];
                            }
                            if(usePerDir(DirMinusY)){
                                virtualChild[childIndex(0,0,1)] = &bottomborder[1];
                                virtualChild[childIndex(1,0,1)] = &bottomborder[1];
                            }
                            else{
                                virtualChild[childIndex(0,0,1)] = &angleborderlf[1];
                                virtualChild[childIndex(0,1,1)] = &angleborderlf[1];
                            }
                            kernels[0]->M2M( &centerZXFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusX) ){
                            neighbors[neighIndex(1,0,-2)] = &centerZXFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusX) ){
                            neighbors[neighIndex(-1,0,-2)] = &centerZXFace;
                            counter += 1;
                        }
                    }
                }

                // M2L for border
                kernels[0]->M2L( &upperCells[0] , neighbors, counter, 2);

                // L2L from border M2L to top of tree
                CellClass* virtualChild[8];
                memset(virtualChild, 0, sizeof(CellClass*) * 8);
                virtualChild[childIndex(0,0,0)] = &upperCells[1];
                kernels[0]->L2L( &upperCells[0], virtualChild, 2);

                // dealloc
                delete[] leftborder;
                delete[] bottomborder;
                delete[] frontborder;
                delete[] angleborderlb;
                delete[] angleborderfb;
                delete[] angleborderlf;
                delete[] angleborder;
            }

        }

        // Finally L2L until level 0
        {
            CellClass* virtualChild[8];
            memset(virtualChild, 0, sizeof(CellClass*) * 8);
            for(int idxLevel = 1 ; idxLevel <= nbLevelsAboveRoot  ; ++idxLevel){
                virtualChild[childIndex(1,1,1)] = &upperCells[idxLevel+1];
                kernels[0]->L2L( &upperCells[idxLevel], virtualChild, idxLevel + 2);
            }
        }

        // L2L from 0 to level 1
        rootCellFromProc = upperCells[nbLevelsAboveRoot+1];

        delete[] upperCells;

        delete[] cellsXAxis;
        delete[] cellsYAxis;
        delete[] cellsZAxis;
        delete[] cellsXYAxis;
        delete[] cellsYZAxis;
        delete[] cellsXZAxis;
    }
};






#endif //FFMMALGORITHMTHREAD_HPP


