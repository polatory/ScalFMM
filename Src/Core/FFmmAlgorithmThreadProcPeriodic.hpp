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


#include "../Utils/FAssert.hpp"
#include "../Utils/FLog.hpp"
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

#include "FCoreCommon.hpp"

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
class FFmmAlgorithmThreadProcPeriodic : public FAbstractAlgorithm {

    static const int MaxSizePerCell = 2048;

    OctreeClass* const tree;                 //< The octree to work on
    KernelClass** kernels;                   //< The kernels    

    const FMpi::FComm& comm;                 //< MPI comm

    CellClass rootCellFromProc;     //< root of tree needed by the periodicity
    const int nbLevelsAboveRoot;    //< The nb of level the user ask to go above the tree (>= -1)
    const int offsetRealTree;       //< nbLevelsAboveRoot GetFackLevel

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
                                    const int inUpperLevel = 2)
        : tree(inTree) , kernels(0), comm(inComm), nbLevelsAboveRoot(inUpperLevel), offsetRealTree(GetFackLevel(inUpperLevel)),
           numberOfLeafs(0),
          MaxThreads(omp_get_max_threads()), nbProcess(inComm.processCount()), idProcess(inComm.processId()),
          OctreeHeight(tree->getHeight()),intervals(new Interval[inComm.processCount()]),
          workingIntervalsPerLevel(new Interval[inComm.processCount() * tree->getHeight()]) {

        FAssertLF(tree, "tree cannot be null");
        FAssertLF(-1 <= inUpperLevel, "inUpperLevel cannot be < -1");

        FLOG(FLog::Controller << "FFmmAlgorithmThreadProcPeriodic\n");
        FLOG(FLog::Controller << "Max threads = "  << MaxThreads << ", Procs = " << nbProcess << ", I am " << idProcess << ".\n");
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
    void execute(const unsigned operationsToProceed = FFmmNearAndFarFields){
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
            FAssertLF(iterArray, "iterArray bad alloc");

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
        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if((operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P)) directPass();

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
        FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
        FLOG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        int leafs = 0;
        do{
            iterArray[leafs++] = octreeIterator;
        } while(octreeIterator.moveRight());

        FLOG(FTic computationCounter);
        #pragma omp parallel
        {
            KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for nowait
            for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentListSrc());
            }
        }
        FLOG(computationCounter.tac());


        FLOG( FLog::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );

    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FLOG( FLog::Controller.write("\tStart Upward Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);
        FLOG(FTic prepareCounter);
        FLOG(FTic waitCounter);

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

            FLOG(prepareCounter.tic());
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
            FLOG(prepareCounter.tac());
            FTRACE( regionTrace.end() );

            // Compute
            const int endIndex = (hasToReceive?numberOfCells-1:numberOfCells);
            FLOG(computationCounter.tic());
            #pragma omp parallel
            {
                KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
                #pragma omp for nowait
                for( int idxCell = cellsToSend + 1 ; idxCell < endIndex ; ++idxCell){
                    myThreadkernels.M2M( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), fackLevel);
                }
            }
            FLOG(computationCounter.tac());

            // Are we sending or waiting anything?
            if(iterRequests){
                FLOG(waitCounter.tic());
                MPI_Waitall( iterRequests, requests, status);
                FLOG(waitCounter.tac());

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

                            FAssertLF(!currentChild[position], "Already has a cell here");

                            recvBufferCells[position].deserializeUp(recvBuffer);
                            currentChild[position] = (CellClass*) &recvBufferCells[position];

                            state >>= 1;
                            ++position;
                        }
                    }

                    // Finally compute
                    FLOG(computationCounter.tic());
                    (*kernels[0]).M2M( iterArray[numberOfCells - 1].getCurrentCell() , currentChild, fackLevel);
                    FLOG(computationCounter.tac());


                    firstProcThatSend = endProcThatSend - 1;
                }
            }

            sendBuffer.reset();
            recvBuffer.seek(0);
        }


        FLOG( FLog::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        FLOG( FLog::Controller << "\t\t Prepare : " << prepareCounter.cumulated() << " s\n" );
        FLOG( FLog::Controller << "\t\t Wait : " << waitCounter.cumulated() << " s\n" );

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
                        FAssertLF(!currentChild[position], "Already has a cell here");

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

        FLOG( FLog::Controller.write("\tStart Downward Pass (M2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);
        FLOG(FTic sendCounter);
        FLOG(FTic receiveCounter);
        FLOG(FTic prepareCounter);
        FLOG(FTic gatherCounter);

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
            FLOG(prepareCounter.tic());

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
                                                                        neighborsIndexes, neighborsPosition, AllDirs);

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
            FLOG(prepareCounter.tac());

            delete[] alreadySent;
        }

        //////////////////////////////////////////////////////////////////
        // Gather this information
        //////////////////////////////////////////////////////////////////

        FLOG(gatherCounter.tic());
        // All process say to each others
        // what the will send to who
        int*const globalReceiveMap = new int[nbProcess * nbProcess * OctreeHeight];
        memset(globalReceiveMap, 0, sizeof(int) * nbProcess * nbProcess * OctreeHeight);
        FMpi::MpiAssert( MPI_Allgather( indexToSend, nbProcess * OctreeHeight, MPI_INT, globalReceiveMap, nbProcess * OctreeHeight, MPI_INT, comm.getComm()),  __LINE__ );
        FLOG(gatherCounter.tac());

        //////////////////////////////////////////////////////////////////
        // Send and receive for real
        //////////////////////////////////////////////////////////////////

        FLOG(sendCounter.tic());
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
        FLOG(sendCounter.tac());

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

                FLOG(computationCounter.tic());
                #pragma omp parallel
                {
                    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                    const CellClass* neighbors[343];

                    #pragma omp for  schedule(dynamic) nowait
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        const int counter = tree->getPeriodicInteractionNeighbors(neighbors, iterArray[idxCell].getCurrentGlobalCoordinate(),idxLevel, AllDirs);
                        if(counter) myThreadkernels->M2L( iterArray[idxCell].getCurrentCell() , neighbors, counter, fackLevel);
                    }

                    myThreadkernels->finishedLevelM2L(fackLevel);
                }
                FLOG(computationCounter.tac());
            }
        }

        //////////////////////////////////////////////////////////////////
        // Wait received data and compute
        //////////////////////////////////////////////////////////////////

        // Wait to receive every things (and send every things)
        MPI_Waitall(iterRequest, requests, status);

        {
            FTRACE( FTrace::FRegion regionTrace("Compute Received data", __FUNCTION__ , __FILE__ , __LINE__) );
            FLOG(receiveCounter.tic());
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
                FLOG(computationCounter.tic());
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
                                                                         neighborsIndex, neighborsPosition, AllDirs);
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
                FLOG(computationCounter.tac());
            }
            FLOG(receiveCounter.tac());
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

        FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        FLOG( FLog::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
        FLOG( FLog::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
        FLOG( FLog::Controller << "\t\t Gather : " << gatherCounter.cumulated() << " s\n" );
        FLOG( FLog::Controller << "\t\t Prepare : " << prepareCounter.cumulated() << " s\n" );
    }

    //////////////////////////////////////////////////////////////////
    // ---------------- L2L ---------------
    //////////////////////////////////////////////////////////////////

    void downardPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FLOG( FLog::Controller.write("\tStart Downward Pass (L2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);
        FLOG(FTic prepareCounter);
        FLOG(FTic waitCounter);

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

            FLOG(prepareCounter.tic());

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
            FLOG(prepareCounter.tac());

            FLOG(computationCounter.tic());
            #pragma omp parallel
            {
                KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
                #pragma omp for nowait
                for(int idxCell = firstCellWork + 1 ; idxCell < numberOfCells ; ++idxCell){
                    myThreadkernels.L2L( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), fackLevel);
                }
            }
            FLOG(computationCounter.tac());

            // are we sending or receiving?
            if(iterRequests){
                // process
                FLOG(waitCounter.tic());
                MPI_Waitall( iterRequests, requests, status);
                FLOG(waitCounter.tac());

                if(needToRecv){
                    // Need to compute
                    FLOG(computationCounter.tic());
                    iterArray[firstCellWork].getCurrentCell()->deserializeDown(recvBuffer);

                    kernels[0]->L2L( iterArray[firstCellWork].getCurrentCell() , iterArray[firstCellWork].getCurrentChild(), fackLevel);
                    FLOG(computationCounter.tac());
                }
            }

            sendBuffer.reset();
            recvBuffer.seek(0);
        }

        delete[] requests;
        delete[] status;

        FLOG( FLog::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        FLOG( FLog::Controller << "\t\t Prepare : " << prepareCounter.cumulated() << " s\n" );
        FLOG( FLog::Controller << "\t\t Wait : " << waitCounter.cumulated() << " s\n" );
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
        FLOG( FLog::Controller.write("\tStart Direct Pass\n").write(FLog::Flush); );
        FLOG( FTic counterTime);
        FLOG( FTic prepareCounter);
        FLOG( FTic gatherCounter);
        FLOG( FTic waitCounter);

        ///////////////////////////////////////////////////
        // Prepare data to send receive
        ///////////////////////////////////////////////////
        FLOG(prepareCounter.tic());

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

            FLOG(gatherCounter.tic());
            FMpi::MpiAssert( MPI_Allgather( partsToSend, nbProcess, MPI_INT, globalReceiveMap, nbProcess, MPI_INT, comm.getComm()),  __LINE__ );
            FLOG(gatherCounter.tac());


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
        FLOG(prepareCounter.tac());

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

        FLOG(FTic computationCounter);

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
                                           currentIter.cell->getCoordinate(), LeafIndex, AllDirs);
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
        FLOG(computationCounter.tac());
        FTRACE( regionP2PTrace.end() );

        //////////////////////////////////////////////////////////
        // Wait send receive
        //////////////////////////////////////////////////////////

        FLOG(FTic computation2Counter);

        // Create an octree with leaves from others
        OctreeClass otherP2Ptree( tree->getHeight(), tree->getSubHeight(), tree->getBoxWidth(), tree->getBoxCenter() );
        int complete = 0;
        int*const indexMessage = new int[nbProcess * 2];
        while( complete != iterRequest){
            memset(indexMessage, 0, sizeof(int) * nbProcess * 2);
            int countMessages = 0;
            // Wait data
            FLOG(waitCounter.tic());
            MPI_Waitsome(iterRequest, requests, &countMessages, indexMessage, status);
            FLOG(waitCounter.tac());
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
        FLOG( computation2Counter.tic() );

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
                const int nbNeigh = getNeighborsIndexesPeriodic(currentIter.cell->getCoordinate(), limite, indexesNeighbors, indexArray, AllDirs);

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

        FLOG(computation2Counter.tac());


        FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation L2P + P2P : " << computationCounter.elapsed() << " s\n" );
        FLOG( FLog::Controller << "\t\t Computation P2P 2 : " << computation2Counter.elapsed() << " s\n" );
        FLOG( FLog::Controller << "\t\t Prepare P2P : " << prepareCounter.elapsed() << " s\n" );
        FLOG( FLog::Controller << "\t\t Gather P2P : " << gatherCounter.elapsed() << " s\n" );
        FLOG( FLog::Controller << "\t\t Wait : " << waitCounter.elapsed() << " s\n" );

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


    long long int theoricalRepetition() const {
        if( nbLevelsAboveRoot == -1 ){
            // we know it is 3 (-1;+1)
            return 3;
        }
        // Else we find the repetition in one dir and double it
        const long long int oneDirectionRepetition = (1<<(nbLevelsAboveRoot+2)); // 2^nbLevelsAboveRoot in each dim
        const long long int fullRepetition = 2 * oneDirectionRepetition;
        return fullRepetition;
    }


    void repetitionsIntervals(FTreeCoordinate*const min, FTreeCoordinate*const max) const {
        if( nbLevelsAboveRoot == -1 ){
            // We know it is (-1;1)
            min->setPosition(-1,-1,-1);
            max->setPosition(1,1,1);
        }
        else{
            const int halfRepeated = int(theoricalRepetition()/2);
            min->setPosition(-halfRepeated,-halfRepeated,-halfRepeated);
            // if we repeat the box 8 times, we go from [-4 to 3]
            max->setPosition(halfRepeated-1,halfRepeated-1,halfRepeated-1);
        }
    }


    FReal extendedBoxWidth() const {
        // The simulation box is repeated is repeated 4 times if nbLevelsAboveRoot==-1
        // And then it doubles by two
        return tree->getBoxWidth() * FReal(1<<(nbLevelsAboveRoot+3));
    }

    /** This function has to be used to init the kernel with correct args
      * it return the box cneter seen from a kernel point of view from the periodicity the user ask for
      * this is computed using the originalBoxWidth and originalBoxCenter given in parameter
      * @param originalBoxCenter the real system center
      * @param originalBoxWidth the real system size
      * @return the center the kernel should use
      */
    FPoint extendedBoxCenter() const {
        const FReal originalBoxWidth     = tree->getBoxWidth();
        const FReal originalBoxWidthDiv2 = originalBoxWidth/2.0;
        const FPoint originalBoxCenter   = tree->getBoxCenter();

        const FReal offset = extendedBoxWidth()/FReal(2.0);
        return FPoint( originalBoxCenter.getX() - originalBoxWidthDiv2 + offset,
                       originalBoxCenter.getY() - originalBoxWidthDiv2 + offset,
                       originalBoxCenter.getZ() - originalBoxWidthDiv2 + offset);
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
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FLOG( FLog::Controller.write("\tStart Periodic Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);

        if( nbLevelsAboveRoot != -1 ){
            // we will use offsetRealTree-1 cells but for simplicity allocate offsetRealTree
            // upperCells[offsetRealTree-1] is root cell
            CellClass*const upperCells = new CellClass[offsetRealTree];
            {
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoLeft();
                kernels[0]->M2M( &upperCells[offsetRealTree-1], octreeIterator.getCurrentBox(), offsetRealTree);
            }
            {
                CellClass* virtualChild[8];
                for(int idxLevel = offsetRealTree-1 ; idxLevel > 1  ; --idxLevel){
                    FMemUtils::setall(virtualChild,&upperCells[idxLevel],8);
                    kernels[0]->M2M( &upperCells[idxLevel-1], virtualChild, idxLevel);
                }
            }
            CellClass*const downerCells = new CellClass[offsetRealTree];

            {
                const int idxUpperLevel = 2;

                const CellClass* neighbors[343];
                memset(neighbors, 0, sizeof(CellClass*) * 343);
                int counter = 0;
                for(int idxX = -2 ; idxX <= 1 ; ++idxX){
                    for(int idxY = -2 ; idxY <= 1 ; ++idxY){
                        for(int idxZ = -2 ; idxZ <= 1 ; ++idxZ){
                            if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                                neighbors[neighIndex(idxX,idxY,idxZ)] = &upperCells[idxUpperLevel-1];
                                ++counter;
                            }
                        }
                    }
                }
                // compute M2L
                kernels[0]->M2L( &downerCells[idxUpperLevel-1] , neighbors, counter, idxUpperLevel);
            }

            for(int idxUpperLevel = 3 ; idxUpperLevel <= offsetRealTree ; ++idxUpperLevel){
                const CellClass* neighbors[343];
                memset(neighbors, 0, sizeof(CellClass*) * 343);
                int counter = 0;
                for(int idxX = -2 ; idxX <= 3 ; ++idxX){
                    for(int idxY = -2 ; idxY <= 3 ; ++idxY){
                        for(int idxZ = -2 ; idxZ <= 3 ; ++idxZ){
                            if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                                neighbors[neighIndex(idxX,idxY,idxZ)] = &upperCells[idxUpperLevel-1];
                                ++counter;
                            }
                        }
                    }
                }

                // compute M2L
                kernels[0]->M2L( &downerCells[idxUpperLevel-1] , neighbors, counter, idxUpperLevel);
            }

            {
                CellClass* virtualChild[8];
                memset(virtualChild, 0, sizeof(CellClass*) * 8);
                for(int idxLevel = 2 ; idxLevel <= offsetRealTree-1  ; ++idxLevel){
                    virtualChild[0] = &downerCells[idxLevel];
                    kernels[0]->L2L( &downerCells[idxLevel-1], virtualChild, idxLevel);
                }
            }

            // L2L from 0 to level 1
            {
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoLeft();
                kernels[0]->L2L( &downerCells[offsetRealTree-1], octreeIterator.getCurrentBox(), offsetRealTree);
            }

            delete[] upperCells;
            delete[] downerCells;
        }

        FLOG( FLog::Controller << "\tFinished (@Periodic = "  << counterTime.tacAndElapsed() << "s)\n" );
    }


};






#endif //FFMMALGORITHMTHREAD_HPP


