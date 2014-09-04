// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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
#ifndef FFMMALGORITHMTHREADPROC_HPP
#define FFMMALGORITHMTHREADPROC_HPP

#include <omp.h>

//
#include "../Utils/FAssert.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"

#include "../Utils/FGlobal.hpp"

#include "../Containers/FBoolArray.hpp"
#include "../Containers/FOctree.hpp"
#include "../Containers/FLightOctree.hpp"

#include "../Containers/FBufferWriter.hpp"
#include "../Containers/FBufferReader.hpp"
#include "../Containers/FMpiBufferWriter.hpp"
#include "../Containers/FMpiBufferReader.hpp"

#include "../Utils/FMpi.hpp"
#include <sys/time.h>

#include "FCoreCommon.hpp"

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FFmmAlgorithmThreadProc
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
class FFmmAlgorithmThreadProc : public FAbstractAlgorithm {
    // Can be deleted
    // const static int MaxSizePerCell = CellClass::GetSize();

    OctreeClass* const tree;                 //< The octree to work on
    KernelClass** kernels;                   //< The kernels

    const FMpi::FComm& comm;                 //< MPI comm

    typename OctreeClass::Iterator*     iterArray;  //Will be used to store pointers to cells/leafs to work with
    typename OctreeClass::Iterator* iterArrayComm;  //Will be used to store pointers to cells/leafs to send/rcv
    int numberOfLeafs;                          //< To store the size at the previous level

    const int MaxThreads;               //< the max number of thread allowed by openmp

    const int nbProcess;                //< Number of process
    const int idProcess;                //< Id of current process

    const int OctreeHeight;            //<Height of the tree


    /** An interval is the morton index interval
     * that a proc use (it holds data in this interval)
     */
    struct Interval{
	MortonIndex leftIndex;
	MortonIndex rightIndex;
    };
    /** My interval */
    Interval*const intervals;
    /** All process intervals */
    Interval*const workingIntervalsPerLevel;

    /** Get an interval from proc id and level */
    Interval& getWorkingInterval( int level,  int proc){
	return workingIntervalsPerLevel[OctreeHeight * proc + level];
    }

    const Interval& getWorkingInterval( int level,  int proc) const {
	return workingIntervalsPerLevel[OctreeHeight * proc + level];
    }

    /** To know if a proc has work at a given level (if it hold cells and was responsible of them) */
    bool procHasWorkAtLevel(const int idxLevel , const int idxProc) const {
	return getWorkingInterval(idxLevel, idxProc).leftIndex <= getWorkingInterval(idxLevel, idxProc).rightIndex;
    }

    /** Return true if the idxProc left cell at idxLevel+1 has the same parent as us for our right cell */
    bool procCoversMyRightBorderCell(const int idxLevel , const int idxProc) const {
	return (getWorkingInterval((idxLevel+1) , idProcess).rightIndex>>3) == (getWorkingInterval((idxLevel+1) ,idxProc).leftIndex >>3);
    }

    /** Return true if the idxProc right cell at idxLevel+1 has the same parent as us for our left cell */
    bool procCoversMyLeftBorderCell(const int idxLevel , const int idxProc) const {
	return (getWorkingInterval((idxLevel+1) , idxProc).rightIndex >>3) == (getWorkingInterval((idxLevel+1) , idProcess).leftIndex>>3);
    }

public:
    /** Get current proc interval at level */
    Interval& getWorkingInterval( int level){
	return getWorkingInterval(level, idProcess);
    }

    /** Does the current proc has some work at this level */
    bool hasWorkAtLevel( int level){
	return idProcess == 0 || (getWorkingInterval(level, idProcess - 1).rightIndex) < (getWorkingInterval(level, idProcess).rightIndex);
    }

    /** The constructor need the octree and the kernels used for computation
     * @param inTree the octree to work on
     * @param inKernels the kernels to call
     * An assert is launched if one of the arguments is null
     */
    FFmmAlgorithmThreadProc(const FMpi::FComm& inComm, OctreeClass* const inTree, KernelClass* const inKernels)
	: tree(inTree) , kernels(nullptr), comm(inComm), iterArray(nullptr),iterArrayComm(nullptr),numberOfLeafs(0),
	  MaxThreads(omp_get_max_threads()), nbProcess(inComm.processCount()), idProcess(inComm.processId()),
	  OctreeHeight(tree->getHeight()),intervals(new Interval[inComm.processCount()]),
	  workingIntervalsPerLevel(new Interval[inComm.processCount() * tree->getHeight()])
    {
	FAssertLF(tree, "tree cannot be null");


	this->kernels = new KernelClass*[MaxThreads];
	for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
	    this->kernels[idxThread] = new KernelClass(*inKernels);
	}

	FLOG(FLog::Controller << "FFmmAlgorithmThreadProc\n");
	FLOG(FLog::Controller << "Max threads = "  << MaxThreads << ", Procs = " << nbProcess << ", I am " << idProcess << ".\n");
    }
    /** Default destructor */
    virtual ~FFmmAlgorithmThreadProc(){
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

	    Interval myFullInterval;
	    {//Building the interval with the first and last leaves (and count the number of leaves)
		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.gotoBottomLeft();
		myFullInterval.leftIndex = octreeIterator.getCurrentGlobalIndex();
		do{
		    ++this->numberOfLeafs;
		} while(octreeIterator.moveRight());
		myFullInterval.rightIndex = octreeIterator.getCurrentGlobalIndex();
	    }
	    // Allocate a number to store the pointer of the cells at a level
	    iterArray     = new typename OctreeClass::Iterator[numberOfLeafs];
	    iterArrayComm = new typename OctreeClass::Iterator[numberOfLeafs];
	    FAssertLF(iterArray,     "iterArray     bad alloc");
	    FAssertLF(iterArrayComm, "iterArrayComm bad alloc");

	    // We get the leftIndex/rightIndex indexes from each procs
	    FMpi::MpiAssert( MPI_Allgather( &myFullInterval, sizeof(Interval), MPI_BYTE, intervals, sizeof(Interval), MPI_BYTE, comm.getComm()),  __LINE__ );

	    // Build my intervals for all levels
	    std::unique_ptr<Interval[]> myIntervals(new Interval[OctreeHeight]);
	    // At leaf level we know it is the full interval
	    myIntervals[OctreeHeight - 1] = myFullInterval;

	    // We can estimate the interval for each level by using the parent/child relation
	    for(int idxLevel = OctreeHeight - 2 ; idxLevel >= 0 ; --idxLevel){
		myIntervals[idxLevel].leftIndex = myIntervals[idxLevel+1].leftIndex >> 3;
		myIntervals[idxLevel].rightIndex = myIntervals[idxLevel+1].rightIndex >> 3;
	    }

	    // Process 0 uses the estimates as real intervals, but other processes
	    // should remove cells that belong to others
	    if(idProcess != 0){
		//We test for each level if process on left (idProcess-1) own cell I thought I owned
		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.gotoBottomLeft();
		octreeIterator.moveUp();

		// At h-1 the working limit is the parent of the right cell of the proc on the left
		MortonIndex workingLimitAtLevel = intervals[idProcess-1].rightIndex >> 3;

		// We check if we have no more work to do
		int nullIntervalFromLevel = 0;

		for(int idxLevel = OctreeHeight - 2 ; idxLevel >= 1 && nullIntervalFromLevel == 0 ; --idxLevel){
		    while(octreeIterator.getCurrentGlobalIndex() <= workingLimitAtLevel){
			if( !octreeIterator.moveRight() ){
			    // We cannot move right we are not owner of any more cell
			    nullIntervalFromLevel = idxLevel;
			    break;
			}
		    }
		    // If we are responsible for some cells at this level keep the first index
		    if(nullIntervalFromLevel == 0){
			myIntervals[idxLevel].leftIndex = octreeIterator.getCurrentGlobalIndex();
			octreeIterator.moveUp();
			workingLimitAtLevel >>= 3;
		    }
		}
		// In case we are not responsible for any cells we put the leftIndex = rightIndex+1
		for(int idxLevel = nullIntervalFromLevel ; idxLevel >= 1 ; --idxLevel){
		    myIntervals[idxLevel].leftIndex = myIntervals[idxLevel].rightIndex + 1;
		}
	    }

	    // We get the leftIndex/rightIndex indexes from each procs
	    FMpi::MpiAssert( MPI_Allgather( myIntervals.get(), int(sizeof(Interval)) * OctreeHeight, MPI_BYTE,
					    workingIntervalsPerLevel, int(sizeof(Interval)) * OctreeHeight, MPI_BYTE, comm.getComm()),  __LINE__ );
	}

	// run;
	if(operationsToProceed & FFmmP2M) bottomPass();

	if(operationsToProceed & FFmmM2M) upwardPass();

	if(operationsToProceed & FFmmM2L) transferPass();

	if(operationsToProceed & FFmmL2L) downardPass();

	if((operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P)) directPassOld();


	// delete array
	delete []     iterArray;
	delete [] iterArrayComm;
	iterArray     = nullptr;
	iterArrayComm = nullptr;
    }

private:

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /**
     * P2M Bottom Pass
     * No communication are involved in the P2M.
     * It is similar to multi threaded version.
     */
    void bottomPass(){
	FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
	FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
	FLOG(FTic counterTime);
	FLOG(FTic computationCounter);
	typename OctreeClass::Iterator octreeIterator(tree);

	// Copy the ptr to leaves in array
	octreeIterator.gotoBottomLeft();
	int leafs = 0;
	do{
	    iterArray[leafs++] = octreeIterator;
	} while(octreeIterator.moveRight());

	FLOG(computationCounter.tic());
	#pragma omp parallel
	{
	    // Each thread get its own kernel
	    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
	    // Parallel iteration on the leaves
	    #pragma omp for nowait
	    for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
		myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentListSrc());
	    }
	}
	FLOG(computationCounter.tac());
	FLOG( FLog::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << " s)\n" );
	FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
	const int MaxSizePerCell = CellClass::GetSize();
	FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
	FLOG( FLog::Controller.write("\tStart Upward Pass\n").write(FLog::Flush); );
	FLOG(FTic counterTime);
	FLOG(FTic computationCounter);
	FLOG(FTic singleCounter);
	FLOG(FTic parallelCounter);

	// Start from leal level (height-1)
	typename OctreeClass::Iterator octreeIterator(tree);
	octreeIterator.gotoBottomLeft();
	octreeIterator.moveUp();
	typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

	// The proc to send the shared cells to
	// Starting to the proc on the left this variable will go to 0
	int currentProcIdToSendTo = (idProcess - 1);

	// There are a maximum of 1 sends and 8-1 receptions
	MPI_Request requests[8];
	MPI_Status status[8];

	// Maximum data per message is:
	FMpiBufferWriter sendBuffer(comm.getComm(), 7*MaxSizePerCell + 1);
	const int recvBufferOffset = (7 * MaxSizePerCell + 1);
	FMpiBufferReader recvBuffer(comm.getComm(), 7*recvBufferOffset);
	CellClass recvBufferCells[8];

	// The first proc that send to me a cell
	// This variable will go to nbProcess
	int firstProcThatSend = idProcess + 1;
	FLOG(computationCounter.tic());

	// We work from height-1 to 1
	for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
	    // Does my cells are covered by my neighbors working interval and so I have no more work?
	    const bool noMoreWorkForMe = (idProcess != 0 && !procHasWorkAtLevel(idxLevel+1, idProcess));
	    if(noMoreWorkForMe){
		FAssertLF(procHasWorkAtLevel(idxLevel, idProcess) == false);
		break;
	    }

	    // Copy and count ALL the cells (even the ones outside the working interval)
	    int totalNbCellsAtLevel = 0;
	    do{
		iterArray[totalNbCellsAtLevel++] = octreeIterator;
	    } while(octreeIterator.moveRight());
	    avoidGotoLeftIterator.moveUp();
	    octreeIterator = avoidGotoLeftIterator;

	    int iterMpiRequests   = 0; // The iterator for send/recv requests

	    int nbCellsToSkip     = 0; // The number of cells to send
	    // Skip all the cells that are out of my working interval
	    while(nbCellsToSkip < totalNbCellsAtLevel && iterArray[nbCellsToSkip].getCurrentGlobalIndex() < getWorkingInterval(idxLevel, idProcess).leftIndex){
		++nbCellsToSkip;
	    }

	    // We need to know if we will recv something in order to know if threads skip the last cell
	    int nbCellsForThreads = totalNbCellsAtLevel; // totalNbCellsAtLevel or totalNbCellsAtLevel-1
	    bool hasToReceive = false;
	    if(idProcess != nbProcess-1 && procHasWorkAtLevel(idxLevel , idProcess)){
		// Find the first proc that may send to me
		while(firstProcThatSend < nbProcess && !procHasWorkAtLevel(idxLevel+1, firstProcThatSend) ){
		    firstProcThatSend += 1;
		}
		// Do we have to receive?
		if(firstProcThatSend < nbProcess && procHasWorkAtLevel(idxLevel+1, firstProcThatSend) && procCoversMyRightBorderCell(idxLevel, firstProcThatSend) ){
		    hasToReceive = true;
		    // Threads do not compute the last cell, we will do it once data are received
		    nbCellsForThreads -= 1;
		}
	    }

	    FLOG(parallelCounter.tic());
	    #pragma omp parallel
	    {
		const int threadNumber = omp_get_thread_num();
		KernelClass* myThreadkernels = (kernels[threadNumber]);
		//This single section post and receive the comms, and then do the M2M associated with it.
		#pragma omp single nowait
		{
		    FLOG(singleCounter.tic());
		    // Master proc never send
		    if(idProcess != 0){
			// Skip process that have no work at that level
			while( currentProcIdToSendTo && !procHasWorkAtLevel(idxLevel, currentProcIdToSendTo)  ){
			    --currentProcIdToSendTo;
			}
			// Does the next proc that has work is sharing the parent of my left cell
			if(procHasWorkAtLevel(idxLevel, currentProcIdToSendTo) && procCoversMyLeftBorderCell(idxLevel, currentProcIdToSendTo)){
			    FAssertLF(nbCellsToSkip != 0);

			    char packageFlags = 0;
			    sendBuffer.write(packageFlags);

			    // Only the cell the most on the right out of my working interval should be taken in
			    // consideration (at pos nbCellsToSkip-1) other (x < nbCellsToSkip-1) have already been sent
			    const CellClass* const* const child = iterArray[nbCellsToSkip-1].getCurrentChild();
			    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
				// Check if child exists and it was part of my working interval
				if( child[idxChild] && getWorkingInterval((idxLevel+1), idProcess).leftIndex <= child[idxChild]->getMortonIndex() ){
				    // Add the cell to the buffer
				    child[idxChild]->serializeUp(sendBuffer);
				    packageFlags = char(packageFlags | (0x1 << idxChild));
				}
			    }
			    // Add the flag as first value
			    sendBuffer.writeAt(0,packageFlags);
			    // Post the message
			    MPI_Isend(sendBuffer.data(), sendBuffer.getSize(), MPI_PACKED, currentProcIdToSendTo,
				      FMpi::TagFmmM2M + idxLevel, comm.getComm(), &requests[iterMpiRequests++]);
			}
		    }

		    //Post receive, Datas needed in several parts of the section
		    int nbProcThatSendToMe = 0;

		    if(hasToReceive){
			//Test : if the firstProcThatSend father minimal value in interval is lesser than mine
			int idProcSource = firstProcThatSend;
			// Find the last proc that should send to me
			while( idProcSource < nbProcess
			       && ( !procHasWorkAtLevel(idxLevel+1, idProcSource) || procCoversMyRightBorderCell(idxLevel, idProcSource) )){
			    if(procHasWorkAtLevel(idxLevel+1, idProcSource) && procCoversMyRightBorderCell(idxLevel, idProcSource)){
				MPI_Irecv(&recvBuffer.data()[nbProcThatSendToMe * recvBufferOffset], recvBufferOffset, MPI_PACKED,
					idProcSource, FMpi::TagFmmM2M + idxLevel, comm.getComm(), &requests[iterMpiRequests++]);
				nbProcThatSendToMe += 1;
				FAssertLF(nbProcThatSendToMe <= 7);
			    }
			    ++idProcSource;
			}
		    }

		    //Wait For the comms, and do the work
		    // Are we sending or waiting anything?
		    if(iterMpiRequests){
			FAssertLF(iterMpiRequests <= 8);
			MPI_Waitall( iterMpiRequests, requests, status);
		    }
		    // We had received something so we need to proceed the last M2M
		    if( hasToReceive ){
			FAssertLF(iterMpiRequests != 0);
			CellClass* currentChild[8];
			memcpy(currentChild, iterArray[totalNbCellsAtLevel - 1].getCurrentChild(), 8 * sizeof(CellClass*));

			// Retreive data and merge my child and the child from others
			for(int idxProc = 0 ; idxProc < nbProcThatSendToMe ; ++idxProc){
			    recvBuffer.seek(idxProc * recvBufferOffset);
			    int packageFlags = int(recvBuffer.getValue<char>());

			    int position = 0;
			    while( packageFlags && position < 8){
				while(!(packageFlags & 0x1)){
				    packageFlags >>= 1;
				    ++position;
				}
				FAssertLF(!currentChild[position], "Already has a cell here");
				recvBufferCells[position].deserializeUp(recvBuffer);
				currentChild[position] = (CellClass*) &recvBufferCells[position];

				packageFlags >>= 1;
				++position;
			    }
			}
			// Finally compute
			(*kernels[threadNumber]).M2M( iterArray[totalNbCellsAtLevel - 1].getCurrentCell() , currentChild, idxLevel);
			firstProcThatSend += nbProcThatSendToMe - 1;
		    }
		    // Reset buffer
		    sendBuffer.reset();
		    recvBuffer.seek(0);
		    FLOG(singleCounter.tac());
		}//End Of Single section

		// All threads proceed the M2M
		#pragma omp for nowait
		for( int idxCell = nbCellsToSkip ; idxCell < nbCellsForThreads ; ++idxCell){
		    myThreadkernels->M2M( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
		}
	    }//End of parallel section
	    FLOG(parallelCounter.tac());
	}

	FLOG(counterTime.tac());
	FLOG(computationCounter.tac());
	FLOG( FLog::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.elapsed() << " s)\n" );
	FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );
	FLOG( FLog::Controller << "\t\t Wait : " << singleCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Wait : " << parallelCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downard
    /////////////////////////////////////////////////////////////////////////////


    void transferPass(){
	const int MaxSizePerCell = CellClass::GetSize();
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

	// All process say to each others
	// what the will send to who
	int*const globalReceiveMap = new int[nbProcess * nbProcess * OctreeHeight];
	memset(globalReceiveMap, 0, sizeof(int) * nbProcess * nbProcess * OctreeHeight);

	FMpiBufferWriter**const sendBuffer = new FMpiBufferWriter*[nbProcess * OctreeHeight];
	memset(sendBuffer, 0, sizeof(FMpiBufferWriter*) * nbProcess * OctreeHeight);

	FMpiBufferReader**const recvBuffer = new FMpiBufferReader*[nbProcess * OctreeHeight];
	memset(recvBuffer, 0, sizeof(FMpiBufferReader*) * nbProcess * OctreeHeight);

	#pragma omp parallel
	{
	    #pragma omp master
	    {
		{
		    FTRACE( FTrace::FRegion regionTrace( "Preprocess" , __FUNCTION__ , __FILE__ , __LINE__) );
		    FLOG(prepareCounter.tic());

		    std::unique_ptr<typename OctreeClass::Iterator[]> iterArrayLocal(new typename OctreeClass::Iterator[numberOfLeafs]);

		    // To know if a leaf has been already sent to a proc
		    bool*const alreadySent = new bool[nbProcess];
		    memset(alreadySent, 0, sizeof(bool) * nbProcess);

		    typename OctreeClass::Iterator octreeIterator(tree);
		    octreeIterator.moveDown();
		    typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
		    // for each levels
		    for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
			if(!procHasWorkAtLevel(idxLevel, idProcess)){
			    avoidGotoLeftIterator.moveDown();
			    octreeIterator = avoidGotoLeftIterator;
			    continue;
			}

			int numberOfCells = 0;

			while(octreeIterator.getCurrentGlobalIndex() <  getWorkingInterval(idxLevel , idProcess).leftIndex){
			    octreeIterator.moveRight();
			}

			// for each cells
			do{
			    iterArrayLocal[numberOfCells] = octreeIterator;
			    ++numberOfCells;
			} while(octreeIterator.moveRight());
			avoidGotoLeftIterator.moveDown();
			octreeIterator = avoidGotoLeftIterator;

			leafsNeedOther[idxLevel] = new FBoolArray(numberOfCells);

			// Which cell potentialy needs other data and in the same time
			// are potentialy needed by other
			MortonIndex neighborsIndexes[189];
			for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
			    // Find the M2L neigbors of a cell
			    const int counter = iterArrayLocal[idxCell].getCurrentGlobalCoordinate().getInteractionNeighbors(idxLevel, neighborsIndexes);

			    memset(alreadySent, false, sizeof(bool) * nbProcess);
			    bool needOther = false;
			    // Test each negibors to know which one do not belong to us
			    for(int idxNeigh = 0 ; idxNeigh < counter ; ++idxNeigh){
				if(neighborsIndexes[idxNeigh] < getWorkingInterval(idxLevel , idProcess).leftIndex
					|| (getWorkingInterval(idxLevel , idProcess).rightIndex) < neighborsIndexes[idxNeigh]){
				    int procToReceive = idProcess;
				    while( 0 != procToReceive && neighborsIndexes[idxNeigh] < getWorkingInterval(idxLevel , procToReceive).leftIndex ){
					--procToReceive;
				    }
				    while( procToReceive != nbProcess -1 && (getWorkingInterval(idxLevel , procToReceive).rightIndex) < neighborsIndexes[idxNeigh]){
					++procToReceive;
				    }
				    // Maybe already sent to that proc?
				    if( !alreadySent[procToReceive]
					    && getWorkingInterval(idxLevel , procToReceive).leftIndex <= neighborsIndexes[idxNeigh]
					    && neighborsIndexes[idxNeigh] <= getWorkingInterval(idxLevel , procToReceive).rightIndex){

					alreadySent[procToReceive] = true;

					needOther = true;

					toSend[idxLevel * nbProcess + procToReceive].push(iterArrayLocal[idxCell]);
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

		for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
		    for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
			const int toSendAtProcAtLevel = indexToSend[idxLevel * nbProcess + idxProc];
			if(toSendAtProcAtLevel != 0){
			    sendBuffer[idxLevel * nbProcess + idxProc] = new FMpiBufferWriter(comm.getComm(),toSendAtProcAtLevel * SizeOfCellToSend);

			    for(int idxLeaf = 0 ; idxLeaf < toSendAtProcAtLevel; ++idxLeaf){
				const MortonIndex cellIndex = toSend[idxLevel * nbProcess + idxProc][idxLeaf].getCurrentGlobalIndex();
				sendBuffer[idxLevel * nbProcess + idxProc]->write(cellIndex);
				toSend[idxLevel * nbProcess + idxProc][idxLeaf].getCurrentCell()->serializeUp(*sendBuffer[idxLevel * nbProcess + idxProc]);
			    }

			    FMpi::MpiAssert( MPI_Isend( sendBuffer[idxLevel * nbProcess + idxProc]->data(),
					     sendBuffer[idxLevel * nbProcess + idxProc]->getSize(),MPI_PACKED, idxProc,
				    FMpi::TagLast + idxLevel, comm.getComm(), &requests[iterRequest++]) , __LINE__ );
			}

			const int toReceiveFromProcAtLevel = globalReceiveMap[(idxProc * nbProcess * OctreeHeight) + idxLevel * nbProcess + idProcess];
			if(toReceiveFromProcAtLevel){
			    recvBuffer[idxLevel * nbProcess + idxProc] = new FMpiBufferReader(comm.getComm(),toReceiveFromProcAtLevel * SizeOfCellToSend);

			    FMpi::MpiAssert( MPI_Irecv(recvBuffer[idxLevel * nbProcess + idxProc]->data(),
					     recvBuffer[idxLevel * nbProcess + idxProc]->getCapacity(), MPI_PACKED,idxProc,
				    FMpi::TagLast + idxLevel, comm.getComm(), &requests[iterRequest++]) , __LINE__ );
			}
		    }
		}

		//////////////////////////////////////////////////////////////////
		// Wait received data and compute
		//////////////////////////////////////////////////////////////////

		// Wait to receive every things (and send every things)
		MPI_Waitall(iterRequest, requests, status);

		delete[] requests;
		delete[] status;

		FLOG(sendCounter.tac());
	    }

	    //////////////////////////////////////////////////////////////////
	    // Do M2L
	    //////////////////////////////////////////////////////////////////

	    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
	    const CellClass* neighbors[343];

	    #pragma omp single nowait
	    {
		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.moveDown();
		typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
		// Now we can compute all the data
		// for each levels
		for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
		    if(!procHasWorkAtLevel(idxLevel, idProcess)){
			avoidGotoLeftIterator.moveDown();
			octreeIterator = avoidGotoLeftIterator;
			continue;
		    }

		    int numberOfCells = 0;
		    while(octreeIterator.getCurrentGlobalIndex() <  getWorkingInterval(idxLevel , idProcess).leftIndex){
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
		    {
			const int chunckSize = FMath::Max(1, numberOfCells/(omp_get_num_threads()*omp_get_num_threads()));
			for(int idxCell = 0 ; idxCell < numberOfCells ; idxCell += chunckSize){
			    #pragma omp task
			    {
				const int nbCellToCompute = FMath::Min(chunckSize, numberOfCells-idxCell);
				for(int idxCellToCompute = idxCell ; idxCellToCompute < idxCell+nbCellToCompute ; ++idxCellToCompute){
				    const int counter = tree->getInteractionNeighbors(neighbors,  iterArray[idxCellToCompute].getCurrentGlobalCoordinate(), idxLevel);
				    if(counter) myThreadkernels->M2L( iterArray[idxCellToCompute].getCurrentCell() , neighbors, counter, idxLevel);
				}
			    }
			}
		    }

		    #pragma omp taskwait

		    for(int idxThread = 0 ; idxThread < omp_get_num_threads() ; ++idxThread){
			#pragma omp task
			{
			    kernels[idxThread]->finishedLevelM2L(idxLevel);
			}
		    }

		    FLOG(computationCounter.tac());
		}
	    }
	}


	{
	    FTRACE( FTrace::FRegion regionTrace("Compute Received data", __FUNCTION__ , __FILE__ , __LINE__) );
	    FLOG(receiveCounter.tic());
	    typename OctreeClass::Iterator octreeIterator(tree);
	    octreeIterator.moveDown();
	    typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
	    // compute the second time
	    // for each levels
	    for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
		if(!procHasWorkAtLevel(idxLevel, idProcess)){
		    avoidGotoLeftIterator.moveDown();
		    octreeIterator = avoidGotoLeftIterator;
		    continue;
		}

		// put the received data into a temporary tree
		FLightOctree<CellClass> tempTree;
		for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
		    const int toReceiveFromProcAtLevel = globalReceiveMap[(idxProc * nbProcess * OctreeHeight) + idxLevel * nbProcess + idProcess];

		    for(int idxCell = 0 ; idxCell < toReceiveFromProcAtLevel ; ++idxCell){
			const MortonIndex cellIndex = recvBuffer[idxLevel * nbProcess + idxProc]->FMpiBufferReader::getValue<MortonIndex>();

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

		while(octreeIterator.getCurrentGlobalIndex() <  getWorkingInterval(idxLevel , idProcess).leftIndex){
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
		leafsNeedOther[idxLevel] = nullptr;

		// Compute this cells
		FLOG(computationCounter.tic());
		#pragma omp parallel
		{
		    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
		    MortonIndex neighborsIndex[189];
		    int neighborsPosition[189];
		    const CellClass* neighbors[343];

		    #pragma omp for schedule(static) nowait
		    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
			// compute indexes
			memset(neighbors, 0, 343 * sizeof(CellClass*));
			const int counterNeighbors = iterArray[idxCell].getCurrentGlobalCoordinate().getInteractionNeighbors(idxLevel, neighborsIndex, neighborsPosition);

			int counter = 0;
			// does we receive this index from someone?
			for(int idxNeig = 0 ;idxNeig < counterNeighbors ; ++idxNeig){
			    if(neighborsIndex[idxNeig] < (getWorkingInterval(idxLevel , idProcess).leftIndex)
				    || (getWorkingInterval(idxLevel , idProcess).rightIndex) < neighborsIndex[idxNeig]){

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
			    myThreadkernels->M2L( iterArray[idxCell].getCurrentCell() , neighbors, counter, idxLevel);
			}
		    }

		    myThreadkernels->finishedLevelM2L(idxLevel);
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


	FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << " s)\n" );
	FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Gather : " << gatherCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Prepare : " << prepareCounter.cumulated() << " s\n" );

    }

    //////////////////////////////////////////////////////////////////
    // ---------------- L2L ---------------
    //////////////////////////////////////////////////////////////////

    void downardPass(){ // second L2L
	const int MaxSizePerCell = CellClass::GetSize();
	FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
	FLOG( FLog::Controller.write("\tStart Downward Pass (L2L)\n").write(FLog::Flush); );
	FLOG(FTic counterTime);
	FLOG(FTic computationCounter);
	FLOG(FTic prepareCounter);
	FLOG(FTic waitCounter);

	// Start from leal level - 1
	typename OctreeClass::Iterator octreeIterator(tree);
	octreeIterator.moveDown();
	typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

	// Max 1 receive and 7 send (but 7 times the same data)
	MPI_Request*const requests = new MPI_Request[8];
	MPI_Status*const status = new MPI_Status[8];

	const int heightMinusOne = OctreeHeight - 1;

	FMpiBufferWriter sendBuffer(comm.getComm(),MaxSizePerCell);
	FMpiBufferReader recvBuffer(comm.getComm(),MaxSizePerCell);

	int righestProcToSendTo   = nbProcess - 1;

	// for each levels exepted leaf level
	for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
	    // If nothing to do in the next level skip the current one
	    if(idProcess != 0 && !procHasWorkAtLevel(idxLevel+1, idProcess) ){
		avoidGotoLeftIterator.moveDown();
		octreeIterator = avoidGotoLeftIterator;
		continue;
	    }

	    // Copy all the cells in an array even the one that are out of my working interval
	    int totalNbCellsAtLevel = 0;
	    do{
		iterArray[totalNbCellsAtLevel++] = octreeIterator;
	    } while(octreeIterator.moveRight());
	    avoidGotoLeftIterator.moveDown();
	    octreeIterator = avoidGotoLeftIterator;

	    // Count the number of cells that are out of my working interval
	    int nbCellsToSkip = 0;
	    while(nbCellsToSkip < totalNbCellsAtLevel && iterArray[nbCellsToSkip].getCurrentGlobalIndex() < getWorkingInterval(idxLevel , idProcess).leftIndex){
		nbCellsToSkip += 1;
	    }

	    // Check if someone will send a cell to me
	    bool hasToReceive = false;
	    int idxProcToReceive = idProcess - 1;
	    if(idProcess != 0 && nbCellsToSkip){
		// Starting from my left neighbor stop at the first proc that has work to do (not null interval)
		while(idxProcToReceive && !procHasWorkAtLevel(idxLevel, idxProcToReceive) ){
		    idxProcToReceive -= 1;
		}
		// Check if we find such a proc and that it share a cell with us on the border
		if(procHasWorkAtLevel(idxLevel, idxProcToReceive) && procCoversMyLeftBorderCell(idxLevel, idxProcToReceive)){
		    hasToReceive = true;
		}
	    }

	    #pragma omp parallel
	    {
		int threadNumber = omp_get_thread_num();
		KernelClass* myThreadkernels = (kernels[threadNumber]);
		#pragma omp single nowait
		{
		    FLOG(prepareCounter.tic());
		    int iterRequests = 0;
		    // Post the receive
		    if(hasToReceive){
			FMpi::MpiAssert( MPI_Irecv( recvBuffer.data(), recvBuffer.getCapacity(), MPI_PACKED, idxProcToReceive,
				   FMpi::TagFmmL2L + idxLevel, comm.getComm(), &requests[iterRequests++]), __LINE__ );
		    }

		    // We have to be sure that we are not sending if we have no work in the current level
		    if(idProcess != nbProcess - 1 && idProcess < righestProcToSendTo && procHasWorkAtLevel(idxLevel, idProcess)){
			int idxProcSend = idProcess + 1;
			int nbMessageSent = 0;
			// From the proc on the right to righestProcToSendTo, check if we have to send something
			while(idxProcSend <= righestProcToSendTo && ( !procHasWorkAtLevel(idxLevel+1, idxProcSend) || procCoversMyRightBorderCell(idxLevel, idxProcSend)) ){
			    // We know that if the proc has work at the next level it share a cell with us due to the while condition
			    if(procHasWorkAtLevel(idxLevel+1, idxProcSend)){
				FAssertLF(procCoversMyRightBorderCell(idxLevel, idxProcSend));
				// If first message then serialize the cell to send
				if( nbMessageSent == 0 ){
				    // We send our last cell
				    iterArray[totalNbCellsAtLevel - 1].getCurrentCell()->serializeDown(sendBuffer);
				}
				// Post the send message
				FMpi::MpiAssert( MPI_Isend(sendBuffer.data(), sendBuffer.getSize(), MPI_PACKED, idxProcSend,
					  FMpi::TagFmmL2L + idxLevel, comm.getComm(), &requests[iterRequests++]), __LINE__);
				// Inc and check the counter
				nbMessageSent += 1;
				FAssertLF(nbMessageSent <= 7);
			    }
			    idxProcSend += 1;
			}
			// Next time we will not need to go further than idxProcSend
			righestProcToSendTo = idxProcSend;
		    }
		    // Finalize the communication
		    if(iterRequests){
			FLOG(waitCounter.tic());
			FAssertLF(iterRequests <= 8);
			FMpi::MpiAssert(MPI_Waitall( iterRequests, requests, status), __LINE__);
			FLOG(waitCounter.tac());
		    }
		    // If we receive something proceed the L2L
		    if(hasToReceive){
			FAssertLF(iterRequests != 0);
			// In this case we know that we have to perform the L2L with the last cell that are
			// exclude from our working interval nbCellsToSkip-1
			iterArray[nbCellsToSkip-1].getCurrentCell()->deserializeDown(recvBuffer);
			kernels[threadNumber]->L2L( iterArray[nbCellsToSkip-1].getCurrentCell() , iterArray[nbCellsToSkip-1].getCurrentChild(), idxLevel);
		    }
		    FLOG(prepareCounter.tac());
		}

		#pragma omp single nowait
		{
		    FLOG(computationCounter.tic());
		}
		// Threads are working on all the cell of our working interval at that level
		#pragma omp for nowait
		for(int idxCell = nbCellsToSkip ; idxCell < totalNbCellsAtLevel ; ++idxCell){
		    myThreadkernels->L2L( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
		}
	    }
	    FLOG(computationCounter.tac());

	    sendBuffer.reset();
	    recvBuffer.seek(0);
	}

	delete[] requests;
	delete[] status;

	FLOG( FLog::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << " s)\n" );
	FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Prepare : " << prepareCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Wait : " << waitCounter.cumulated() << " s\n" );
    }


    /////////////////////////////////////////////////////////////////////////////
    // Direct
    /////////////////////////////////////////////////////////////////////////////
    struct LeafData{
	FTreeCoordinate coord;
	CellClass* cell;
	ContainerClass* targets;
	ContainerClass* sources;
    };


    /** P2P */
    void directPassOld(){
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

	FMpiBufferWriter**const sendBuffer = new FMpiBufferWriter*[nbProcess];
	memset(sendBuffer, 0, sizeof(FMpiBufferWriter*) * nbProcess);

	FMpiBufferReader**const recvBuffer = new FMpiBufferReader*[nbProcess];
	memset(recvBuffer, 0, sizeof(FMpiBufferReader*) * nbProcess);

	/* This a nbProcess x nbProcess matrix of integer
     * let U and V be id of processes :
     * globalReceiveMap[U*nbProcess + V] == size of information needed by V and own by U
     */
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

	    // Number of cells max
	    //const int limite = 1 << (this->OctreeHeight - 1);
	    // pointer to send
	    FVector<typename OctreeClass::Iterator>*const toSend = new FVector<typename OctreeClass::Iterator>[nbProcess];

	    // array that will be send to other processus for them to build the globalReceiveMap
	    int partsToSend[nbProcess];
	    memset(partsToSend, 0, sizeof(int) * nbProcess);

	    // To know if a leaf has been already sent to a proc
	    int alreadySent[nbProcess];

	    //Will store the indexes of the neighbors of current cell
	    MortonIndex indexesNeighbors[26];
	    //Obviously unused
	    //int uselessIndexArray[26];

	    for(int idxLeaf = 0 ; idxLeaf < this->numberOfLeafs ; ++idxLeaf){
		memset(alreadySent, 0, sizeof(int) * nbProcess);
		bool needOther = false;
		//Get the neighbors of current cell in indexesNeighbors, and their number in neighCount
		const int neighCount = (iterArray[idxLeaf].getCurrentGlobalCoordinate()).getNeighborsIndexes(OctreeHeight,indexesNeighbors);
		//Loop over the neighbor leafs
		for(int idxNeigh = 0 ; idxNeigh < neighCount ; ++idxNeigh){
		    //Test if leaf belongs to someone else (false if it's mine)
		    if(indexesNeighbors[idxNeigh] < (intervals[idProcess].leftIndex) || (intervals[idProcess].rightIndex) < indexesNeighbors[idxNeigh]){
			needOther = true;

			// find the proc that will need current leaf
			int procToReceive = idProcess;
			while( procToReceive != 0 && indexesNeighbors[idxNeigh] < intervals[procToReceive].leftIndex){
			    --procToReceive; //scroll process "before" current process
			}

			while( procToReceive != nbProcess - 1 && (intervals[procToReceive].rightIndex) < indexesNeighbors[idxNeigh]){
			    ++procToReceive;//scroll process "after" current process
			}
			//  Test : Not Already Send && USELESS TEST ?
			if( !alreadySent[procToReceive] && intervals[procToReceive].leftIndex <= indexesNeighbors[idxNeigh] && indexesNeighbors[idxNeigh] <= intervals[procToReceive].rightIndex){

			    alreadySent[procToReceive] = 1;
			    toSend[procToReceive].push( iterArray[idxLeaf] );
			    partsToSend[procToReceive] += iterArray[idxLeaf].getCurrentListSrc()->getSavedSize();
			    partsToSend[procToReceive] += int(sizeof(MortonIndex));
			}
		    }
		}

		if(needOther){ //means that something need to be sent (or received)
		    leafsNeedOther.set(idxLeaf,true);
		    ++countNeedOther;
		}
	    }

	    // No idea why it is mandatory there, could it be a few line before,
	    for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
		if(partsToSend[idxProc]){
		    partsToSend[idxProc] += int(sizeof(int));
		}
	    }

	    //Share to all processus globalReceiveMap
	    FLOG(gatherCounter.tic());
	    FMpi::MpiAssert( MPI_Allgather( partsToSend, nbProcess, MPI_INT, globalReceiveMap, nbProcess, MPI_INT, comm.getComm()),  __LINE__ );
	    FLOG(gatherCounter.tac());

	    //Prepare receive
	    for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
		if(globalReceiveMap[idxProc * nbProcess + idProcess]){ //if idxProc has sth for me.
		    //allocate buffer of right size
		    recvBuffer[idxProc] = new FMpiBufferReader(comm.getComm(),globalReceiveMap[idxProc * nbProcess + idProcess]);
		    FMpi::MpiAssert( MPI_Irecv(recvBuffer[idxProc]->data(), recvBuffer[idxProc]->getCapacity(), MPI_PACKED,
					       idxProc, FMpi::TagFmmP2P, comm.getComm(), &requests[iterRequest++]) , __LINE__ );
		}
	    }

	    nbMessagesToRecv = iterRequest;
	    // Prepare send
	    for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
		if(toSend[idxProc].getSize() != 0){
		    sendBuffer[idxProc] = new FMpiBufferWriter(comm.getComm(),globalReceiveMap[idProcess*nbProcess+idxProc]);
		    // << is equivalent to write().
		    (*sendBuffer[idxProc]) << toSend[idxProc].getSize();
		    for(int idxLeaf = 0 ; idxLeaf < toSend[idxProc].getSize() ; ++idxLeaf){
			(*sendBuffer[idxProc]) << toSend[idxProc][idxLeaf].getCurrentGlobalIndex();
			toSend[idxProc][idxLeaf].getCurrentListSrc()->save(*sendBuffer[idxProc]);
		    }
		    //TEST BERENGER
		    //if(sendBuffer[idxProc]->getSize() != partsToSend[idxProc]){
		    FMpi::MpiAssert( MPI_Isend( sendBuffer[idxProc]->data(), sendBuffer[idxProc]->getSize() , MPI_PACKED ,
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

		leafsDataArray[startPosAtShape[shapePosition]].coord = myLeafs[idxInArray].getCurrentGlobalCoordinate();
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
	    int previous = 0;

	    for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
		const int endAtThisShape = shapeLeaf[idxShape] + previous;

#pragma omp for
		for(int idxLeafs = previous ; idxLeafs < endAtThisShape ; ++idxLeafs){
		    LeafData& currentIter = leafsDataArray[idxLeafs];
		    myThreadkernels.L2P(currentIter.cell, currentIter.targets);

		    // need the current particles and neighbors particles
		    const int counter = tree->getLeafsNeighbors(neighbors, currentIter.coord, LeafIndex);
		    myThreadkernels.P2P( currentIter.coord,currentIter.targets,
					 currentIter.sources, neighbors, counter);

		}
		previous = endAtThisShape;
	    }
	}
	FLOG(computationCounter.tac());
	FTRACE( regionP2PTrace.end() );

	//////////////////////////////////////////////////////////
	// Waitsend receive
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
		    recvBuffer[idxProc] = nullptr;
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
	    const int nbLeafToProceed = leafsNeedOtherData.getSize();

#pragma omp for
	    for(int idxLeafs = 0 ; idxLeafs < nbLeafToProceed ; ++idxLeafs){
		LeafData currentIter = leafsNeedOtherData[idxLeafs];

		// need the current particles and neighbors particles
		int counter = 0;
		memset( neighbors, 0, sizeof(ContainerClass*) * 27);

		// Take possible data
		const int nbNeigh = currentIter.coord.getNeighborsIndexes(OctreeHeight, indexesNeighbors, indexArray);

		for(int idxNeigh = 0 ; idxNeigh < nbNeigh ; ++idxNeigh){
		    if(indexesNeighbors[idxNeigh] < (intervals[idProcess].leftIndex) || (intervals[idProcess].rightIndex) < indexesNeighbors[idxNeigh]){
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


	FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << " s)\n" );
	FLOG( FLog::Controller << "\t\t Computation L2P + P2P : " << computationCounter.elapsed() << " s\n" );
	FLOG( FLog::Controller << "\t\t Computation P2P 2 : " << computation2Counter.elapsed() << " s\n" );
	FLOG( FLog::Controller << "\t\t Prepare P2P : " << prepareCounter.elapsed() << " s\n" );
	FLOG( FLog::Controller << "\t\t Gather P2P : " << gatherCounter.elapsed() << " s\n" );
	FLOG( FLog::Controller << "\t\t Wait : " << waitCounter.elapsed() << " s\n" );

    }

    /** P2P */
    void directPass(){
	FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
	FLOG( FLog::Controller.write("\tStart Direct Pass\n").write(FLog::Flush); );
	FLOG( FTic counterTime);
	FLOG( FTic prepareCounter);
	FLOG( FTic gatherCounter);
	FLOG( FTic waitCounter);
	FLOG(FTic computation2Counter);
	FLOG(FTic computationCounter);

	///////////////////////////////////////////////////
	// Prepar data to send receive
	///////////////////////////////////////////////////
	FLOG(prepareCounter.tic());

	// To send in asynchronous way
	MPI_Request requests[2 * nbProcess];
	MPI_Status    status[2 * nbProcess];
	int iterRequest = 0;
	int nbMessagesToRecv = 0;

	FMpiBufferWriter**const sendBuffer = new FMpiBufferWriter*[nbProcess];
	memset(sendBuffer, 0, sizeof(FMpiBufferWriter*) * nbProcess);

	FMpiBufferReader**const recvBuffer = new FMpiBufferReader*[nbProcess];
	memset(recvBuffer, 0, sizeof(FMpiBufferReader*) * nbProcess);

	/* This a nbProcess x nbProcess matrix of integer
     * let U and V be id of processes :
     * globalReceiveMap[U*nbProcess + V] == size of information needed by V and own by U
     */
	int*const globalReceiveMap = new int[nbProcess * nbProcess];
	memset(globalReceiveMap, 0, sizeof(int) * nbProcess * nbProcess);

	FVector<LeafData> * leafsNeedOtherData;
	LeafData* leafsDataArray;
	OctreeClass* otherP2Ptree;

	FBoolArray leafsNeedOther(this->numberOfLeafs);
	int countNeedOther = 0;

	///////////////////////////////////////////////////
	// Prepare data for thread P2P
	///////////////////////////////////////////////////

	// init
	const int LeafIndex = OctreeHeight - 1;
	const int SizeShape = 3*3*3;

	int shapeLeaf[SizeShape];
	memset(shapeLeaf,0,SizeShape*sizeof(int));

	leafsDataArray = new LeafData[this->numberOfLeafs];

	leafsNeedOtherData = new FVector<LeafData>(countNeedOther);


	// This first part is sequential, we split the datas between
	// colors to avoid writing concurrency later with omp threads

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

		leafsDataArray[startPosAtShape[shapePosition]].coord  = myLeafs[idxInArray].getCurrentGlobalCoordinate();
		leafsDataArray[startPosAtShape[shapePosition]].cell      = myLeafs[idxInArray].getCurrentCell();
		leafsDataArray[startPosAtShape[shapePosition]].targets = myLeafs[idxInArray].getCurrentListTargets();
		leafsDataArray[startPosAtShape[shapePosition]].sources = myLeafs[idxInArray].getCurrentListSrc();
		if( leafsNeedOther.get(idxLeaf) ) leafsNeedOtherData->push(leafsDataArray[startPosAtShape[shapePosition]]);

		++startPosAtShape[shapePosition];
	    }

	    delete[] shapeType;
	    delete[] myLeafs;
	}

	//At this point, we start with the parallel section
	//One thread will be in charge of communication
	//Two comm : AllGather then iSend and IRecv
#pragma omp parallel
	{
#pragma omp single nowait
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

		// Number of cells max
		//const int limite = 1 << (this->OctreeHeight - 1);
		// pointer to send
		FVector<typename OctreeClass::Iterator>*const toSend = new FVector<typename OctreeClass::Iterator>[nbProcess];

		// array that will be send to other processus for them to build the globalReceiveMap
		int partsToSend[nbProcess];
		memset(partsToSend, 0, sizeof(int) * nbProcess);

		// To know if a leaf has been already sent to a proc
		int alreadySent[nbProcess];

		//Will store the indexes of the neighbors of current cell
		MortonIndex indexesNeighbors[26];
		//Obviously unused
		//int uselessIndexArray[26];

		for(int idxLeaf = 0 ; idxLeaf < this->numberOfLeafs ; ++idxLeaf){
		    memset(alreadySent, 0, sizeof(int) * nbProcess);
		    bool needOther = false;
		    //Get the neighbors of current cell in indexesNeighbors, and their number in neighCount
		    const int neighCount = (iterArray[idxLeaf].getCurrentGlobalCoordinate()).getNeighborsIndexes(OctreeHeight,indexesNeighbors);
		    //Loop over the neighbor leafs
		    for(int idxNeigh = 0 ; idxNeigh < neighCount ; ++idxNeigh){
			//Test if leaf belongs to someone else (false if it's mine)
			if(indexesNeighbors[idxNeigh] < (intervals[idProcess].leftIndex) || (intervals[idProcess].rightIndex) < indexesNeighbors[idxNeigh]){
			    needOther = true;

			    // find the proc that will need current leaf
			    int procToReceive = idProcess;
			    while( procToReceive != 0 && indexesNeighbors[idxNeigh] < intervals[procToReceive].leftIndex){
				--procToReceive; //scroll process "before" current process
			    }

			    while( procToReceive != nbProcess - 1 && (intervals[procToReceive].rightIndex) < indexesNeighbors[idxNeigh]){
				++procToReceive;//scroll process "after" current process
			    }
			    //  Test : Not Already Send && USELESS TEST ?
			    if( !alreadySent[procToReceive] && intervals[procToReceive].leftIndex <= indexesNeighbors[idxNeigh] && indexesNeighbors[idxNeigh] <= intervals[procToReceive].rightIndex){

				alreadySent[procToReceive] = 1;
				toSend[procToReceive].push( iterArray[idxLeaf] );
				partsToSend[procToReceive] += iterArray[idxLeaf].getCurrentListSrc()->getSavedSize();
				partsToSend[procToReceive] += int(sizeof(MortonIndex));
			    }
			}
		    }

		    if(needOther){ //means that something need to be sent (or received)
			leafsNeedOther.set(idxLeaf,true);
			++countNeedOther;
		    }
		}

		// No idea why it is mandatory there, could it be a few line before,
		for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
		    if(partsToSend[idxProc]){
			partsToSend[idxProc] += int(sizeof(int));
		    }
		}

		//Share to all processus globalReceiveMap
		FLOG(gatherCounter.tic());
		FMpi::MpiAssert( MPI_Allgather( partsToSend, nbProcess, MPI_INT, globalReceiveMap, nbProcess, MPI_INT, comm.getComm()),  __LINE__ );
		FLOG(gatherCounter.tac());

		//Prepare receive
		for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
		    if(globalReceiveMap[idxProc * nbProcess + idProcess]){ //if idxProc has sth for me.
			//allocate buffer of right size
			recvBuffer[idxProc] = new FMpiBufferReader(comm.getComm(),globalReceiveMap[idxProc * nbProcess + idProcess]);
			FMpi::MpiAssert( MPI_Irecv(recvBuffer[idxProc]->data(), recvBuffer[idxProc]->getCapacity(), MPI_PACKED,
						   idxProc, FMpi::TagFmmP2P, comm.getComm(), &requests[iterRequest++]) , __LINE__ );
		    }
		}

		nbMessagesToRecv = iterRequest;
		// Prepare send
		for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
		    if(toSend[idxProc].getSize() != 0){
			sendBuffer[idxProc] = new FMpiBufferWriter(comm.getComm(),globalReceiveMap[idProcess*nbProcess+idxProc]);
			// << is equivalent to write().
			(*sendBuffer[idxProc]) << toSend[idxProc].getSize();
			for(int idxLeaf = 0 ; idxLeaf < toSend[idxProc].getSize() ; ++idxLeaf){
			    (*sendBuffer[idxProc]) << toSend[idxProc][idxLeaf].getCurrentGlobalIndex();
			    toSend[idxProc][idxLeaf].getCurrentListSrc()->save(*sendBuffer[idxProc]);
			}
			//TEST BERENGER
			//if(sendBuffer[idxProc]->getSize() != partsToSend[idxProc]){
			FMpi::MpiAssert( MPI_Isend( sendBuffer[idxProc]->data(), sendBuffer[idxProc]->getSize() , MPI_PACKED ,
						    idxProc, FMpi::TagFmmP2P, comm.getComm(), &requests[iterRequest++]) , __LINE__ );

		    }
		}

		delete[] toSend;
		//////////////////////////////////////////////////////////
		// Waitsend receive
		//////////////////////////////////////////////////////////

		FLOG(computation2Counter.tic());

		// Create an octree with leaves from others
		otherP2Ptree = new OctreeClass( tree->getHeight(), tree->getSubHeight(), tree->getBoxWidth(), tree->getBoxCenter() );
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
				otherP2Ptree->createLeaf(leafIndex)->getSrc()->restore((*recvBuffer[idxProc]));
			    }
			    delete recvBuffer[idxProc];
			    recvBuffer[idxProc] = 0;
			}
		    }
		}
		delete[] indexMessage;

	    }//End single section

	    FLOG(prepareCounter.tac());

	    //////////////////////////////////////////////////////////
	    // Computation P2P that DO NOT need others data
	    //////////////////////////////////////////////////////////
	    FTRACE( FTrace::FRegion regionP2PTrace("Compute P2P", __FUNCTION__ , __FILE__ , __LINE__) );

	    FLOG(computationCounter.tic());

	    KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
	    // There is a maximum of 26 neighbors
	    ContainerClass* neighbors[27];
	    int previous = 0;

	    for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
		const int endAtThisShape = shapeLeaf[idxShape] + previous;

#pragma omp for schedule(auto)
		for(int idxLeafs = previous ; idxLeafs < endAtThisShape ; ++idxLeafs){
		    LeafData& currentIter = leafsDataArray[idxLeafs];
		    myThreadkernels.L2P(currentIter.cell, currentIter.targets);

		    // need the current particles and neighbors particles
		    const int counter = tree->getLeafsNeighbors(neighbors, currentIter.coord, LeafIndex);
		    myThreadkernels.P2P( currentIter.coord,currentIter.targets,
					 currentIter.sources, neighbors, counter);
		}

		previous = endAtThisShape;
	    }
	    //}


	    FLOG(computationCounter.tac());
	    FTRACE( regionP2PTrace.end() );


	    //////////////////////////////////////////////////////////
	    // Computation P2P that need others data
	    //////////////////////////////////////////////////////////

	    FTRACE( FTrace::FRegion regionOtherTrace("Compute P2P Other", __FUNCTION__ , __FILE__ , __LINE__) );

	    {
		/*KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);*/
		// There is a maximum of 26 neighbors
		/*ContainerClass* neighbors[27];*/
		MortonIndex indexesNeighbors[27];
		int indexArray[26];
		// Box limite
		const int nbLeafToProceed = leafsNeedOtherData->getSize();

#pragma omp for schedule(auto) nowait
		for(int idxLeafs = 0 ; idxLeafs < nbLeafToProceed ; ++idxLeafs){
		    LeafData currentIter = (*leafsNeedOtherData)[idxLeafs];

		    // need the current particles and neighbors particles
		    int counter = 0;
		    memset( neighbors, 0, sizeof(ContainerClass*) * 27);

		    // Take possible data
		    const int nbNeigh = currentIter.coord.getNeighborsIndexes(OctreeHeight, indexesNeighbors, indexArray);

		    for(int idxNeigh = 0 ; idxNeigh < nbNeigh ; ++idxNeigh){
			if(indexesNeighbors[idxNeigh] < (intervals[idProcess].leftIndex) || (intervals[idProcess].rightIndex) < indexesNeighbors[idxNeigh]){
			    ContainerClass*const hypotheticNeighbor = otherP2Ptree->getLeafSrc(indexesNeighbors[idxNeigh]);
			    if(hypotheticNeighbor){
				neighbors[ indexArray[idxNeigh] ] = hypotheticNeighbor;
				++counter;
			    }
			}
		    }

		    myThreadkernels.P2PRemote( currentIter.cell->getCoordinate(), currentIter.targets,
					       currentIter.sources, neighbors, counter);
		}//End For

	    }
	}//End parallel section

	for(int idxProc = 0 ; idxProc < nbProcess ; ++idxProc){
	    delete sendBuffer[idxProc];
	    delete recvBuffer[idxProc];
	}
	delete[] globalReceiveMap;
	delete[] leafsDataArray;

	FLOG(computation2Counter.tac());


	FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << " s)\n" );
	FLOG( FLog::Controller << "\t\t Computation L2P + P2P : " << computationCounter.elapsed() << " s\n" );
	FLOG( FLog::Controller << "\t\t Computation P2P 2 : " << computation2Counter.elapsed() << " s\n" );
	FLOG( FLog::Controller << "\t\t Prepare P2P : " << prepareCounter.elapsed() << " s\n" );
	FLOG( FLog::Controller << "\t\t Gather P2P : " << gatherCounter.elapsed() << " s\n" );
	FLOG( FLog::Controller << "\t\t Wait : " << waitCounter.elapsed() << " s\n" );

    }

};






#endif //FFMMALGORITHMTHREAD_HPP
