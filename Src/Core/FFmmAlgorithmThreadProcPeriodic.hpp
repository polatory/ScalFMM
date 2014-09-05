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

#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"
#include "../Utils/FMemUtils.hpp"

#include "../Containers/FBoolArray.hpp"
#include "../Containers/FOctree.hpp"
#include "../Containers/FLightOctree.hpp"

#include "../Containers/FBufferWriter.hpp"
#include "../Containers/FBufferReader.hpp"
#include "../Containers/FMpiBufferWriter.hpp"
#include "../Containers/FMpiBufferReader.hpp"

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
	MortonIndex leftIndex;
	MortonIndex rightIndex;
    };
    Interval*const intervals;
    Interval*const workingIntervalsPerLevel;


    Interval& getWorkingInterval(const int level, const int proc)const {
	return workingIntervalsPerLevel[OctreeHeight * proc + level];
    }

    static int GetFackLevel(const int inLevelAboveRequiered){
	if( inLevelAboveRequiered == -1 ) return 1;
	if( inLevelAboveRequiered == 0  ) return 2;
	return inLevelAboveRequiered + 3;
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
	return idProcess == 0 || getWorkingInterval(level, idProcess - 1).rightIndex < getWorkingInterval(level, idProcess).rightIndex;
    }

    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThreadProcPeriodic(const FMpi::FComm& inComm, OctreeClass* const inTree,
				    const int inUpperLevel = 2)
	: tree(inTree) , kernels(nullptr), comm(inComm), nbLevelsAboveRoot(inUpperLevel), offsetRealTree(GetFackLevel(inUpperLevel)),
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

	// Count leaf
	this->numberOfLeafs = 0;
	{

	    Interval myLastInterval;
	    {
		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.gotoBottomLeft();
		myLastInterval.leftIndex = octreeIterator.getCurrentGlobalIndex();
		do{
		    ++this->numberOfLeafs;
		} while(octreeIterator.moveRight());
		myLastInterval.rightIndex = octreeIterator.getCurrentGlobalIndex();
	    }
	    iterArray = new typename OctreeClass::Iterator[numberOfLeafs];
	    FAssertLF(iterArray, "iterArray bad alloc");

	    // We get the min/max indexes from each procs
	    FMpi::MpiAssert( MPI_Allgather( &myLastInterval, sizeof(Interval), MPI_BYTE, intervals, sizeof(Interval), MPI_BYTE, comm.getComm()),  __LINE__ );

	    Interval*const myIntervals = new Interval[OctreeHeight];
	    myIntervals[OctreeHeight - 1] = myLastInterval;
	    for(int idxLevel = OctreeHeight - 2 ; idxLevel >= 0 ; --idxLevel){
		myIntervals[idxLevel].leftIndex = myIntervals[idxLevel+1].leftIndex >> 3;
		myIntervals[idxLevel].rightIndex = myIntervals[idxLevel+1].rightIndex >> 3;
	    }
	    if(idProcess != 0){
		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.gotoBottomLeft();
		octreeIterator.moveUp();

		MortonIndex currentLimit = intervals[idProcess-1].rightIndex >> 3;

		for(int idxLevel = OctreeHeight - 2 ; idxLevel >= 1 ; --idxLevel){
		    while(octreeIterator.getCurrentGlobalIndex() <= currentLimit){
			if( !octreeIterator.moveRight() ) break;
		    }
		    myIntervals[idxLevel].leftIndex = octreeIterator.getCurrentGlobalIndex();
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

	if(operationsToProceed & FFmmM2L) transferPassOld();

	if(operationsToProceed & FFmmL2L) downardPass();

	if((operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P)) directPass();

	// delete array
	delete [] iterArray;
	iterArray = nullptr;
    }

private:

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M Bottom Pass */
    void bottomPass(){
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
	const int MaxSizePerCell = CellClass::GetSize();
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
	    const int fackLevel = idxLevel + offsetRealTree;
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
			(*kernels[threadNumber]).M2M( iterArray[totalNbCellsAtLevel - 1].getCurrentCell() , currentChild, fackLevel);
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
		    myThreadkernels->M2M( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), fackLevel);
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


	//////////////////////////////////////////////////////////////////
	//Periodicity
	//////////////////////////////////////////////////////////////////

	octreeIterator = typename OctreeClass::Iterator(tree);

	if( idProcess == 0){
	    int iterRequests = 0;

	    CellClass* currentChild[8];
	    memcpy(currentChild, octreeIterator.getCurrentBox(), 8 * sizeof(CellClass*));

	    for(int idxProc = 1 ; idxProc < nbProcess ; ++idxProc ){
		if( getWorkingInterval(1, idxProc - 1).rightIndex < getWorkingInterval(1, idxProc).rightIndex ){
		    MPI_Irecv(&recvBuffer.data()[idxProc * recvBufferOffset], recvBufferOffset, MPI_BYTE, idxProc,
			      FMpi::TagFmmM2M, comm.getComm(), &requests[iterRequests++]);
		}
	    }

	    MPI_Waitall( iterRequests, requests, MPI_STATUSES_IGNORE);

	    // retreive data and merge my child and the child from others
	    for(int idxProc = 1 ; idxProc < nbProcess ; ++idxProc){
		if( getWorkingInterval(1, idxProc - 1).rightIndex < getWorkingInterval(1, idxProc).rightIndex ){
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
		const int firstChild = getWorkingInterval(1, idProcess).leftIndex & 7;
		const int lastChild = getWorkingInterval(1, idProcess).rightIndex & 7;

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

    void transferPass(){
	const int MaxSizePerCell = CellClass::GetSize();
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
		    FLOG(prepareCounter.tic());

		    std::unique_ptr<typename OctreeClass::Iterator[]> iterArrayLocal(new typename OctreeClass::Iterator[numberOfLeafs]);

		    // To know if a leaf has been already sent to a proc
		    bool*const alreadySent = new bool[nbProcess];
		    memset(alreadySent, 0, sizeof(bool) * nbProcess);

		    typename OctreeClass::Iterator octreeIterator(tree);
		    //octreeIterator.moveDown();
		    typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
		    // for each levels
		    for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
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
			int neighborsPosition[189];
			MortonIndex neighborsIndexes[189];
			for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
			    // Find the M2L neigbors of a cell
			    const int counter =getPeriodicInteractionNeighbors(iterArray[idxCell].getCurrentGlobalCoordinate(),
									       idxLevel,
									       neighborsIndexes, neighborsPosition, AllDirs);

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

		for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
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
		//octreeIterator.moveDown();
		typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
		// Now we can compute all the data
		// for each levels
		for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
		    const int fackLevel = idxLevel + offsetRealTree;
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
				    const int counter = tree->
					getPeriodicInteractionNeighbors(neighbors,
									iterArray[idxCellToCompute].getCurrentGlobalCoordinate(),
									idxLevel, AllDirs);
				    if(counter)
					myThreadkernels->M2L( iterArray[idxCellToCompute].getCurrentCell() , neighbors, counter, fackLevel);
				}
			    }
			}
		    }

#pragma omp taskwait

		    for(int idxThread = 0 ; idxThread < omp_get_num_threads() ; ++idxThread){
#pragma omp task
			{
			    kernels[idxThread]->finishedLevelM2L(fackLevel);
			}
		    }

		    FLOG(computationCounter.tac());
		}
	    }
	}


	{
	    FLOG(receiveCounter.tic());
	    typename OctreeClass::Iterator octreeIterator(tree);
	    octreeIterator.moveDown();
	    typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
	    // compute the second time
	    // for each levels
	    for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
		const int fackLevel = idxLevel + offsetRealTree;
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


	FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << " s)\n" );
	FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Gather : " << gatherCounter.cumulated() << " s\n" );
	FLOG( FLog::Controller << "\t\t Prepare : " << prepareCounter.cumulated() << " s\n" );

    }
    /** M2L  */
    void transferPassOld(){
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
	    FLOG(prepareCounter.tic());

	    // To know if a leaf has been already sent to a proc
	    bool*const alreadySent = new bool[nbProcess];
	    memset(alreadySent, 0, sizeof(bool) * nbProcess);

	    typename OctreeClass::Iterator octreeIterator(tree);
	    typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
	    // for each levels
	    for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
		if(idProcess != 0
			&& getWorkingInterval(idxLevel, idProcess).rightIndex <= getWorkingInterval(idxLevel, idProcess - 1).rightIndex){
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
			if(neighborsIndexes[idxNeigh] < getWorkingInterval(idxLevel , idProcess).leftIndex
				|| getWorkingInterval(idxLevel , idProcess).rightIndex < neighborsIndexes[idxNeigh]){
			    int procToReceive = idProcess;
			    while( 0 != procToReceive && neighborsIndexes[idxNeigh] < getWorkingInterval(idxLevel , procToReceive).leftIndex ){
				--procToReceive;
			    }
			    while( procToReceive != nbProcess -1 && getWorkingInterval(idxLevel , procToReceive).rightIndex < neighborsIndexes[idxNeigh]){
				++procToReceive;
			    }

			    // Maybe already sent to that proc?
			    if( !alreadySent[procToReceive]
				    && getWorkingInterval(idxLevel , procToReceive).leftIndex <= neighborsIndexes[idxNeigh]
				    && neighborsIndexes[idxNeigh] <= getWorkingInterval(idxLevel , procToReceive).rightIndex){

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
	    typename OctreeClass::Iterator octreeIterator(tree);
	    typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
	    // Now we can compute all the data
	    // for each levels
	    for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
		const int fackLevel = idxLevel + offsetRealTree;
		if(idProcess != 0
			&& getWorkingInterval(idxLevel, idProcess).rightIndex <= getWorkingInterval(idxLevel, idProcess - 1).rightIndex){

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
	    FLOG(receiveCounter.tic());
	    typename OctreeClass::Iterator octreeIterator(tree);
	    typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
	    // compute the second time
	    // for each levels
	    for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
		const int fackLevel = idxLevel + offsetRealTree;
		if(idProcess != 0
			&& getWorkingInterval(idxLevel, idProcess).rightIndex <= getWorkingInterval(idxLevel, idProcess - 1).rightIndex){

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

		    #pragma omp for schedule(dynamic) nowait
		    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
			// compute indexes
			memset(neighbors, 0, 343 * sizeof(CellClass*));
			const int counterNeighbors = getPeriodicInteractionNeighbors(iterArray[idxCell].getCurrentGlobalCoordinate(),
										     idxLevel,
										     neighborsIndex, neighborsPosition, AllDirs);
			int counter = 0;
			// does we receive this index from someone?
			for(int idxNeig = 0 ;idxNeig < counterNeighbors ; ++idxNeig){

			    if(neighborsIndex[idxNeig] < getWorkingInterval(idxLevel , idProcess).leftIndex
				    || getWorkingInterval(idxLevel , idProcess).rightIndex < neighborsIndex[idxNeig]){

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

    void downardPass(){ // second L2L
	const int MaxSizePerCell = CellClass::GetSize();
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
	for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
	    const int fackLevel = idxLevel + offsetRealTree;
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
			kernels[threadNumber]->L2L( iterArray[nbCellsToSkip-1].getCurrentCell() , iterArray[nbCellsToSkip-1].getCurrentChild(), fackLevel);
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
		    myThreadkernels->L2L( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), fackLevel);
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
	MortonIndex index;
	CellClass* cell;
	ContainerClass* targets;
	ContainerClass* sources;
    };
    /** P2P */
    void directPass(){
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
		    if(indexesNeighbors[idxNeigh] < intervals[idProcess].leftIndex || intervals[idProcess].rightIndex < indexesNeighbors[idxNeigh]){
			needOther = true;

			// find the proc that need this information
			int procToReceive = idProcess;
			while( procToReceive != 0 && indexesNeighbors[idxNeigh] < intervals[procToReceive].leftIndex){
			    --procToReceive;
			}

			while( procToReceive != nbProcess - 1 && intervals[procToReceive].rightIndex < indexesNeighbors[idxNeigh]){
			    ++procToReceive;
			}

			if( !alreadySent[procToReceive] && intervals[procToReceive].leftIndex <= indexesNeighbors[idxNeigh] && indexesNeighbors[idxNeigh] <= intervals[procToReceive].rightIndex){

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
				neighbors[idxNeig] = nullptr;
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
		    recvBuffer[idxProc] = nullptr;
		}
	    }
	}
	delete[] indexMessage;

	//////////////////////////////////////////////////////////
	// Computation P2P that need others data
	//////////////////////////////////////////////////////////

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
		    if(indexesNeighbors[idxNeigh] < intervals[idProcess].leftIndex || intervals[idProcess].rightIndex < indexesNeighbors[idxNeigh]){
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
