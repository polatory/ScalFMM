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
  
  typename OctreeClass::Iterator* iterArray;  //
  int numberOfLeafs;                          //< To store the size at the previous level

  const int MaxThreads;               //< the max number of thread allowed by openmp
  
  const int nbProcess;                //< Number of process
  const int idProcess;                //< Id of current process
  
  const int OctreeHeight;            //<Height of the tree

  
  /** An interval is the morton index interval
   * that a proc use (it holds data in this interval)
   */
  struct Interval{
    MortonIndex min;
    MortonIndex max;
  };
  /** My interval */
  Interval*const intervals;
  /** All process intervals */
  Interval*const workingIntervalsPerLevel;

  /** Get an interval from proc id and level */
  Interval& getWorkingInterval( int level,  int proc){
    return workingIntervalsPerLevel[OctreeHeight * proc + level];
  }

public:
  /** Get current proc interval at level */
  Interval& getWorkingInterval( int level){
    return getWorkingInterval(level, idProcess);
  }

  /** Does the current proc has some work at this level */
  bool hasWorkAtLevel( int level){
    return idProcess == 0 || (getWorkingInterval(level, idProcess - 1).max) < (getWorkingInterval(level, idProcess).max);
  }

  /** The constructor need the octree and the kernels used for computation
   * @param inTree the octree to work on
   * @param inKernels the kernels to call
   * An assert is launched if one of the arguments is null
   */
  FFmmAlgorithmThreadProc(const FMpi::FComm& inComm, OctreeClass* const inTree, KernelClass* const inKernels)
    : tree(inTree) , kernels(0), comm(inComm), iterArray(nullptr),numberOfLeafs(0),
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
      printf("Proc::%d From leaf %lld to leaf %lld\n",idProcess,myLastInterval.min,myLastInterval.max);
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
    FLOG(FTic computationCounter);
    typename OctreeClass::Iterator octreeIterator(tree);

    // Iterate on leafs
    octreeIterator.gotoBottomLeft();
    int leafs = 0;
    do{
      iterArray[leafs++] = octreeIterator;
    } while(octreeIterator.moveRight());

    FLOG(computationCounter.tic());
#pragma omp parallel
    {
      KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
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
    FMpiBufferWriter sendBuffer(comm.getComm(),7*MaxSizePerCell);
    const int recvBufferOffset = (8 * MaxSizePerCell + 1);
    FMpiBufferReader recvBuffer(comm.getComm(), nbProcess*recvBufferOffset);
    CellClass recvBufferCells[8];
    
    int firstProcThatSend = idProcess + 1;
    FLOG(computationCounter.tic());
    //Loop for work
    for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
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
      
      int cellsToSend = -1;
      int iterRequests = 0;
     
      int endIndex = numberOfCells;
      //Test if i'm not the last, and I need st to compute my last M2M
      if((idProcess != nbProcess-1) && (getWorkingInterval(idxLevel+1,idProcess+1)).min <= getWorkingInterval(idxLevel+1,idProcess).max){
	endIndex--;
      }

      while(iterArray[cellsToSend+1].getCurrentGlobalIndex() < getWorkingInterval(idxLevel, idProcess).min){
	++cellsToSend;
      }
      FLOG(parallelCounter.tic());
#pragma omp parallel
      {
	//This single section is supposed post and receive the comms, and then do the M2M associated with its.
#pragma omp single
	{
	  FLOG(singleCounter.tic());
	  //Datas needed in several parts of the section
	  bool hasToReceive;
	  int endProcThatSend;
	  
	  //Post Send
	  if(idProcess != 0
	     && (getWorkingInterval((idxLevel+1), idProcess).min >>3) <= (getWorkingInterval((idxLevel+1), idProcess - 1).max >>3)){
	    
	    char state = 0;
	    sendBuffer.write(state);
	    
	    const CellClass* const* const child = iterArray[cellsToSend].getCurrentChild();
	    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
	      if( child[idxChild] && getWorkingInterval((idxLevel+1), idProcess).min <= child[idxChild]->getMortonIndex() ){
		child[idxChild]->serializeUp(sendBuffer);
		state = char(state | (0x1 << idxChild));
		
	      }
	    }
	    sendBuffer.writeAt(0,state);
	    
	    while( sendToProc && iterArray[cellsToSend].getCurrentGlobalIndex() <= getWorkingInterval(idxLevel , sendToProc - 1).max){
	      --sendToProc;
	    }
	    
	    MPI_Isend(sendBuffer.data(), sendBuffer.getSize(), MPI_PACKED, sendToProc, 
		      FMpi::TagFmmM2M, comm.getComm(), &requests[iterRequests++]);
	    
	  }
	  //Post receive
	  {
	    // We may need to receive something
	    hasToReceive = false;
	    endProcThatSend = firstProcThatSend;
	    
	    if(idProcess != nbProcess - 1){ // if I'm the last one (idProcess == nbProcess-1), I shall not receive anything in a M2M
	      while(firstProcThatSend < nbProcess
		    && (getWorkingInterval((idxLevel+1), firstProcThatSend).max) <= (getWorkingInterval((idxLevel+1), idProcess).max)){
		// Second condition :: while firstProcThatSend max morton index is < to myself max interval
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
		    MPI_Irecv(&recvBuffer.data()[idxProc * recvBufferOffset], recvBufferOffset, MPI_PACKED,
			      idxProc, FMpi::TagFmmM2M, comm.getComm(), &requests[iterRequests++]);
		  }
		}
	      }
	    }
	  }
	  //Wait For the comms, and do the work
	  {
	    // Are we sending or waiting anything?
	    if(iterRequests){
	      MPI_Waitall( iterRequests, requests, status);
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
		(*kernels[0]).M2M( iterArray[numberOfCells - 1].getCurrentCell() , currentChild, idxLevel);
		
		firstProcThatSend = endProcThatSend - 1;
	      }
	    }
	  }
	  sendBuffer.reset();
	  recvBuffer.seek(0);
	  FLOG(singleCounter.tac());
	}//End Of Single section
		
	KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
#pragma omp for nowait
	for( int idxCell = cellsToSend+1 ; idxCell < endIndex ; ++idxCell){
	  myThreadkernels.M2M( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
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

  /** M2L  */
  void transferPass(){
    const int MaxSizePerCell = CellClass::GetSize();
    FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

    FLOG( FLog::Controller.write("\tStart Downward Pass (M2L)\n").write(FLog::Flush); );
    FLOG(FTic counterTime);
    FLOG(FTic computationCounter);
    FLOG(FTic singleCounter);
    FLOG(FTic gatherCounter);
    FLOG(FTic m2lSelf);
    FLOG(FTic m2lFar);
    FLOG(FTic sendCounter);

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

    //Variables needed in multiple omp sections
    int iterRequest;
    MPI_Request * requests;
    MPI_Status * status;
    int* globalReceiveMap;
    FMpiBufferWriter** sendBuffer;
    FMpiBufferReader** recvBuffer;

#pragma omp parallel
    {
#pragma omp single nowait 
      {
	FTRACE( FTrace::FRegion regionTrace( "Preprocess" , __FUNCTION__ , __FILE__ , __LINE__) );
	FLOG(singleCounter.tic());
	
	// To know if a leaf has been already sent to a proc
	bool*const alreadySent = new bool[nbProcess];
	memset(alreadySent, 0, sizeof(bool) * nbProcess);
	
	typename OctreeClass::Iterator octreeIterator(tree);
	octreeIterator.moveDown();
	typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
	// for each levels
	for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
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
	  MortonIndex neighborsIndexes[189];
	  for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
	    // Find the M2L neigbors of a cell
	    const int counter = iterArray[idxCell].getCurrentGlobalCoordinate().getInteractionNeighbors(idxLevel, neighborsIndexes);
	    
	    memset(alreadySent, false, sizeof(bool) * nbProcess);
	    bool needOther = false;
	    // Test each negibors to know which one do not belong to us
	    for(int idxNeigh = 0 ; idxNeigh < counter ; ++idxNeigh){
	      if(neighborsIndexes[idxNeigh] < getWorkingInterval(idxLevel , idProcess).min
		 || (getWorkingInterval(idxLevel , idProcess).max) < neighborsIndexes[idxNeigh]){
		int procToReceive = idProcess;
		while( 0 != procToReceive && neighborsIndexes[idxNeigh] < getWorkingInterval(idxLevel , procToReceive).min ){
		  --procToReceive;
		}
		while( procToReceive != nbProcess -1 && (getWorkingInterval(idxLevel , procToReceive).max) < neighborsIndexes[idxNeigh]){
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

	delete[] alreadySent;
      

	//////////////////////////////////////////////////////////////////
	// Gather this information
	//////////////////////////////////////////////////////////////////
	
	FLOG(gatherCounter.tic());
	// All process say to each others
	// what the will send to who
	globalReceiveMap = new int[nbProcess * nbProcess * OctreeHeight];
	memset(globalReceiveMap, 0, sizeof(int) * nbProcess * nbProcess * OctreeHeight);
	FMpi::MpiAssert( MPI_Allgather( indexToSend, nbProcess * OctreeHeight, MPI_INT, globalReceiveMap, nbProcess * OctreeHeight, MPI_INT, comm.getComm()),  __LINE__ );
	FLOG(gatherCounter.tac());
	

	//////////////////////////////////////////////////////////////////
	// Send and receive for real
	//////////////////////////////////////////////////////////////////
	
	FLOG(sendCounter.tic());
	// Then they can send and receive (because they know what they will receive)
	// To send in asynchrone way
	requests = new MPI_Request[2 * nbProcess * OctreeHeight];
	status = new MPI_Status[2 * nbProcess * OctreeHeight];
	iterRequest = 0;

	const int SizeOfCellToSend = int(sizeof(MortonIndex) + sizeof(int) + MaxSizePerCell);
	
	sendBuffer = new FMpiBufferWriter*[nbProcess * OctreeHeight];
	memset(sendBuffer, 0, sizeof(FMpiBufferWriter*) * nbProcess * OctreeHeight);
	
	recvBuffer = new FMpiBufferReader*[nbProcess * OctreeHeight];
	memset(recvBuffer, 0, sizeof(FMpiBufferReader*) * nbProcess * OctreeHeight);

    
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
	FLOG(sendCounter.tac());
	{
	  //////////////////////////////////////////////////////////////////
	  // Wait received data and compute
	  //////////////////////////////////////////////////////////////////
	  
	  // Wait to receive every things (and send every things)
	  MPI_Waitall(iterRequest, requests, status);
	}
	FLOG(singleCounter.tac());
      }//End of Single section
      
      //////////////////////////////////////////////////////////////////
      // Do M2L SELF
      //////////////////////////////////////////////////////////////////
      FLOG(m2lSelf.tic());
      {
	FTRACE( FTrace::FRegion regionTrace("Compute", __FUNCTION__ , __FILE__ , __LINE__) );
	typename OctreeClass::Iterator octreeIterator(tree);
	octreeIterator.moveDown();
	typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
	// Now we can compute all the data
	// for each levels
	for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
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

	  {
	    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
	    const CellClass* neighbors[343];
	    
#pragma omp for schedule(dynamic) //nowait 
	    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
	      const int counter = tree->getInteractionNeighbors(neighbors,  iterArray[idxCell].getCurrentGlobalCoordinate(), idxLevel);
	      if(counter) myThreadkernels->M2L( iterArray[idxCell].getCurrentCell() , neighbors, counter, idxLevel);
	    }
	    myThreadkernels->finishedLevelM2L(idxLevel);
	  }
	  FLOG(computationCounter.tac());
	}
      }
      FLOG(m2lSelf.tac());
    }
        
    
    FTRACE( FTrace::FRegion regionTrace("Compute Received data", __FUNCTION__ , __FILE__ , __LINE__) );
    
    typename OctreeClass::Iterator octreeIterator(tree);
    octreeIterator.moveDown();
    typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);
    // compute the second time
    // for each levels
    for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
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
      FLOG(m2lFar.tic());
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
	  const int counterNeighbors = iterArray[idxCell].getCurrentGlobalCoordinate().getInteractionNeighbors(idxLevel, neighborsIndex, neighborsPosition);

	  int counter = 0;
	  // does we receive this index from someone?
	  for(int idxNeig = 0 ;idxNeig < counterNeighbors ; ++idxNeig){
	    if(neighborsIndex[idxNeig] < (getWorkingInterval(idxLevel , idProcess).min)
	       || (getWorkingInterval(idxLevel , idProcess).max) < neighborsIndex[idxNeig]){

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
      FLOG(m2lFar.tac());
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

    FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << " s)\n" );
    FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    FLOG( FLog::Controller << "\t\t Single : " << singleCounter.cumulated() << " s\n" );
    FLOG( FLog::Controller << "\t\t M2L Self : " << m2lSelf.cumulated() << " s\n" );
    FLOG( FLog::Controller << "\t\t M2L Far : " << m2lFar.cumulated() << " s\n" );
    FLOG( FLog::Controller << "\t\t M2L Gather : " << gatherCounter.elapsed() << " s\n" );
    FLOG( FLog::Controller << "\t\t M2L Send : " << sendCounter.elapsed() << " s\n" );
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
    
    MPI_Request*const requests = new MPI_Request[nbProcess];
    MPI_Status*const status = new MPI_Status[nbProcess];
    
    const int heightMinusOne = OctreeHeight - 1;

    FMpiBufferWriter sendBuffer(comm.getComm(),MaxSizePerCell);
    FMpiBufferReader recvBuffer(comm.getComm(),MaxSizePerCell);

    // for each levels exepted leaf level
    for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
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

#pragma omp parallel
      {
#pragma omp single
	{
	  bool needToRecv = false;
	  int iterRequests = 0;

	  FLOG(prepareCounter.tic());
	  
	  // do we need to receive one or zeros cell
	  if(idProcess != 0
	     && (getWorkingInterval((idxLevel + 1) , idProcess).min >> 3 ) <= (getWorkingInterval((idxLevel+1) , idProcess - 1).max >> 3 ) ){
	    needToRecv = true;


	    MPI_Irecv( recvBuffer.data(), recvBuffer.getCapacity(), MPI_PACKED, MPI_ANY_SOURCE,
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

		MPI_Isend(sendBuffer.data(), sendBuffer.getSize(), MPI_PACKED, idxProc,
			  FMpi::TagFmmL2L, comm.getComm(), &requests[iterRequests++]);
	      }
	      
	    }
	  }
	  
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
	      
	      kernels[0]->L2L( iterArray[firstCellWork].getCurrentCell() , iterArray[firstCellWork].getCurrentChild(), idxLevel);
	      FLOG(computationCounter.tac());
	    }
	  }
	  
	}//End Of Single Section
	FLOG(prepareCounter.tac());

	FLOG(computationCounter.tic());
	//#pragma omp parallel
	KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
#pragma omp for nowait
	for(int idxCell = firstCellWork + 1 ; idxCell < numberOfCells ; ++idxCell){
	  myThreadkernels.L2L( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
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

	leafsDataArray[startPosAtShape[shapePosition]].coord = myLeafs[idxInArray].getCurrentGlobalCoordinate();
	leafsDataArray[startPosAtShape[shapePosition]].cell = myLeafs[idxInArray].getCurrentCell();
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
	    if(indexesNeighbors[idxNeigh] < (intervals[idProcess].min) || (intervals[idProcess].max) < indexesNeighbors[idxNeigh]){
	      needOther = true;
	    
	      // find the proc that will need current leaf
	      int procToReceive = idProcess;
	      while( procToReceive != 0 && indexesNeighbors[idxNeigh] < intervals[procToReceive].min){
		--procToReceive; //scroll process "before" current process
	      }
	    
	      while( procToReceive != nbProcess - 1 && (intervals[procToReceive].max) < indexesNeighbors[idxNeigh]){
		++procToReceive;//scroll process "after" current process
	      }
	      //  Test : Not Already Send && USELESS TEST ?
	      if( !alreadySent[procToReceive] && intervals[procToReceive].min <= indexesNeighbors[idxNeigh] && indexesNeighbors[idxNeigh] <= intervals[procToReceive].max){
	      
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
      
	{//TODO : remove 
	  //Print the globalReceiveMap for Process 0
	  // if(idProcess == 0)
// 	    {
// 	      printf("\n Proc 0 :: \n");
// 	      for(int u = 0 ; u < nbProcess ; ++u){
// 	        for(int v = 0 ; v < nbProcess ; ++v){
// 	  	printf("\t %d",globalReceiveMap[u*nbProcess+v]);
// 	        }
// 	        printf("\n");
// 	      }
// 	    }
	
	}
      
      
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
	// #pragma omp single
	// 	{
	// 	  char file[5];
	// 	  sprintf(file,"re%d",idProcess);
	// 	  FILE * fd = fopen(file,"a+");
	// 	  fprintf(fd,"Proc %d \t, Color : %d, from %d to %d, total %d \n",idProcess,idxShape,previous,endAtThisShape,endAtThisShape-previous);
	// 	  fclose(fd);
	// 	}
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
	    if(indexesNeighbors[idxNeigh] < (intervals[idProcess].min) || (intervals[idProcess].max) < indexesNeighbors[idxNeigh]){
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
