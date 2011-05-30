#ifndef FFMMALGORITHMTHREADPROC_HPP
#define FFMMALGORITHMTHREADPROC_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"

#include "../Containers/FBoolArray.hpp"
#include "../Containers/FOctree.hpp"
#include "../Containers/FBufferVector.hpp"

#include "../Utils/FMpi.hpp"

#include <omp.h>

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
template<template< class ParticleClass, class CellClass, int OctreeHeight> class KernelClass,
class ParticleClass, class CellClass,
template<class ParticleClass> class LeafClass,
int OctreeHeight, int SubtreeHeight>
        class FFmmAlgorithmThreadProc : protected FAssertable {
    // To reduce the size of variable type based on foctree in this file
    typedef FOctree<ParticleClass, CellClass, LeafClass, OctreeHeight, SubtreeHeight> Octree;
    typedef typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator OctreeIterator;
    typedef KernelClass<ParticleClass, CellClass, OctreeHeight> Kernel;

    FMpi& app;                          //< The app to communicate

    Octree* const tree;                 //< The octree to work on
    Kernel** kernels;                   //< The kernels

    OctreeIterator* iterArray;          //< To store the iterator
    OctreeIterator* previousIterArray;  //< To store the previous iterator

    int leafLeft;                   //< To store the left limit at the previous level
    int leafRight;                  //< To store the right limit at the previous level
    int numberOfLeafs;               //< To store the size at the previous level

    int leftOffsets[OctreeHeight];      //< the right limit at different level
    int rightOffsets[OctreeHeight];     //< the left limit at different level

    const int MaxThreads;               //< the max number of thread allowed by openmp

    const int nbProcess;                //< Number of process
    const int idPorcess;                //< Id of current process

    const static int BufferSize = 2000;      //< To know max of the buffer we receive
    FBufferVector<BufferSize> * sendBuffer;  //< To put data to send into a buffer

    /** To swap between two arrays
      * the current and the previous
      */
    void swapArray(){
        OctreeIterator* const temp = iterArray;
        iterArray = previousIterArray;
        previousIterArray = temp;
    }

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThreadProc(FMpi& inApp, Octree* const inTree, Kernel* const inKernels)
        : app(inApp), tree(inTree) , kernels(0), iterArray(0),
        previousIterArray(0), leafLeft(0),leafRight(0), numberOfLeafs(0),
        MaxThreads(omp_get_max_threads()), nbProcess(inApp.processCount()), idPorcess(inApp.processId()),
        sendBuffer(0) {

        assert(tree, "tree cannot be null", __LINE__, __FILE__);

        this->kernels = new Kernel*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass<ParticleClass, CellClass, OctreeHeight>(*inKernels);
        }

        this->sendBuffer = new FBufferVector<BufferSize>[nbProcess];

        FDEBUG(FDebug::Controller << "FFmmAlgorithmThreadProc\n");
        FDEBUG(FDebug::Controller << "Max threads = "  << MaxThreads << ", Procs = " << app.processCount() << ".\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmThreadProc(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;

        delete [] this->sendBuffer;
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );

        // Count leaf
        int leafs = 0;
        OctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            ++leafs;
        } while(octreeIterator.moveRight());

        iterArray = new OctreeIterator[leafs];
        assert(iterArray, "iterArray bad alloc", __LINE__, __FILE__);

        previousIterArray = new OctreeIterator[leafs];
        assert(previousIterArray, "previousIterArray bad alloc", __LINE__, __FILE__);

        // init offsets
        for(int idxOff = 0 ; idxOff < OctreeHeight ; ++idxOff){
            leftOffsets[idxOff] = 0;
            rightOffsets[idxOff] = 0;
        }

        // run
        bottomPass();

        upwardPass();

        downardPass();

        directPass();

        // delete array
        delete [] iterArray;
        iterArray = 0;
        delete [] previousIterArray;
        previousIterArray = 0;

        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Utils functions
    /////////////////////////////////////////////////////////////////////////////

    int getLeft(const int inSize) const {
        const float step = (float(inSize) / nbProcess);
        return int(FMath::Ceil(step * idPorcess));
    }

    int getRight(const int inSize) const {
        const float step = (float(inSize) / nbProcess);
        const int res = int(FMath::Ceil(step * (idPorcess+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    int getOtherRight(const int inSize,const int other) const {
        const float step = (float(inSize) / nbProcess);
        const int res = int(FMath::Ceil(step * (other+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    int getProc(const int position, const int inSize) const {
        const float step = (float(inSize) / nbProcess);
        return int(position/step);
    }

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M */
    void bottomPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);

        OctreeIterator octreeIterator(tree);

        int leafs = 0;

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            iterArray[leafs] = octreeIterator;
            ++leafs;
        } while(octreeIterator.moveRight());

        const int startIdx = getLeft(leafs);
        const int endIdx = getRight(leafs);

        FDEBUG(FTic computationCounter);
        #pragma omp parallel
        {
            Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for nowait
            for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                // We need the current cell that represent the leaf
                // and the list of particles
                myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentListSrc());
            }
        }
        FDEBUG(computationCounter.tac());

        swapArray();
        this->leafLeft = startIdx;
        this->leafRight = endIdx - 1;
        this->numberOfLeafs = leafs;

        FDEBUG( FDebug::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic computationCounter);
        FDEBUG(FTic sendCounter);
        FDEBUG(FTic receiveCounter);

        // Start from leal level - 1
        OctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();
        OctreeIterator avoidGotoLeftIterator(octreeIterator);

        int previousLeft = this->leafLeft;
        int previousRight = this->leafRight;
        int previousSize = this->numberOfLeafs;

        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){

            int numberOfCells = 0;
            // for each cells
            do{
                iterArray[numberOfCells] = octreeIterator;
                ++numberOfCells;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();

            const int startIdx = getLeft(numberOfCells);
            const int endIdx = getRight(numberOfCells);

            if(startIdx < numberOfCells){
                FDEBUG(sendCounter.tic());
                int leftOffset = 0;
                {
                    const MortonIndex MostLeftChild = iterArray[startIdx].getCurrentGlobalIndex() << 3;
                    const MortonIndex leftChildIter = previousIterArray[previousLeft].getCurrentGlobalIndex();

                    if(leftChildIter < MostLeftChild){
                        int parentOffset = startIdx - 1;
                        MortonIndex parentIndex = iterArray[parentOffset].getCurrentGlobalIndex();
                        MortonIndex childIndex = 0;
                        while( (childIndex = previousIterArray[previousLeft+leftOffset].getCurrentGlobalIndex()) < MostLeftChild){
                            childIndex >>= 3;

                            while(childIndex != parentIndex){
                                if(childIndex < parentIndex) --parentOffset;
                                else ++parentOffset;

                                parentIndex = iterArray[parentOffset].getCurrentGlobalIndex();

                            }

                            const int idxReceiver = getProc(parentOffset,numberOfCells);
                            if(!this->sendBuffer[idxReceiver].addEnoughSpace(previousIterArray[previousLeft+leftOffset].getCurrentCell()->bytesToSendUp())){
                                app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                                this->sendBuffer[idxReceiver].clear();
                            }
                            this->sendBuffer[idxReceiver].addDataUp(previousLeft+leftOffset,*previousIterArray[previousLeft+leftOffset].getCurrentCell());

                            ++leftOffset;
                        }
                    }
                    else if(previousLeft > 0 && leftChildIter > MostLeftChild){
                        while( previousIterArray[previousLeft+leftOffset - 1].getCurrentGlobalIndex() >= MostLeftChild){
                            --leftOffset;
                        }
                    }
                    leftOffsets[idxLevel+1] = leftOffset;
                }

                int rightOffset = 0;
                {
                    const MortonIndex MostRightChild = (iterArray[endIdx-1].getCurrentGlobalIndex() << 3) | 7;
                    const MortonIndex rightChildIter = previousIterArray[previousRight].getCurrentGlobalIndex();

                    if(previousRight < previousSize - 1 && rightChildIter < MostRightChild){
                        while( previousIterArray[previousRight-rightOffset+1].getCurrentGlobalIndex() <= MostRightChild){
                            --rightOffset;
                        }
                    }
                    else if(rightChildIter > MostRightChild){
                        int parentOffset = endIdx;
                        MortonIndex parentIndex = iterArray[parentOffset].getCurrentGlobalIndex();
                        MortonIndex childIndex = 0;
                        while( (childIndex = previousIterArray[previousRight-rightOffset].getCurrentGlobalIndex()) > MostRightChild){
                            childIndex >>= 3;
                            while(childIndex != parentIndex){
                                if(childIndex < parentIndex) --parentOffset;
                                else ++parentOffset;
                                parentIndex = iterArray[parentOffset].getCurrentGlobalIndex();
                            }
                            const int idxReceiver = getProc(parentOffset,numberOfCells);
                            if(!this->sendBuffer[idxReceiver].addEnoughSpace(previousIterArray[previousRight-rightOffset].getCurrentCell()->bytesToSendUp())){
                                app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                                this->sendBuffer[idxReceiver].clear();
                            }
                            this->sendBuffer[idxReceiver].addDataUp(previousRight-rightOffset,*previousIterArray[previousRight-rightOffset].getCurrentCell());

                            ++rightOffset;
                        }
                    }
                    rightOffsets[idxLevel+1] = rightOffset;
                }
                FDEBUG(sendCounter.tac());

                #pragma omp parallel
                {
                    // send computed data
                    #pragma omp single nowait
                    {
                        if(rightOffset > 0 || leftOffset > 0){
                            for(int idxProc = 0 ; idxProc < this->nbProcess ; ++idxProc){
                                if(this->sendBuffer[idxProc].getSize()){
                                    app.sendData(idxProc,this->sendBuffer[idxProc].getSize(),this->sendBuffer[idxProc].getData(),idxLevel);
                                    this->sendBuffer[idxProc].clear();
                                }
                            }
                        }
                    }

                    // received computed data
                    #pragma omp single
                    {
                        FDEBUG(receiveCounter.tic());
                        int needToReceive = FMath::Max(0,-rightOffset) + FMath::Max(0,-leftOffset);

                        int source = 0, filled = 0;
                        int position;
                        char buffer[BufferSize];

                        while(needToReceive){
                            app.receiveData( BufferSize, idxLevel, buffer, &source, &filled);
                            for(int idxBuff = 0 ; idxBuff < filled;){
                                memcpy(&position,&buffer[idxBuff],sizeof(int));
                                idxBuff += sizeof(int);
                                idxBuff += previousIterArray[position].getCurrentCell()->readUp(&buffer[idxBuff],filled-idxBuff);
                                --needToReceive;
                            }
                        }
                        FDEBUG(receiveCounter.tac());
                    }

                    #pragma omp single nowait
                    {
                        FDEBUG(computationCounter.tic());
                    }

                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                    #pragma omp for nowait
                    for(int idxCell = startIdx ; idxCell < endIdx ; ++idxCell){
                        myThreadkernels->M2M( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
                    }

                    #pragma omp single nowait
                    {
                        FDEBUG(computationCounter.tac());
                    }

                }
            }
            else {
                int parentOffset = numberOfCells - 1;
                MortonIndex parentIndex = iterArray[parentOffset].getCurrentGlobalIndex();

                for(int idxCell = previousRight ; idxCell >= previousLeft ; --idxCell){
                    const MortonIndex childIndex = previousIterArray[idxCell].getCurrentGlobalIndex() >> 3;
                    while(childIndex != parentIndex){
                        --parentOffset;
                        parentIndex = iterArray[parentOffset].getCurrentGlobalIndex();
                    }
                    const int idxReceiver = getProc(parentOffset,numberOfCells);
                    if(!this->sendBuffer[idxReceiver].addEnoughSpace(previousIterArray[idxCell].getCurrentCell()->bytesToSendUp())){
                        app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                        this->sendBuffer[idxReceiver].clear();
                    }
                    this->sendBuffer[idxReceiver].addDataUp(idxCell,*previousIterArray[idxCell].getCurrentCell());
                }

                for(int idxProc = 0 ; idxProc < this->nbProcess ; ++idxProc){
                    if(this->sendBuffer[idxProc].getSize()){
                        app.sendData(idxProc,this->sendBuffer[idxProc].getSize(),this->sendBuffer[idxProc].getData(),idxLevel);
                        this->sendBuffer[idxProc].clear();
                    }
                }

                leftOffsets[idxLevel+1] = (previousRight-previousLeft) + 1;
            }

            swapArray();
            previousLeft = startIdx;
            previousRight = endIdx - 1;
            previousSize = numberOfCells;

        }

        app.processBarrier();

        FDEBUG( FDebug::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downard
    /////////////////////////////////////////////////////////////////////////////

    /** M2L L2L */
    void downardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );

        { // first M2L
            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            FDEBUG(FTic computationCounter);
            FDEBUG(FTic sendCounter);
            FDEBUG(FTic receiveCounter);
            FDEBUG(FTic waitingToReceiveCounter);
            FDEBUG(FTic waitSendCounter);
            FDEBUG(FTic findCounter);

            OctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();
            OctreeIterator avoidGotoLeftIterator(octreeIterator);

            FBoolArray* const indexToReceive = new FBoolArray(1 << (3*(OctreeHeight-1)));

            struct LimitCell { int counter; const CellClass* neighbors[208]; };
            LimitCell ** const unfinishedCells = new LimitCell*[this->leafRight - this->leafLeft + 1];

            FBoolArray* alreadySent[this->nbProcess];
            MortonIndex upperLimitForProc[this->nbProcess];
            for(int idxProc = 0 ; idxProc < nbProcess; ++idxProc){
                alreadySent[idxProc] = new FBoolArray(this->leafRight - this->leafLeft + 1);
            }

            // for each levels
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){

                int numberOfCells = 0;
                // for each cells
                do{
                    iterArray[numberOfCells] = octreeIterator;
                    ++numberOfCells;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                const int startIdx = getLeft(numberOfCells);
                const int endIdx = getRight(numberOfCells);

                // Get limit in term of morton index
                const MortonIndex startIdxIndex = iterArray[startIdx].getCurrentGlobalIndex();
                const MortonIndex endIdxIndex = iterArray[endIdx-1].getCurrentGlobalIndex();

                // reset arrays
                for(int idxProc = 0 ; idxProc < nbProcess; ++idxProc){
                    alreadySent[idxProc]->setToZeros();
                    upperLimitForProc[idxProc] = iterArray[getOtherRight(numberOfCells,idxProc)-1].getCurrentGlobalIndex();
                }
                indexToReceive->setToZeros();


                // to count missing data
                int needToReceive = 0;

                #pragma omp parallel
                {
                    const CellClass* neighbors[208];
                    MortonIndex neighborsIndexes[208];
                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];

                    #pragma omp single nowait
                    {
                        FDEBUG(sendCounter.tic());
                        FDEBUG(computationCounter.tic());
                    }

                    #pragma omp for //schedule(dynamic)
                    for(int idxCell = startIdx ; idxCell < endIdx ; ++idxCell){
                        const int neighborsCounter = tree->getDistantNeighborsWithIndex(neighbors, neighborsIndexes, iterArray[idxCell].getCurrentGlobalIndex(),idxLevel);
                        bool needData = false;


                        for(int idxNeighbor = 0 ; idxNeighbor < neighborsCounter ; ++idxNeighbor){
                            // Get Morton index from this neigbor cell
                            const MortonIndex indexCell = neighborsIndexes[idxNeighbor];
                            // Is this cell out of our interval?
                            if(indexCell < startIdxIndex || endIdxIndex < indexCell){
                                FDEBUG(findCounter.tic());
                                // Yes we need to send the center of computation
                                // but who this cell belong to?
                                int idxReceiver = this->idPorcess;

                                if(indexCell < startIdxIndex){
                                    --idxReceiver;
                                    while(idxReceiver && indexCell <= upperLimitForProc[idxReceiver-1]){
                                        --idxReceiver;
                                    }
                                }
                                else {
                                    ++idxReceiver;
                                    while(upperLimitForProc[idxReceiver] < indexCell){
                                        ++idxReceiver;
                                    }
                                }

                                FDEBUG(findCounter.tac());

                                FDEBUG(computationCounter.tac());
                                #pragma omp critical(CheckToSend)
                                {
                                    if(!alreadySent[idxReceiver]->get(idxCell-startIdx)){
                                        alreadySent[idxReceiver]->set(idxCell-startIdx,true);
                                        needData = true;

                                        if(!this->sendBuffer[idxReceiver].addEnoughSpace(iterArray[idxCell].getCurrentCell()->bytesToSendUp())){
                                            FDEBUG(waitSendCounter.tic());
                                            app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                                            FDEBUG(waitSendCounter.tac());
                                            this->sendBuffer[idxReceiver].clear();
                                        }
                                        this->sendBuffer[idxReceiver].addDataUp(idxCell,*iterArray[idxCell].getCurrentCell());
                                    }
                                }
                                #pragma omp critical(CheckToReceive)
                                {
                                    if(!indexToReceive->get(indexCell)){
                                        ++needToReceive;
                                        indexToReceive->set(indexCell,true);
                                        needData = true;
                                    }
                                }
                                FDEBUG(computationCounter.tic());
                            }
                        }
                        FDEBUG(computationCounter.tac());
                        // if need data we can not compute now
                        if(needData){
                            const int currentCell = idxCell - startIdx;
                            unfinishedCells[currentCell] = new LimitCell();
                            unfinishedCells[currentCell]->counter = neighborsCounter;
                            memcpy(unfinishedCells[currentCell]->neighbors,neighbors,sizeof(CellClass*)*neighborsCounter);
                            alreadySent[idPorcess]->set(idxCell-startIdx,true);
                        }
                        // we can compute now !
                        else if(neighborsCounter){
                            FDEBUG(computationCounter.tic());
                            myThreadkernels->M2L(  iterArray[idxCell].getCurrentCell() , neighbors, neighborsCounter, idxLevel);
                            FDEBUG(computationCounter.tac());
                        }
                        FDEBUG(computationCounter.tic());
                    }
                    FDEBUG(computationCounter.tac());

                    #pragma omp single nowait
                    {
                        FDEBUG(waitSendCounter.tic());
                        for(int idxProc = 0 ; idxProc < this->nbProcess ; ++idxProc){
                            if(this->sendBuffer[idxProc].getSize()){
                                app.sendData(idxProc,this->sendBuffer[idxProc].getSize(),this->sendBuffer[idxProc].getData(),idxLevel);
                                this->sendBuffer[idxProc].clear();
                            }
                        }
                        FDEBUG(waitSendCounter.tac());

                        FDEBUG(sendCounter.tac());
                    }

                    // received computed data
                    #pragma omp single
                    {
                        FDEBUG(receiveCounter.tic());


                        //FDEBUG( FDebug::Controller << "\t\tNeed to receive "  << needToReceive << " cells.\n" );

                        int source = 0, filled = 0;
                        int position;
                        char buffer[BufferSize];

                        while(needToReceive){
                            FDEBUG(waitingToReceiveCounter.tic());
                            app.receiveData( BufferSize, idxLevel, buffer, &source, &filled);
                            FDEBUG(waitingToReceiveCounter.tac());
                            for(int idxBuff = 0 ; idxBuff < filled;){
                                memcpy(&position,&buffer[idxBuff],sizeof(int));
                                idxBuff += sizeof(int);
                                idxBuff += iterArray[position].getCurrentCell()->readUp(&buffer[idxBuff],filled-idxBuff);
                                --needToReceive;
                            }
                        }

                        FDEBUG(receiveCounter.tac());
                    }

                    #pragma omp single nowait
                    {
                        FDEBUG(computationCounter.tic());
                    }

                    #pragma omp for nowait
                    for(int idxCell = startIdx ; idxCell < endIdx ; ++idxCell){
                        if(alreadySent[idPorcess]->get(idxCell-startIdx)){
                            myThreadkernels->M2L(  iterArray[idxCell].getCurrentCell() , unfinishedCells[idxCell-startIdx]->neighbors, unfinishedCells[idxCell-startIdx]->counter, idxLevel);
                            delete unfinishedCells[idxCell-startIdx];
                        }
                    }

                    #pragma omp single nowait
                    {
                        FDEBUG(computationCounter.tac());
                    }
                }
            }

            for(int idxProc = 0 ; idxProc < nbProcess; ++idxProc){
                delete alreadySent[idxProc];
            }

            delete [] unfinishedCells;
            delete indexToReceive;

            app.processBarrier();

            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Wait data to Receive : " << waitingToReceiveCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t\tTotal time to send in mpi "  << waitSendCounter.cumulated() << " s.\n" );
            FDEBUG( FDebug::Controller << "\t\t\tTotal time to find "  << findCounter.cumulated() << " s.\n" );
        }


        { // second L2L
            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            FDEBUG(FTic computationCounter);
            FDEBUG(FTic sendCounter);
            FDEBUG(FTic receiveCounter);

            OctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();

            OctreeIterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = OctreeHeight - 1;

            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){

                int numberOfCells = 0;
                // for each cells
                do{
                    iterArray[numberOfCells] = octreeIterator;
                    ++numberOfCells;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                const int startIdx = getLeft(numberOfCells);
                const int endIdx = getRight(numberOfCells);

                const int currentLeft = startIdx;
                const int currentRight = endIdx -1;

                #pragma omp parallel
                {
                    // send computed data
                    #pragma omp single nowait
                    {
                        FDEBUG(sendCounter.tic());
                        const int leftOffset = -leftOffsets[idxLevel];
                        for(int idxCell = 1 ; idxCell <= leftOffset ; ++idxCell){
                            const int idxReceiver = getProc((currentLeft-idxCell),numberOfCells);
                            if(!this->sendBuffer[idxReceiver].addEnoughSpace(iterArray[currentLeft-idxCell].getCurrentCell()->bytesToSendDown())){
                                app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                                this->sendBuffer[idxReceiver].clear();
                            }
                            this->sendBuffer[idxReceiver].addDataDown(currentLeft-idxCell,*iterArray[currentLeft-idxCell].getCurrentCell());
                        }
                        const int rightOffset = -rightOffsets[idxLevel];
                        for(int idxCell = 1 ; idxCell <= rightOffset ; ++idxCell){
                            const int idxReceiver = getProc((currentRight+idxCell),numberOfCells);
                            if(!this->sendBuffer[idxReceiver].addEnoughSpace(iterArray[currentRight+idxCell].getCurrentCell()->bytesToSendDown())){
                                app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                                this->sendBuffer[idxReceiver].clear();
                            }
                            this->sendBuffer[idxReceiver].addDataDown(currentRight+idxCell,*iterArray[currentRight+idxCell].getCurrentCell());
                        }

                        for(int idxProc = 0 ; idxProc < this->nbProcess ; ++idxProc){
                            if(this->sendBuffer[idxProc].getSize()){
                                app.sendData(idxProc,this->sendBuffer[idxProc].getSize(),this->sendBuffer[idxProc].getData(),idxLevel);
                                this->sendBuffer[idxProc].clear();
                            }
                        }

                        FDEBUG(sendCounter.tac());
                    }

                    // received computed data
                    #pragma omp single
                    {
                        FDEBUG(receiveCounter.tic());
                        int needToReceive = FMath::Max(0,rightOffsets[idxLevel]) + FMath::Max(0,leftOffsets[idxLevel]);
                        //FDEBUG( FDebug::Controller << "\t\tNeed to receive "  << needToReceive << " cells.\n" );

                        int source = 0, filled = 0;
                        int position;
                        char buffer[BufferSize];

                        while(needToReceive){
                            app.receiveData( BufferSize, idxLevel, buffer, &source, &filled);
                            for(int idxBuff = 0 ; idxBuff < filled;){
                                memcpy(&position,&buffer[idxBuff],sizeof(int));
                                idxBuff += sizeof(int);
                                idxBuff += iterArray[position].getCurrentCell()->readDown(&buffer[idxBuff],filled-idxBuff);
                                --needToReceive;
                            }
                        }

                        FDEBUG(receiveCounter.tac());
                    }
                }

                if(idxLevel != heightMinusOne){
                    FDEBUG(computationCounter.tic());
                    #pragma omp parallel
                    {
                        Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                        #pragma omp for
                        for(int idxCell = startIdx ; idxCell < endIdx ; ++idxCell){
                            myThreadkernels->L2L( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
                        }
                    }
                    FDEBUG(computationCounter.tac());
                }
            }

            app.processBarrier();

            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
        }

        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Direct
    /////////////////////////////////////////////////////////////////////////////

    /** P2P */
    void directPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);

        // init
        const int LeafIndex = OctreeHeight - 1;
        const int SizeShape = 3*3*3;
        int shapeLeaf[SizeShape];
        OctreeIterator* shapeArray[SizeShape];
        for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
            shapeLeaf[idxShape] = 0;
        }

        // split data
        {
            OctreeIterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();

            // remove useless leafs
            for(int idxLeaf = 0 ; idxLeaf < this->leafLeft ; ++idxLeaf){
                octreeIterator.moveRight();
            }

            // to store which shape for each leaf
            int* const shapeType = new int [this->leafRight - this->leafLeft + 1];

            for(int idxLeaf = this->leafLeft ; idxLeaf <= this->leafRight ; ++idxLeaf){
                iterArray[idxLeaf] = octreeIterator;

                const MortonIndex index = octreeIterator.getCurrentGlobalIndex();
                FTreeCoordinate coord;
                coord.setPositionFromMorton(index, LeafIndex);
                const int shape = (coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3);
                shapeType[idxLeaf-this->leafLeft] = shape;
                ++shapeLeaf[shape];

                octreeIterator.moveRight();
            }

            // init iter array
            int countShape[SizeShape];
            for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
                shapeArray[idxShape] = new OctreeIterator[shapeLeaf[idxShape]];
                countShape[idxShape] = 0;
            }

            // store leafs
            for(int idxLeaf = this->leafLeft ; idxLeaf <= this->leafRight ; ++idxLeaf){
                const int idxShape = shapeType[idxLeaf - this->leafLeft];
                shapeArray[idxShape][countShape[idxShape]++] = iterArray[idxLeaf];
            }

            delete[] shapeType;
        }

        FDEBUG(FTic computationCounter);

        const int startIdx = this->leafLeft;

        #pragma omp parallel
        {
            Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
            // There is a maximum of 26 neighbors
            FList<ParticleClass>* neighbors[26];
            MortonIndex neighborsIndex[26];

            for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
                const int leafAtThisShape = shapeLeaf[idxShape];

                #pragma omp for
                for(int idxLeafs = startIdx ; idxLeafs < leafAtThisShape ; ++idxLeafs){
                    OctreeIterator currentIter = shapeArray[idxShape][idxLeafs];
                    myThreadkernels->L2P(currentIter.getCurrentCell(), currentIter.getCurrentListTargets());
                    // need the current particles and neighbors particles
                    const int counter = tree->getLeafsNeighborsWithIndex(neighbors, neighborsIndex, currentIter.getCurrentGlobalIndex(),LeafIndex);
                    myThreadkernels->P2P( currentIter.getCurrentGlobalIndex(), currentIter.getCurrentListTargets(), currentIter.getCurrentListSrc() , neighbors, neighborsIndex, counter);
                }
            }
        }
        FDEBUG(computationCounter.tac());

        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation L2P + P2P : " << computationCounter.elapsed() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

};






#endif //FFMMALGORITHMTHREAD_HPP

// [--LICENSE--]
