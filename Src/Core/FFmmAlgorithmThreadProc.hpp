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
* schedule(runtime)
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

    int previousLeft;                   //< To store the left limit at the previous level
    int previousRight;                  //< To store the right limit at the previous level
    int previousSize;                   //< To store the size at the previous level

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
        previousIterArray(0), previousLeft(0),previousRight(0), previousSize(0),
        MaxThreads(omp_get_max_threads()), nbProcess(inApp.processCount()), idPorcess(inApp.processId()),
        sendBuffer(0) {

        assert(tree, "tree cannot be null", __LINE__, __FILE__);

        this->kernels = new Kernel*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass<ParticleClass, CellClass, OctreeHeight>(*inKernels);
        }

        this->sendBuffer = new FBufferVector<2000>[nbProcess];

        FDEBUG(FDebug::Controller << "FFmmAlgorithmThreadProc\n");
        FDEBUG(FDebug::Controller << "Max threads = "  << MaxThreads << " .\n");
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
        this->previousLeft = startIdx;
        this->previousRight = endIdx - 1;
        this->previousSize = leafs;


        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.tacAndElapsed() << "s)\n" );
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

        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){

            int leafs = 0;
            // for each cells
            do{
                iterArray[leafs] = octreeIterator;
                ++leafs;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();

            const int startIdx = getLeft(leafs);
            const int endIdx = getRight(leafs);

            if(startIdx < leafs){
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

                            const int idxReceiver = getProc(parentOffset,leafs);
                            if(!this->sendBuffer[idxReceiver].addEnoughSpace(previousIterArray[this->previousLeft+leftOffset].getCurrentCell()->bytesToSendUp())){
                                app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                                this->sendBuffer[idxReceiver].clear();
                            }
                            this->sendBuffer[idxReceiver].addDataUp(previousLeft+leftOffset,*previousIterArray[this->previousLeft+leftOffset].getCurrentCell());

                            ++leftOffset;
                        }
                        for(int idxProc = 0 ; idxProc < this->nbProcess ; ++idxProc){
                            if(this->sendBuffer[idxProc].getSize()){
                                app.sendData(idxProc,this->sendBuffer[idxProc].getSize(),this->sendBuffer[idxProc].getData(),idxLevel);
                                this->sendBuffer[idxProc].clear();
                            }
                        }
                    }
                    else if(this->previousLeft > 0 && leftChildIter > MostLeftChild){
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

                    if(this->previousRight < this->previousSize - 1 && rightChildIter < MostRightChild){
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
                            const int idxReceiver = getProc(parentOffset,leafs);
                            if(!this->sendBuffer[idxReceiver].addEnoughSpace(previousIterArray[this->previousRight-rightOffset].getCurrentCell()->bytesToSendUp())){
                                app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                                this->sendBuffer[idxReceiver].clear();
                            }
                            this->sendBuffer[idxReceiver].addDataUp(previousRight-rightOffset,*previousIterArray[this->previousRight-rightOffset].getCurrentCell());

                            ++rightOffset;
                        }
                        for(int idxProc = 0 ; idxProc < this->nbProcess ; ++idxProc){
                            if(this->sendBuffer[idxProc].getSize()){
                                app.sendData(idxProc,this->sendBuffer[idxProc].getSize(),this->sendBuffer[idxProc].getData(),idxLevel);
                                this->sendBuffer[idxProc].clear();
                            }
                        }
                    }
                    rightOffsets[idxLevel+1] = rightOffset;
                }
                FDEBUG(sendCounter.tac());

                #pragma omp parallel
                {
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
                    for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                        myThreadkernels->M2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentChild(), idxLevel);
                    }

                    #pragma omp single nowait
                    {
                        FDEBUG(computationCounter.tac());
                    }

                }
            }
            else {
                int parentOffset = leafs - 1;
                MortonIndex parentIndex = iterArray[parentOffset].getCurrentGlobalIndex();

                for(int idxLeafs = previousRight ; idxLeafs >= previousLeft ; --idxLeafs){
                    const MortonIndex childIndex = previousIterArray[idxLeafs].getCurrentGlobalIndex() >> 3;
                    while(childIndex != parentIndex){
                        --parentOffset;
                        parentIndex = iterArray[parentOffset].getCurrentGlobalIndex();
                    }
                    const int idxReceiver = getProc(parentOffset,leafs);
                    if(!this->sendBuffer[idxReceiver].addEnoughSpace(previousIterArray[idxLeafs].getCurrentCell()->bytesToSendUp())){
                        app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                        this->sendBuffer[idxReceiver].clear();
                    }
                    this->sendBuffer[idxReceiver].addDataUp(idxLeafs,*previousIterArray[idxLeafs].getCurrentCell());
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
            this->previousLeft = startIdx;
            this->previousRight = endIdx - 1;
            this->previousSize = leafs;

        }

        app.processBarrier();

        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.tacAndElapsed() << "s)\n" );
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

            OctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();
            OctreeIterator avoidGotoLeftIterator(octreeIterator);

            // for each levels
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                int leafs = 0;
                // for each cells
                do{
                    iterArray[leafs] = octreeIterator;
                    ++leafs;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                const int startIdx = getLeft(leafs);
                const int endIdx = getRight(leafs);

                struct LimitCell {
                    int counter;
                    CellClass* neighbors[208];
                };
                LimitCell ** unfinishedCells = new LimitCell*[endIdx - startIdx];

                // Get limit in term of morton index
                const MortonIndex startIdxIndex = iterArray[startIdx].getCurrentGlobalIndex();
                const MortonIndex endIdxIndex = iterArray[endIdx-1].getCurrentGlobalIndex();

                int needToReceive = 0;
                FBoolArray* alreadySent[nbProcess];
                for(int idxProc = 0 ; idxProc < nbProcess; ++idxProc){
                    alreadySent[idxProc] = new FBoolArray(leafs);
                }
                #pragma omp parallel
                {
                    CellClass* neighbors[208];
                    MortonIndex neighborsIndexes[208];
                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];

                    #pragma omp single nowait
                    {
                        FDEBUG(sendCounter.tic());
                    }

                    #pragma omp for
                    for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                        const int neighborsCounter = tree->getDistantNeighborsWithIndex(neighbors, neighborsIndexes, iterArray[idxLeafs].getCurrentGlobalIndex(),idxLevel);
                        bool needData = false;

                        for(int idxNeighbor = 0 ; idxNeighbor < neighborsCounter ; ++idxNeighbor){
                            // Get Morton index from this neigbor cell
                            const MortonIndex indexCell = neighborsIndexes[idxNeighbor];
                            // Is this cell out of our interval?
                            if(indexCell < startIdxIndex || endIdxIndex < indexCell){
                                // Yes we need to send the center of computation
                                // but who this cell belong to?
                                int cellPositionInArray = 0;
                                const CellClass* const cell = neighbors[idxNeighbor];
                                if(indexCell < startIdxIndex){
                                    // This cell is on the left
                                    for(int idxFind = startIdx - 1; idxFind >= 0 ; --idxFind){
                                        if(iterArray[idxFind].getCurrentCell() == cell){
                                            cellPositionInArray = idxFind;
                                            break;
                                        }
                                    }
                                }
                                else {
                                    // This cell is on the right
                                    for(int idxFind = endIdx ; idxFind < leafs ; ++idxFind){
                                        if(iterArray[idxFind].getCurrentCell() == cell){
                                            cellPositionInArray = idxFind;
                                            break;
                                        }
                                    }
                                }

                                // Find receiver and send him the cell
                                const int idxReceiver = getProc(cellPositionInArray,leafs);
                                #pragma omp critical(CheckToSend)
                                {
                                    if(!alreadySent[idxReceiver]->get(idxLeafs)){
                                        alreadySent[idxReceiver]->set(idxLeafs,true);
                                        needData = true;

                                        if(!this->sendBuffer[idxReceiver].addEnoughSpace(iterArray[idxLeafs].getCurrentCell()->bytesToSendUp())){
                                            app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                                            this->sendBuffer[idxReceiver].clear();
                                        }
                                        this->sendBuffer[idxReceiver].addDataUp(idxLeafs,*iterArray[idxLeafs].getCurrentCell());
                                    }
                                }
                                #pragma omp critical(CheckToReceive)
                                {
                                    if(!alreadySent[idPorcess]->get(cellPositionInArray)){
                                        ++needToReceive;
                                        alreadySent[idPorcess]->set(cellPositionInArray,true);
                                        needData = true;
                                    }
                                }
                            }
                        }
                        // if need data we can not compute now
                        if(needData){
                            const int currentCell = idxLeafs - startIdx;
                            unfinishedCells[currentCell] = new LimitCell();
                            unfinishedCells[currentCell]->counter = neighborsCounter;
                            memcpy(unfinishedCells[currentCell]->neighbors,neighbors,sizeof(CellClass*)*neighborsCounter);
                            alreadySent[idPorcess]->set(idxLeafs,true);
                        }
                        // we can compute now !
                        else if(neighborsCounter){
                            myThreadkernels->M2L(  iterArray[idxLeafs].getCurrentCell() , neighbors, neighborsCounter, idxLevel);
                        }
                    }

                    #pragma omp single nowait
                    {
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
                        FDEBUG( FDebug::Controller << "\t\tNeed to receive "  << needToReceive << " cells.\n" );
                        FDEBUG(receiveCounter.tic());

                        int source = 0, filled = 0;
                        int position;
                        char buffer[BufferSize];

                        while(needToReceive){
                            app.receiveData( BufferSize, idxLevel, buffer, &source, &filled);
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
                    for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                        if(alreadySent[idPorcess]->get(idxLeafs)){
                            myThreadkernels->M2L(  iterArray[idxLeafs].getCurrentCell() , unfinishedCells[idxLeafs-startIdx]->neighbors, unfinishedCells[idxLeafs-startIdx]->counter, idxLevel);
                            delete unfinishedCells[idxLeafs-startIdx];
                        }
                    }

                    #pragma omp single nowait
                    {
                        FDEBUG(computationCounter.tac());
                    }
                }

                delete [] unfinishedCells;

                for(int idxProc = 0 ; idxProc < nbProcess; ++idxProc){
                    delete alreadySent[idxProc];
                }

            }
            app.processBarrier();

            FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Send : " << sendCounter.cumulated() << " s\n" );
            FDEBUG( FDebug::Controller << "\t\t Receive : " << receiveCounter.cumulated() << " s\n" );
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

                int leafs = 0;
                // for each cells
                do{
                    iterArray[leafs] = octreeIterator;
                    ++leafs;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                const int startIdx = getLeft(leafs);
                const int endIdx = getRight(leafs);

                const int currentLeft = startIdx;
                const int currentRight = endIdx -1;

                #pragma omp parallel
                {
                    // send computed data
                    #pragma omp single nowait
                    {
                        FDEBUG(sendCounter.tic());
                        const int leftOffset = -leftOffsets[idxLevel];
                        for(int idxLeafs = 1 ; idxLeafs <= leftOffset ; ++idxLeafs){
                            const int idxReceiver = getProc((currentLeft-idxLeafs),leafs);
                            if(!this->sendBuffer[idxReceiver].addEnoughSpace(iterArray[currentLeft-idxLeafs].getCurrentCell()->bytesToSendDown())){
                                app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                                this->sendBuffer[idxReceiver].clear();
                            }
                            this->sendBuffer[idxReceiver].addDataDown(currentLeft-idxLeafs,*iterArray[currentLeft-idxLeafs].getCurrentCell());
                        }
                        const int rightOffset = -rightOffsets[idxLevel];
                        for(int idxLeafs = 1 ; idxLeafs <= rightOffset ; ++idxLeafs){
                            const int idxReceiver = getProc((currentRight+idxLeafs),leafs);
                            if(!this->sendBuffer[idxReceiver].addEnoughSpace(iterArray[currentRight+idxLeafs].getCurrentCell()->bytesToSendDown())){
                                app.sendData(idxReceiver,this->sendBuffer[idxReceiver].getSize(),this->sendBuffer[idxReceiver].getData(),idxLevel);
                                this->sendBuffer[idxReceiver].clear();
                            }
                            this->sendBuffer[idxReceiver].addDataDown(currentRight+idxLeafs,*iterArray[currentRight+idxLeafs].getCurrentCell());
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
                        for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                            myThreadkernels->L2L( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentChild(), idxLevel);
                        }
                    }
                    FDEBUG(computationCounter.tac());
                }
            }

            app.processBarrier();

            FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.tacAndElapsed() << "s)\n" );
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

        const int LeafIndex = OctreeHeight - 1;
        int leafs = 0;
        {
            OctreeIterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            // for each leafs
            do{
                iterArray[leafs] = octreeIterator;
                ++leafs;
            } while(octreeIterator.moveRight());
        }

        FDEBUG(FTic computationCounter);

        const int startIdx = getLeft(leafs);
        const int endIdx = getRight(leafs);

        #pragma omp parallel
        {
            Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
            // There is a maximum of 26 neighbors
            FList<ParticleClass*>* neighbors[26];

            #pragma omp for
            for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                myThreadkernels->L2P(iterArray[idxLeafs].getCurrentCell(), iterArray[idxLeafs].getCurrentListTargets());
                // need the current particles and neighbors particles
                const int counter = tree->getLeafsNeighbors(neighbors, iterArray[idxLeafs].getCurrentGlobalIndex(),LeafIndex);
                myThreadkernels->P2P( iterArray[idxLeafs].getCurrentListTargets(), iterArray[idxLeafs].getCurrentListSrc() , neighbors, counter);
            }
        }
        FDEBUG(computationCounter.tac());


        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Test function
    /////////////////////////////////////////////////////////////////////////////

    /** This function test the octree to be sure that the fmm algorithm
      * has worked completly.
      */
    void ValidateFMMAlgoProc(Octree* const valideTree){
        std::cout << "Check Result\n";
        {
            OctreeIterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();

            OctreeIterator octreeIteratorValide(valideTree);
            octreeIteratorValide.gotoBottomLeft();

            for(int level = OctreeHeight - 1 ; level > 0 ; --level){
                int NbLeafs = 0;
                do{
                    ++NbLeafs;
                } while(octreeIterator.moveRight());
                octreeIterator.gotoLeft();

                const int startIdx = getLeft(NbLeafs);
                const int endIdx = getRight(NbLeafs);
                // Check that each particle has been summed with all other

                for(int idx = 0 ; idx < startIdx ; ++idx){
                    octreeIterator.moveRight();
                    octreeIteratorValide.moveRight();
                }

                for(int idx = startIdx ; idx < endIdx ; ++idx){
                    if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                        std::cout << "Error index are not equal!" << std::endl;
                    }
                    else{
                        if(octreeIterator.getCurrentCell()->getDataUp() != octreeIteratorValide.getCurrentCell()->getDataUp()){
                            std::cout << "M2M error at level " << level << " up bad " << octreeIterator.getCurrentCell()->getDataUp()
                                    << " good " << octreeIteratorValide.getCurrentCell()->getDataUp() << " idx " << idx << std::endl;
                        }
                        if(octreeIterator.getCurrentCell()->getDataDown() != octreeIteratorValide.getCurrentCell()->getDataDown()){
                            std::cout << "L2L error at level " << level << " down bad " << octreeIterator.getCurrentCell()->getDataDown()
                                    << " good " << octreeIteratorValide.getCurrentCell()->getDataDown() << " idx " << idx << std::endl;
                        }
                    }

                    octreeIterator.moveRight();
                    octreeIteratorValide.moveRight();
                }

                octreeIterator.moveUp();
                octreeIterator.gotoLeft();

                octreeIteratorValide.moveUp();
                octreeIteratorValide.gotoLeft();
            }
        }
        {
            int NbPart = 0;
            int NbLeafs = 0;
            { // Check that each particle has been summed with all other
                OctreeIterator octreeIterator(tree);
                octreeIterator.gotoBottomLeft();
                do{
                    NbPart += octreeIterator.getCurrentListSrc()->getSize();
                    ++NbLeafs;
                } while(octreeIterator.moveRight());
                std::cout << "There is " << NbPart << " particles on " << NbLeafs << " Leafs" << std::endl;
            }
            {
                const int startIdx = getLeft(NbLeafs);
                const int endIdx = getRight(NbLeafs);
                // Check that each particle has been summed with all other
                OctreeIterator octreeIterator(tree);
                octreeIterator.gotoBottomLeft();

                for(int idx = 0 ; idx < startIdx ; ++idx){
                    octreeIterator.moveRight();
                }

                for(int idx = startIdx ; idx < endIdx ; ++idx){
                    typename FList<ParticleClass*>::BasicIterator iter(*octreeIterator.getCurrentListTargets());

                    const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());

                    while( iter.isValide() ){
                        // If a particles has been impacted by less than NbPart - 1 (the current particle)
                        // there is a problem
                        if( (!isUsingTsm && iter.value()->getDataDown() != NbPart - 1) ||
                            (isUsingTsm && iter.value()->getDataDown() != NbPart) ){
                            std::cout << "Problem L2P + P2P, value on particle is : " << iter.value()->getDataDown() << "\n";
                        }
                        iter.progress();
                    }
                    octreeIterator.moveRight();
                }
            }
        }

        std::cout << "Done\n";
    }

    /** To print an octree
      * used to debug and understand how the values were passed
      */
    void print(Octree* const valideTree){
        OctreeIterator octreeIterator(valideTree);
        for(int idxLevel = OctreeHeight - 1 ; idxLevel > 1 ; --idxLevel ){
            do{
                std::cout << "[" << octreeIterator.getCurrentGlobalIndex() << "] up:" << octreeIterator.getCurrentCell()->getDataUp() << " down:" << octreeIterator.getCurrentCell()->getDataDown() << "\t";
            } while(octreeIterator.moveRight());
            std::cout << "\n";
            octreeIterator.gotoLeft();
            octreeIterator.moveDown();
        }
    }

};






#endif //FFMMALGORITHMTHREAD_HPP

// [--LICENSE--]
