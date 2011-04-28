#ifndef FFMMALGORITHMTHREADPROC_HPP
#define FFMMALGORITHMTHREADPROC_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"
#include "../Utils/FBoolArray.hpp"

#include "../Containers/FOctree.hpp"


//================================================================================================
#ifdef SCALFMM_USE_MPI
// Compile by mpic++ testApplication.cpp ../Utils/FAssertable.cpp -o testApplication.exe
// run by mpirun -np 4 ./testApplication.exe
#include "../Utils/FMpiApplication.hpp"
typedef FMpiApplication ApplicationImplementation;
#else
// Compile by g++ testApplication.cpp ../Utils/FAssertable.cpp -o testApplication.exe
#include "../Utils/FSingleApplication.hpp"
typedef FSingleApplication ApplicationImplementation;
#endif
//================================================================================================

#include <omp.h>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmThreadProc
* @brief
* Please read the license
*
* This class is a threaded FMM algorithm
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
class FFmmAlgorithmThreadProc : protected FAssertable, protected ApplicationImplementation{
    // To reduce the size of variable type based on foctree in this file
    typedef FOctree<ParticleClass, CellClass, LeafClass, OctreeHeight, SubtreeHeight> Octree;
    typedef typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator OctreeIterator;
    typedef KernelClass<ParticleClass, CellClass, OctreeHeight> Kernel;

    Octree* const tree;                  //< The octree to work on
    Kernel* kernels[FThreadNumbers];          //< The kernels

    FDEBUG(FTic counterTime);                //< In case of debug: to count the elapsed time
    FDEBUG(FTic computationCounter);     //< In case of debug: to  count computation time

    OctreeIterator* iterArray;
    OctreeIterator* previousIterArray;

    int previousLeft;
    int previousRight;
    int previousSize;

    int leftOffsets[OctreeHeight];
    int rightOffsets[OctreeHeight];

    void run(){}

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
    FFmmAlgorithmThreadProc(Octree* const inTree, Kernel* const inKernels, const int inArgc, char ** const inArgv )
                      : ApplicationImplementation(inArgc,inArgv), tree(inTree) , iterArray(0),
                        previousIterArray(0), previousLeft(0),previousRight(0), previousSize(0) {

        assert(tree, "tree cannot be null", __LINE__, __FILE__);
        assert(kernels, "kernels cannot be null", __LINE__, __FILE__);

        for(int idxThread = 0 ; idxThread < FThreadNumbers ; ++idxThread){
            this->kernels[idxThread] = new KernelClass<ParticleClass, CellClass, OctreeHeight>(*inKernels);
        }

        FDEBUG(FDebug::Controller << "FFmmAlgorithmThreadProc\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmThreadProc(){
        for(int idxThread = 0 ; idxThread < FThreadNumbers ; ++idxThread){
            delete this->kernels[idxThread];
        }
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

        // init kernels
        for(int idxThread = 0 ; idxThread < FThreadNumbers ; ++idxThread){
            this->kernels[idxThread]->init();
        }

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

    /** P2M */
    void bottomPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG( counterTime.tic() );

        OctreeIterator octreeIterator(tree);
        const int nbProcess = processCount();
        const int idPorcess = processId();
        int leafs = 0;

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            iterArray[leafs] = octreeIterator;
            ++leafs;
        } while(octreeIterator.moveRight());

        const float stepIdx = float(leafs) / nbProcess;
        const int startIdx = idPorcess*stepIdx;
        const int endIdx = (idPorcess+1)*stepIdx;

        this->previousLeft = startIdx;
        this->previousRight = endIdx - 1;
        this->previousSize = leafs;

        FDEBUG(computationCounter.tic());
        #pragma omp parallel num_threads(FThreadNumbers)
        {
            Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for
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

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

        // Start from leal level - 1
        OctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();
        OctreeIterator avoidGotoLeftIterator(octreeIterator);

        const int nbProcess = processCount();
        const int idPorcess = processId();

        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
            //print();

            int leafs = 0;
            // for each cells
            do{
                iterArray[leafs] = octreeIterator;
                ++leafs;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();

            const float stepIdx = float(leafs) / nbProcess;
            const int startIdx = idPorcess*stepIdx;
            const int endIdx = (idPorcess+1)*stepIdx;

            //std::cout << idPorcess << ">>--startIdx " << (startIdx) << " endIdx " << (endIdx) << std::endl;
            //std::cout << idPorcess << ">>--previousLeft " << (previousLeft) << " previousRight " << (previousRight) << std::endl;

            int leftOffset = 0;
            {
                const MortonIndex MostLeftChild = iterArray[startIdx].getCurrentGlobalIndex() << 3;
                const MortonIndex leftChildIter = previousIterArray[previousLeft].getCurrentGlobalIndex();
                //std::cout << idPorcess << ">>--MostLeftChild " << (MostLeftChild) << " leftChildIter " << (leftChildIter) << std::endl;
                if(leftChildIter < MostLeftChild){
                    while( previousIterArray[previousLeft+leftOffset].getCurrentGlobalIndex() < MostLeftChild){
                        ++leftOffset;
                    }
                }
                else if(this->previousLeft > 0 && leftChildIter > MostLeftChild){
                    while( previousIterArray[previousLeft+leftOffset - 1].getCurrentGlobalIndex() >= MostLeftChild){
                        --leftOffset;
                    }
                }
            }
            //std::cout << idPorcess << ">>--leftOffset " << (leftOffset) << std::endl;

            int rightOffset = 0;
            {
                const MortonIndex MostRightChild = (iterArray[endIdx-1].getCurrentGlobalIndex() << 3) | 7;
                const MortonIndex rightChildIter = previousIterArray[previousRight].getCurrentGlobalIndex();
                //std::cout << idPorcess << ">>--MostRightChild " << (MostRightChild) << " rightChildIter " << (rightChildIter) << std::endl;
                if(this->previousRight < this->previousSize - 1 && rightChildIter < MostRightChild){
                    while( previousIterArray[previousRight-rightOffset+1].getCurrentGlobalIndex() <= MostRightChild){
                        --rightOffset;
                    }
                }
                else if(rightChildIter > MostRightChild){
                    while( previousIterArray[previousRight-rightOffset].getCurrentGlobalIndex() > MostRightChild){
                        ++rightOffset;
                    }
                }
            }
            //std::cout << idPorcess << ">>--rightOffset " << (rightOffset) << std::endl;

            leftOffsets[idxLevel+1] = leftOffset;
            rightOffsets[idxLevel+1] = rightOffset;

            #pragma omp parallel num_threads(FThreadNumbers)
            {
                // send computed data nowait
                #pragma omp single
                {
                    //std::cout << idPorcess << ">>--Will send " << (leftOffset) << " and " << (rightOffset) << std::endl;
                    for(int idxLeafs = 0 ; idxLeafs < leftOffset ; ++idxLeafs){
                        const int idxReceiver = ((previousLeft+idxLeafs)/int(this->previousSize/nbProcess)) - 1;
                        sendData(idxReceiver,sizeof(CellClass),previousIterArray[previousLeft+idxLeafs].getCurrentCell(),previousLeft+idxLeafs);
                        //std::cout << idPorcess << "\t>>-- sends left to " << (idxReceiver) << " index " << (previousLeft+idxLeafs) << std::endl;
                    }
                    for(int idxLeafs = 0 ; idxLeafs < rightOffset ; ++idxLeafs){
                        const int idxReceiver = FMath::Min(((previousRight+idxLeafs)/int(this->previousSize/nbProcess)) + 1, nbProcess - 1);
                        sendData(idxReceiver,sizeof(CellClass),previousIterArray[previousRight+idxLeafs].getCurrentCell(),previousRight+idxLeafs);
                        //std::cout << idPorcess << "\t>>-- sends right to " << (idxReceiver) << " index " << (previousRight+idxLeafs) << std::endl;
                    }
                    //std::cout << idPorcess <<  ">>--All sent--" << std::endl;
                }

                // received computed data
                #pragma omp single
                {
                    int needToReceive = FMath::Max(0,-rightOffset) + FMath::Max(0,-leftOffset);
                    CellClass tempCell;
                    int source = 0, tag = 0, filled = 0;

                    //std::cout << idPorcess <<  ">>--Will receive " << needToReceive << std::endl;

                    while(needToReceive){
                        receiveData(sizeof(CellClass),&tempCell,&source,&tag,&filled);
                        if(filled){
                            *previousIterArray[tag].getCurrentCell() = tempCell;
                        }
                        --needToReceive;
                        //std::cout << idPorcess <<  ">>receive tag " << (tag) << " tempCell.up " << tempCell.getDataUp() << std::endl;
                    }
                    //std::cout << idPorcess <<  ">>--All receive--" << std::endl;
                }

                #pragma omp single nowait
                {
                    FDEBUG(computationCounter.tic());
                }

                Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                #pragma omp for
                for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                    myThreadkernels->M2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentChild(), idxLevel);
                }

                #pragma omp single nowait
                {
                    FDEBUG(computationCounter.tac());
                    FDEBUG(totalComputation += computationCounter.elapsed());
                }

            }

            swapArray();
            this->previousLeft = startIdx;
            this->previousRight = endIdx - 1;
            this->previousSize = leafs;
        }

        //print();

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** M2L L2L */
    void downardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

        { // first M2L
            OctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();
            OctreeIterator avoidGotoLeftIterator(octreeIterator);

            const int nbProcess = processCount();
            const int idPorcess = processId();


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

                const float stepIdx = float(leafs) / nbProcess;
                const int startIdx = idPorcess*stepIdx;
                const int endIdx = (idPorcess+1)*stepIdx;

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

                //std::cout << "There are " << leafs << " leafs" << std::endl;

                print();

                #pragma omp parallel num_threads(FThreadNumbers)
                {
                    CellClass* neighbors[208];
                    MortonIndex neighborsIndexes[208];
                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];

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
                                const int idxReceiver = FMath::Min( cellPositionInArray / int( (leafs+1) / nbProcess) , nbProcess - 1 );
                                #pragma omp critical(CheckToSend)
                                {
                                    if(!alreadySent[idxReceiver]->get(idxLeafs)){
                                        std::cout << idPorcess << ">>--idxLeafs " << (idxLeafs) << " idxReceiver " << (idxReceiver)
                                                << " cellPositionInArray " << (cellPositionInArray)  << " indexCell " << indexCell<< std::endl;
                                        sendData(idxReceiver,sizeof(CellClass),iterArray[idxLeafs].getCurrentCell(),idxLeafs);
                                        alreadySent[idxReceiver]->set(idxLeafs,true);
                                        needData = true;
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
                        if(needData){
                            std::cout << idPorcess <<  ">>this cell need data " << idxLeafs << " index " << iterArray[idxLeafs].getCurrentGlobalIndex() << " neighborsCounter " << neighborsCounter << std::endl;
                            const int currentCell = idxLeafs - startIdx;
                            unfinishedCells[currentCell] = new LimitCell();
                            unfinishedCells[currentCell]->counter = neighborsCounter;
                            memcpy(unfinishedCells[currentCell]->neighbors,neighbors,sizeof(CellClass*)*neighborsCounter);
                            alreadySent[idPorcess]->set(idxLeafs,true);
                        }
                        else if(neighborsCounter){
                            std::cout << idPorcess <<  ">>compute directly " << idxLeafs << " index " << iterArray[idxLeafs].getCurrentGlobalIndex() << std::endl;
                            myThreadkernels->M2L(  iterArray[idxLeafs].getCurrentCell() , neighbors, neighborsCounter, idxLevel);
                        }
                    }

                    // received computed data
                    #pragma omp single
                    {
                        std::cout << idPorcess << ">>--needToReceive " << (needToReceive) << std::endl;

                        CellClass tempCell;
                        int source = 0, tag = 0, filled = 0;

                        while(needToReceive){
                            receiveData(sizeof(CellClass),&tempCell,&source,&tag,&filled);
                            if(filled){
                                *iterArray[tag].getCurrentCell() = tempCell;
                            }
                            --needToReceive;

                            std::cout << idPorcess <<  ">>receive tag " << (tag) << " tempCell.up " << tempCell.getDataUp() << std::endl;
                        }
                    }

                    #pragma omp for
                    for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                        if(alreadySent[idPorcess]->get(idxLeafs)){
                            std::cout << idPorcess <<  ">>finish to compute " << idxLeafs << " index " << iterArray[idxLeafs].getCurrentGlobalIndex() << std::endl;
                            myThreadkernels->M2L(  iterArray[idxLeafs].getCurrentCell() , unfinishedCells[idxLeafs-startIdx]->neighbors, unfinishedCells[idxLeafs-startIdx]->counter, idxLevel);
                            delete unfinishedCells[idxLeafs-startIdx];
                        }
                    }
                }

                delete [] unfinishedCells;

                for(int idxProc = 0 ; idxProc < nbProcess; ++idxProc){
                    delete alreadySent[idxProc];
                }

            }
        }

        //print();

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );

        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( totalComputation = 0 );
        { // second L2L
            OctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();

            OctreeIterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = OctreeHeight - 1;
            const int nbProcess = processCount();
            const int idPorcess = processId();

            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                //print();

                // keep data
                /*swapArray();
                this->previousLeft = startIdx;
                this->previousRight = endIdx - 1;
                this->previousSize = leafs;*/

                 int leafs = 0;
                 // for each cells
                 do{
                     iterArray[leafs] = octreeIterator;
                     ++leafs;
                 } while(octreeIterator.moveRight());
                 avoidGotoLeftIterator.moveDown();
                 octreeIterator = avoidGotoLeftIterator;

                 const float stepIdx = float(leafs) / nbProcess;
                 const int startIdx = idPorcess*stepIdx;
                 const int endIdx = (idPorcess+1)*stepIdx;

                 const int currentLeft = startIdx;
                 const int currentRight = endIdx -1;

                #pragma omp parallel num_threads(FThreadNumbers)
                {
                    // send computed data
                    #pragma omp single //nowait
                    {
                        const int leftOffset = -leftOffsets[idxLevel];
                        for(int idxLeafs = 1 ; idxLeafs <= leftOffset ; ++idxLeafs){
                            const int idxReceiver = (currentLeft-idxLeafs)/int(leafs/nbProcess);
                            sendData(idxReceiver,sizeof(CellClass),iterArray[currentLeft-idxLeafs].getCurrentCell(),currentLeft-idxLeafs);
                            //std::cout << idPorcess << "\t>>-- sends (1) to " << (idxReceiver) << " index " << (currentLeft-idxLeafs) << std::endl;
                        }
                        const int rightOffset = -rightOffsets[idxLevel];
                        for(int idxLeafs = 1 ; idxLeafs <= rightOffset ; ++idxLeafs){
                            const int idxReceiver = FMath::Min((currentRight+idxLeafs)/int(leafs/nbProcess),nbProcess-1);
                            sendData(idxReceiver,sizeof(CellClass),iterArray[currentRight+idxLeafs].getCurrentCell(),currentRight+idxLeafs);
                            //std::cout << idPorcess << "\t>>-- sends (2) to " << (idxReceiver) << " index " << (currentRight+idxLeafs) << " currentRight " << currentRight << std::endl;
                        }
                        //std::cout << idPorcess << ">>--Will send " << (leftOffset) << " and " << (rightOffset) << std::endl;
                    }

                    // received computed data
                    #pragma omp single
                    {
                        int needToReceive = FMath::Max(0,rightOffsets[idxLevel]) + FMath::Max(0,leftOffsets[idxLevel]);
                        CellClass tempCell;
                        int source = 0, tag = 0, filled = 0;

                        //std::cout << idPorcess << ">>--needToReceive " << (needToReceive) << std::endl;

                        while(needToReceive){
                            receiveData(sizeof(CellClass),&tempCell,&source,&tag,&filled);
                            if(filled){
                                iterArray[tag].getCurrentCell()->addCell(tempCell);
                            }
                            --needToReceive;
                        }
                        //std::cout << "all received" << std::endl;
                    }
                }

                #pragma omp parallel num_threads(FThreadNumbers)
                {
                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                    #pragma omp for
                    for(int idxLeafs = startIdx ; idxLeafs < endIdx ; ++idxLeafs){
                        myThreadkernels->L2L( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentChild(), idxLevel);
                    }
                }

            }
        }

        //print();

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** P2P */
    void directPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );

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

        FDEBUG(computationCounter.tic());

        const float stepIdx = float(leafs) / processCount();
        const int startIdx = processId()*stepIdx;
        const int endIdx = (processId()+1)*stepIdx;

        #pragma omp parallel num_threads(FThreadNumbers)
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

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** This function test the octree to be sure that the fmm algorithm
      * has worked completly.
      */
    void ValidateFMMAlgoProc(){
        std::cout << "Check Result\n";
        int NbPart = 0;
        int NbLeafs = 0;
        { // Check that each particle has been summed with all other
            OctreeIterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            do{
                NbPart += octreeIterator.getCurrentListSrc()->getSize();
                ++NbLeafs;
            } while(octreeIterator.moveRight());
        }
        {
            const float stepIdx = float(NbLeafs) / processCount();
            const int startIdx = processId()*stepIdx;
            const int endIdx = (processId()+1)*stepIdx;
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

        std::cout << "Done\n";
    }

    void print(){
        OctreeIterator octreeIterator(tree);
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
