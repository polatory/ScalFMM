#ifndef FFMMALGORITHMTHREAD_HPP
#define FFMMALGORITHMTHREAD_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"

#include "../Containers/FOctree.hpp"


#include <omp.h>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmThread
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
* When using this algorithm the P2P is thread safe.
*/
template<template< class ParticleClass, class CellClass> class KernelClass,
        class ParticleClass, class CellClass,
        template<class ParticleClass> class LeafClass>
class FFmmAlgorithmThread : protected FAssertable{
    // To reduce the size of variable type based on foctree in this file
    typedef FOctree<ParticleClass, CellClass, LeafClass> Octree;
    typedef typename FOctree<ParticleClass, CellClass,LeafClass>::Iterator OctreeIterator;
    typedef KernelClass<ParticleClass, CellClass> Kernel;

    Octree* const tree;                  //< The octree to work on
    Kernel** kernels;                    //< The kernels

    OctreeIterator* iterArray;
    int leafsNumber;

    static const int SizeShape = 3*3*3;
    int shapeLeaf[SizeShape];    

    const int MaxThreads;

    const int OctreeHeight;

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThread(Octree* const inTree, Kernel* const inKernels)
                      : tree(inTree) , kernels(0), iterArray(0), leafsNumber(0),
                        MaxThreads(omp_get_max_threads()), OctreeHeight(tree->getHeight()) {

        assert(tree, "tree cannot be null", __LINE__, __FILE__);

        this->kernels = new Kernel*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass<ParticleClass, CellClass>(*inKernels);
        }

        FDEBUG(FDebug::Controller << "FFmmAlgorithmThread (Max Thread " << omp_get_max_threads() << ")\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmThread(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );

        for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
            this->shapeLeaf[idxShape] = 0;
        }
        //const int LeafIndex = OctreeHeight - 1;

        // Count leaf
        leafsNumber = 0;
        OctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            ++leafsNumber;
            const FTreeCoordinate& coord = octreeIterator.getCurrentCell()->getCoordinate();
            //const MortonIndex index = octreeIterator.getCurrentGlobalIndex();
            //coord.setPositionFromMorton(index, LeafIndex);
            ++this->shapeLeaf[(coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3)];

        } while(octreeIterator.moveRight());
        iterArray = new OctreeIterator[leafsNumber];
        assert(iterArray, "iterArray bad alloc", __LINE__, __FILE__);

        bottomPass();
        upwardPass();

        downardPass();

        directPass();

        delete [] iterArray;
        iterArray = 0;

        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
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

        FDEBUG(FTic computationCounter);
        #pragma omp parallel
        {
            Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for nowait
            for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                // We need the current cell that represent the leaf
                // and the list of particles
                myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentListSrc());
            }
        }
        FDEBUG(computationCounter.tac() );

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

        // Start from leal level - 1
        OctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();
        OctreeIterator avoidGotoLeftIterator(octreeIterator);

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

            FDEBUG(computationCounter.tic());
            #pragma omp parallel
            {
                Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                #pragma omp for
                for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                    // We need the current cell and the child
                    // child is an array (of 8 child) that may be null
                    myThreadkernels->M2M( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
                }
            }

            FDEBUG(computationCounter.tac());
        }


        FDEBUG( FDebug::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
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

            OctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();
            OctreeIterator avoidGotoLeftIterator(octreeIterator);

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

                FDEBUG(computationCounter.tic());
                #pragma omp parallel
                {
                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                    const CellClass* neighbors[208];

                    #pragma omp for
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        const int counter = tree->getDistantNeighbors(neighbors,  iterArray[idxCell].getCurrentGlobalCoordinate(),idxLevel);
                        if(counter) myThreadkernels->M2L( iterArray[idxCell].getCurrentCell() , neighbors, counter, idxLevel);
                    }
                }
                FDEBUG(computationCounter.tac());
            }

            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        }

        { // second L2L
            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            FDEBUG(FTic computationCounter);

            OctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();

            OctreeIterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = OctreeHeight - 1;
            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                int numberOfCells = 0;
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
                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                    #pragma omp for
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        myThreadkernels->L2L( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
                    }
                }
                FDEBUG(computationCounter.tac());
            }

            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
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
        FDEBUG(FTic computationCounter);
        FDEBUG(FTic computationCounterP2P);

        int countShape[SizeShape];
        for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
            countShape[idxShape] = 0;
        }

        struct LeafData{
            MortonIndex index;
            CellClass* cell;
            FList<ParticleClass>* targets;
            FList<ParticleClass>* sources;
        };
        LeafData* const leafsDataArray = new LeafData[this->leafsNumber];

        const int LeafIndex = OctreeHeight - 1;
        int leafs = 0;
        {
            int startPosAtShape[SizeShape];
            startPosAtShape[0] = 0;
            for(int idxShape = 1 ; idxShape < SizeShape ; ++idxShape){
                startPosAtShape[idxShape] = startPosAtShape[idxShape-1] + this->shapeLeaf[idxShape-1];
            }

            OctreeIterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            // for each leafs
            do{
                //iterArray[leafs] = octreeIterator;
                ++leafs;
                const MortonIndex index = octreeIterator.getCurrentGlobalIndex();
                FTreeCoordinate coord;
                coord.setPositionFromMorton(index, LeafIndex);
                const int shapePosition = (coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3);

                leafsDataArray[startPosAtShape[shapePosition]].index = octreeIterator.getCurrentGlobalIndex();
                leafsDataArray[startPosAtShape[shapePosition]].cell = octreeIterator.getCurrentCell();
                leafsDataArray[startPosAtShape[shapePosition]].targets = octreeIterator.getCurrentListTargets();
                leafsDataArray[startPosAtShape[shapePosition]].sources = octreeIterator.getCurrentListSrc();

                ++startPosAtShape[shapePosition];

            } while(octreeIterator.moveRight());
        }

        FDEBUG(computationCounter.tic());

        #pragma omp parallel
        {
            Kernel& myThreadkernels = (*kernels[omp_get_thread_num()]);
            // There is a maximum of 26 neighbors
            FList<ParticleClass>* neighbors[26];
            MortonIndex neighborsIndex[26];
            int previous = 0;

            for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
                const int endAtThisShape = this->shapeLeaf[idxShape] + previous;

                #pragma omp for
                for(int idxLeafs = previous ; idxLeafs < endAtThisShape ; ++idxLeafs){
                    LeafData& currentIter = leafsDataArray[idxLeafs];
                    myThreadkernels.L2P(currentIter.cell, currentIter.targets);
                    // need the current particles and neighbors particles
                    FDEBUG(if(!omp_get_thread_num()) computationCounterP2P.tic());
                    const int counter = tree->getLeafsNeighborsWithIndex(neighbors, neighborsIndex, currentIter.index, LeafIndex);
                    myThreadkernels.P2P(currentIter.index, currentIter.targets, currentIter.sources , neighbors, neighborsIndex, counter);
                    FDEBUG(if(!omp_get_thread_num()) computationCounterP2P.tac());
                }

                previous = endAtThisShape;
            }
        }
        FDEBUG(computationCounter.tac());

        delete [] leafsDataArray;

        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation L2P + P2P : " << computationCounter.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation P2P : " << computationCounterP2P.cumulated() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

};


#endif //FFMMALGORITHMTHREAD_HPP

// [--LICENSE--]
