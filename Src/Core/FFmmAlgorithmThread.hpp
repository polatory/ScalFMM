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
template<template< class ParticleClass, class CellClass, int OctreeHeight> class KernelClass,
        class ParticleClass, class CellClass,
        template<class ParticleClass> class LeafClass,
        int OctreeHeight, int SubtreeHeight>
class FFmmAlgorithmThread : protected FAssertable{
    // To reduce the size of variable type based on foctree in this file
    typedef FOctree<ParticleClass, CellClass, LeafClass, OctreeHeight, SubtreeHeight> Octree;
    typedef typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator OctreeIterator;
    typedef KernelClass<ParticleClass, CellClass, OctreeHeight> Kernel;

    Octree* const tree;                  //< The octree to work on
    Kernel** kernels;                    //< The kernels

    OctreeIterator* iterArray;

    static const int SizeShape = 3*3*3;
    int shapeLeaf[SizeShape];    

    const int MaxThreads;

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThread(Octree* const inTree, Kernel* const inKernels)
                      : tree(inTree) , kernels(0), iterArray(0), MaxThreads(omp_get_max_threads()) {

        assert(tree, "tree cannot be null", __LINE__, __FILE__);

        this->kernels = new Kernel*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass<ParticleClass, CellClass, OctreeHeight>(*inKernels);
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
        const int LeafIndex = OctreeHeight - 1;

        // Count leaf
        int leafs = 0;
        OctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            ++leafs;
            const MortonIndex index = octreeIterator.getCurrentGlobalIndex();
            FTreeCoordinate coord;
            coord.setPositionFromMorton(index, LeafIndex);
            ++this->shapeLeaf[(coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3)];

        } while(octreeIterator.moveRight());
        iterArray = new OctreeIterator[leafs];
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
            #pragma omp for
            for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                // We need the current cell that represent the leaf
                // and the list of particles
                myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentListSrc());
            }
        }
        FDEBUG(computationCounter.tac() );

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


        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.tacAndElapsed() << "s)\n" );
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
                    CellClass* neighbors[208];
                    #pragma omp for
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        const int counter = tree->getDistantNeighbors(neighbors,  iterArray[idxCell].getCurrentGlobalIndex(),idxLevel);
                        if(counter) myThreadkernels->M2L(  iterArray[idxCell].getCurrentCell() , neighbors, counter, idxLevel);
                    }
                }
                FDEBUG(computationCounter.tac());
            }

            FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.tacAndElapsed() << "s)\n" );
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

            FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.tacAndElapsed() << "s)\n" );
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

        OctreeIterator* shapeArray[SizeShape];
        int countShape[SizeShape];
        for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
            shapeArray[idxShape] = new OctreeIterator[this->shapeLeaf[idxShape]];
            countShape[idxShape] = 0;
        }

        const int LeafIndex = OctreeHeight - 1;
        //int leafs = 0;
        {
            OctreeIterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            // for each leafs
            do{
                //iterArray[leafs] = octreeIterator;
                //++leafs;
                const MortonIndex index = octreeIterator.getCurrentGlobalIndex();
                FTreeCoordinate coord;
                coord.setPositionFromMorton(index, LeafIndex);
                const int shapePosition = (coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3);
                shapeArray[shapePosition][countShape[shapePosition]] = octreeIterator;
                ++countShape[shapePosition];

            } while(octreeIterator.moveRight());
        }

        FDEBUG(computationCounter.tic());
        for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
            const int leafAtThisShape = this->shapeLeaf[idxShape];
            #pragma omp parallel
            {
                Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                // There is a maximum of 26 neighbors
                FList<ParticleClass*>* neighbors[26];

                #pragma omp for
                for(int idxLeafs = 0 ; idxLeafs < leafAtThisShape ; ++idxLeafs){
                    myThreadkernels->L2P(shapeArray[idxShape][idxLeafs].getCurrentCell(), shapeArray[idxShape][idxLeafs].getCurrentListTargets());
                    // need the current particles and neighbors particles
                    const int counter = tree->getLeafsNeighbors(neighbors, shapeArray[idxShape][idxLeafs].getCurrentGlobalIndex(),LeafIndex);
                    myThreadkernels->P2P( shapeArray[idxShape][idxLeafs].getCurrentListTargets(), shapeArray[idxShape][idxLeafs].getCurrentListSrc() , neighbors, counter);
                }
            }
        }
        FDEBUG(computationCounter.tac());

        for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
            delete [] shapeArray[idxShape];
        }

        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

};


#endif //FFMMALGORITHMTHREAD_HPP

// [--LICENSE--]
