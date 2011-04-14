#ifndef FFmmALGORITHMARRAY_HPP
#define FFmmALGORITHMARRAY_HPP
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
* @class FFmmAlgorithmArray
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
class FFmmAlgorithmArray : protected FAssertable{
    // To reduce the size of variable type based on foctree in this file
    typedef FOctree<ParticleClass, CellClass, LeafClass, OctreeHeight, SubtreeHeight> Octree;
    typedef typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator OctreeIterator;
    typedef KernelClass<ParticleClass, CellClass, OctreeHeight> Kernel;

    Octree* const tree;                  //< The octree to work on
    Kernel* kernels[FThreadNumbers];          //< The kernels

    FDEBUG(FTic counterTime);                //< In case of debug: to count the elapsed time
    FDEBUG(FTic computationCounter);     //< In case of debug: to  count computation time

    OctreeIterator* iterArray;

public:	
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmArray(Octree* const inTree, Kernel* const inKernels)
                      : tree(inTree) , iterArray(0) {

        assert(tree, "tree cannot be null", __LINE__, __FILE__);
        assert(kernels, "kernels cannot be null", __LINE__, __FILE__);

        for(int idxThread = 0 ; idxThread < FThreadNumbers ; ++idxThread){
            this->kernels[idxThread] = new KernelClass<ParticleClass, CellClass, OctreeHeight>(*inKernels);
        }

        FDEBUG(FDebug::Controller << "FFmmAlgorithmArray\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmArray(){
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

        for(int idxThread = 0 ; idxThread < FThreadNumbers ; ++idxThread){
            this->kernels[idxThread]->init();
        }

        bottomPass();
        upwardPass();

        downardPass();

        directPass();

        delete [] iterArray;
        iterArray = 0;

        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** P2M */
    void bottomPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG( counterTime.tic() );

        OctreeIterator octreeIterator(tree);
        int leafs = 0;
        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            iterArray[leafs] = octreeIterator;
            ++leafs;
        } while(octreeIterator.moveRight());

        FDEBUG(computationCounter.tic());
        #pragma omp parallel num_threads(FThreadNumbers)
        {
            Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for
            for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                // We need the current cell that represent the leaf
                // and the list of particles
                myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentListSources());
            }
        }
        FDEBUG(computationCounter.tac());

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

            FDEBUG(computationCounter.tic());
            #pragma omp parallel num_threads(FThreadNumbers)
            {
                Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                #pragma omp for
                for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                    // We need the current cell and the child
                    // child is an array (of 8 child) that may be null
                    myThreadkernels->M2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentChild(), idxLevel);
                }
            }
            FDEBUG(computationCounter.tac());
            FDEBUG(totalComputation += computationCounter.elapsed());
        }

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

                FDEBUG(computationCounter.tic());
                #pragma omp parallel num_threads(FThreadNumbers)
                {
                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                    CellClass* neighbors[208];
                    #pragma omp for
                    for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                        const int counter = tree->getDistantNeighbors(neighbors,  iterArray[idxLeafs].getCurrentGlobalIndex(),idxLevel);
                        if(counter) myThreadkernels->M2L(  iterArray[idxLeafs].getCurrentCell() , neighbors, counter, idxLevel);
                    }
                }
                FDEBUG(computationCounter.tac());
                FDEBUG(totalComputation += computationCounter.elapsed());
            }
        }
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
            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                int leafs = 0;
                // for each cells
                do{
                    iterArray[leafs] = octreeIterator;
                    ++leafs;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                FDEBUG(computationCounter.tic());
                #pragma omp parallel num_threads(FThreadNumbers)
                {
                    Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
                    #pragma omp for
                    for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                        myThreadkernels->L2L( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentChild(), idxLevel);
                    }
                }
                FDEBUG(computationCounter.tac());
                FDEBUG(totalComputation += computationCounter.elapsed());
            }
        }

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

        const int heightMinusOne = OctreeHeight - 1;
        FDEBUG(computationCounter.tic());
        #pragma omp parallel num_threads(FThreadNumbers)
        {
            Kernel * const myThreadkernels = kernels[omp_get_thread_num()];
            // There is a maximum of 26 neighbors
            FList<ParticleClass*>* neighbors[26];

            #pragma omp for
            for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                myThreadkernels->L2P(iterArray[idxLeafs].getCurrentCell(), iterArray[idxLeafs].getCurrentListTargets());
                // need the current particles and neighbors particles
                const int counter = tree->getLeafsNeighbors(neighbors, iterArray[idxLeafs].getCurrentGlobalIndex(),heightMinusOne);
                myThreadkernels->P2P( iterArray[idxLeafs].getCurrentListTargets(), iterArray[idxLeafs].getCurrentListSources() , neighbors, counter);
            }
        }
        FDEBUG(computationCounter.tac());

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

};


#endif //FFmmALGORITHMARRAY_HPP

// [--LICENSE--]
