#ifndef FFMMALGORITHMTSM_HPP
#define FFMMALGORITHMTSM_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"

#include "../Containers/FOctree.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmTsm
* @brief
* Please read the license
*
* This class is a basic FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguements.
*/
template<template< class ParticleClass, class CellClass, int OctreeHeight> class KernelClass,
class ParticleClass, class CellClass,
template<class ParticleClass> class LeafClass,
int OctreeHeight, int SubtreeHeight>
        class FFmmAlgorithmTsm : protected FAssertable{
    // To reduce the size of variable type based on foctree in this file
    typedef FOctree<ParticleClass, CellClass, LeafClass, OctreeHeight, SubtreeHeight> Octree;
    typedef typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator FOctreeIterator;

    Octree* const tree;                                                     //< The octree to work on
    KernelClass<ParticleClass, CellClass, OctreeHeight>* const kernels;    //< The kernels

    FDEBUG(FTic counterTime);                                               //< In case of debug: to count the elapsed time
    FDEBUG(FTic computationCounter);                                        //< In case of debug: to  count computation time

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmTsm(Octree* const inTree, KernelClass<ParticleClass,CellClass,OctreeHeight>* const inKernels)
        : tree(inTree) , kernels(inKernels) {

        assert(tree, "tree cannot be null", __LINE__, __FILE__);
        assert(kernels, "kernels cannot be null", __LINE__, __FILE__);

        FDEBUG(FDebug::Controller << "FFmmAlgorithmTsm\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmTsm(){
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );

        kernels->init();

        bottomPass();
        upwardPass();

        downardPass();

        directPass();
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** P2M */
    void bottomPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

        FOctreeIterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            // We need the current cell that represent the leaf
            // and the list of particles
            FDEBUG(computationCounter.tic());
            FList<ParticleClass*>* const sources = octreeIterator.getCurrentListSources();
            if(sources->getSize()){
                octreeIterator.getCurrentCell()->setSourcesChildTrue();
                kernels->P2M( octreeIterator.getCurrentCell() , sources);
            }
            if(octreeIterator.getCurrentListTargets()->getSize()){
                octreeIterator.getCurrentCell()->setTargetsChildTrue();
            }
            FDEBUG(computationCounter.tac());
            FDEBUG(totalComputation += computationCounter.elapsed());
        } while(octreeIterator.moveRight());

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

        // Start from leal level - 1
        FOctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();

        FOctreeIterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
            // for each cells
            do{
                // We need the current cell and the child
                // child is an array (of 8 child) that may be null
                FDEBUG(computationCounter.tic());

                CellClass* potentialChild[8];
                CellClass** const realChild = octreeIterator.getCurrentChild();
                CellClass* const currentCell = octreeIterator.getCurrentCell();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    potentialChild[idxChild] = 0;
                    if(realChild[idxChild]){
                        if(realChild[idxChild]->hasSourcesChild()){
                            currentCell->setSourcesChildTrue();
                            potentialChild[idxChild] = realChild[idxChild];
                        }
                        if(realChild[idxChild]->hasTargetsChild()){
                            currentCell->setTargetsChildTrue();
                        }
                    }
                }
                kernels->M2M( currentCell , potentialChild, idxLevel);

                FDEBUG(computationCounter.tac());
                FDEBUG(totalComputation += computationCounter.elapsed());
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();
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
            FOctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();

            FOctreeIterator avoidGotoLeftIterator(octreeIterator);

            CellClass* neighbors[208];
            // for each levels
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                // for each cells
                do{
                    FDEBUG(computationCounter.tic());
                    CellClass* const currentCell = octreeIterator.getCurrentCell();
                    if(currentCell->hasTargetsChild()){
                        const int counter = tree->getDistantNeighbors(neighbors, octreeIterator.getCurrentGlobalIndex(),idxLevel);
                        int offsetTargetNeighbors = 0;
                        for(int idxRealNeighbors = 0 ; idxRealNeighbors < counter ; ++idxRealNeighbors, ++offsetTargetNeighbors){
                            if(neighbors[idxRealNeighbors]->hasSourcesChild()){
                                if(idxRealNeighbors != offsetTargetNeighbors){
                                    neighbors[offsetTargetNeighbors] = neighbors[idxRealNeighbors];
                                }
                            }
                            else{
                                --offsetTargetNeighbors;
                            }
                        }
                        if(offsetTargetNeighbors){
                            kernels->M2L( currentCell , neighbors, offsetTargetNeighbors, idxLevel);
                        }
                    }
                    FDEBUG(computationCounter.tac());
                    FDEBUG(totalComputation += computationCounter.elapsed());
                } while(octreeIterator.moveRight());

                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;
            }
        }
        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );

        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( totalComputation = 0 );
        { // second L2L
            FOctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();

            FOctreeIterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = OctreeHeight - 1;
            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                // for each cells
                do{
                    FDEBUG(computationCounter.tic());
                    CellClass* potentialChild[8];
                    CellClass** const realChild = octreeIterator.getCurrentChild();
                    CellClass* const currentCell = octreeIterator.getCurrentCell();
                    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                        if(realChild[idxChild] && realChild[idxChild]->hasTargetsChild()){
                            potentialChild[idxChild] = realChild[idxChild];
                        }
                        else{
                            potentialChild[idxChild] = 0;
                        }
                    }
                    kernels->L2L( currentCell , potentialChild, idxLevel);
                    FDEBUG(computationCounter.tac());
                    FDEBUG(totalComputation += computationCounter.elapsed());
                } while(octreeIterator.moveRight());

                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;
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
        FDEBUG( double totalComputation = 0 );

        const int heightMinusOne = OctreeHeight - 1;

        FOctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        // There is a maximum of 26 neighbors
        FList<ParticleClass*>* neighbors[26];
        // for each leafs
        do{
            FDEBUG(computationCounter.tic());
            kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
            // need the current particles and neighbors particles
            const int counter = tree->getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalIndex(),heightMinusOne);
            kernels->P2P( octreeIterator.getCurrentListTargets(), octreeIterator.getCurrentListSources() , neighbors, counter);
            FDEBUG(computationCounter.tac());
            FDEBUG(totalComputation += computationCounter.elapsed());
        } while(octreeIterator.moveRight());

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished ("  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
    }

};


#endif //FFMMALGORITHMTSM_HPP

// [--LICENSE--]
