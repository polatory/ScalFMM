#ifndef FFMMALGORITHM_HPP
#define FFMMALGORITHM_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTic.hpp"

#include "../Containers/FOctree.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFMMAlgorithm
* @brief
* Please read the license
*
* This class is a basic FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguements.
*/
template<template< class ParticuleClass, class CellClass> class KernelClass, class ParticuleClass, class CellClass, int OctreeHeight, int SubtreeHeight>
class FFMMAlgorithm : protected FAssertable{
    // To reduce the size of variable type based on foctree in this file
    typedef FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight> Octree;
    typedef typename FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>::Iterator FOctreeIterator;

    Octree* const tree;                                         //< The octree to work on
    KernelClass<ParticuleClass, CellClass>* const kernels;      //< The kernels

    FDEBUG_TIME(FTic counter);                                  //< In case of debug: to count the time

public:	
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFMMAlgorithm(Octree* const inTree,
                  KernelClass<ParticuleClass, CellClass>* const inKernels)
                      : tree(inTree) , kernels(inKernels) {

        assert(tree, "tree cannot be null", __LINE__, __FILE__);
        assert(kernels, "kernels cannot be null", __LINE__, __FILE__);

        FDEBUG_TRACE(FDebug::Controller << "FFMMAlgorithm\n");
    }

    /** Default destructor */
    virtual ~FFMMAlgorithm(){
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        kernels->init();

        bottomPass();
        upwardPass();

        downardPass();

        directPass();
    }

    /** P2M */
    void bottomPass(){
        FDEBUG_TRACE( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG_TIME( counter.tic() );

        FOctreeIterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            // We need the current cell that represent the leaf
            // and the list of particules
            kernels->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentList());
        } while(octreeIterator.moveRight());

        FDEBUG_TIME( counter.tac() );
        FDEBUG_TRACE( FDebug::Controller << "\tFinished (" FDEBUG_TIME(<< counter.elapsed() <<) "s)\n" );
    }

    /** M2M */
    void upwardPass(){
        FDEBUG_TRACE( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG_TIME( counter.tic() );

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
                kernels->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();
        }

        FDEBUG_TIME( counter.tac() );
        FDEBUG_TRACE( FDebug::Controller << "\tFinished (" FDEBUG_TIME(<< counter.elapsed() <<) "s)\n" );
    }

    /** M2L L2L */
    void downardPass(){
        FDEBUG_TRACE( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
        FDEBUG_TIME( counter.tic() );

        { // first M2L
            FOctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();

            FOctreeIterator avoidGotoLeftIterator(octreeIterator);

            CellClass* neighbors[208];
            // for each levels
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                // for each cells
                do{
                    const int counter = tree->getDistantNeighbors(neighbors, octreeIterator.getCurrentGlobalIndex(),idxLevel);
                    if(counter) kernels->M2L( octreeIterator.getCurrentCell() , neighbors, counter, idxLevel);
                } while(octreeIterator.moveRight());

                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;
            }
        }
        FDEBUG_TIME( counter.tac() );
        FDEBUG_TRACE( FDebug::Controller << "\tFinished (" FDEBUG_TIME(<< counter.elapsed() <<) "s)\n" );

        FDEBUG_TRACE( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
        FDEBUG_TIME( counter.tic() );
        { // second L2L
            FOctreeIterator octreeIterator(tree);
            octreeIterator.moveDown();

            FOctreeIterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = OctreeHeight - 1;
            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                // for each cells
                do{
                    kernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                } while(octreeIterator.moveRight());

                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;
            }
        }

        FDEBUG_TIME( counter.tac() );
        FDEBUG_TRACE( FDebug::Controller << "\tFinished (" FDEBUG_TIME(<< counter.elapsed() <<) "s)\n" );

    }

    /** P2P */
    void directPass(){
        FDEBUG_TRACE( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG_TIME( counter.tic() );

        const int heightMinusOne = OctreeHeight - 1;

        FOctreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        // There is a maximum of 26 neighbors
        FList<ParticuleClass*>* neighbors[26];
        // for each leafs
        do{
            kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentList());
            // need the current particules and neighbors particules
            const int counter = tree->getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalIndex(),heightMinusOne);
            kernels->P2P( octreeIterator.getCurrentList() , neighbors, counter);
        } while(octreeIterator.moveRight());

        FDEBUG_TIME( counter.tac() );
        FDEBUG_TRACE( FDebug::Controller << "\tFinished (" FDEBUG_TIME(<< counter.elapsed() <<) "s)\n" );

    }

};


#endif //FFMMALGORITHM_HPP

// [--LICENSE--]
