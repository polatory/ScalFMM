#ifndef FFMMALGORITHM_HPP
#define FFMMALGORITHM_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"

#include "../Containers/FOctree.hpp"

#include "FAbstractKernels.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFMMAlgorithm
* @brief
* Please read the license
*
* This class is a basic FMM algorithm
* It just iterates on a tree and call the kernls with good arguments
*/
template< class ParticuleClass, class CellClass, int OctreeHeight, int SubtreeHeight>
class FFMMAlgorithm : public FAssertable{
    FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>* const tree;    //< The octree to work on
    FAbstractKernels<ParticuleClass, CellClass>* const kernels;                     //< The kernels

public:	
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree
      * @param inKernels the kernels
      * an assert is launched if one of the arguments is null
      */
    FFMMAlgorithm(FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>* const inTree,
                  FAbstractKernels<ParticuleClass, CellClass>* const inKernels)
        : tree(inTree) , kernels(inKernels) {
        assert(tree, "tree cannot be null", __LINE__, __FILE__);
        assert(kernels, "kernels cannot be null", __LINE__, __FILE__);
    }

    /** Default destructor */
    virtual ~FFMMAlgorithm(){
    }

    /** To execute the fmm algorithm
      * Call this function to run the complete algo
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
        FDEBUG( FDebug::Controller.write("Start Bottom Pass\n") );

        typename FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>::Iterator octreeIterator(tree);
        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            // We need the current cell that represent the leaf
            // and the list of particules
            kernels->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentList());
        } while(octreeIterator.moveRight());

        FDEBUG( FDebug::Controller.write("Finished\n"); )
    }

    /** M2M */
    void upwardPass(){
        FDEBUG( FDebug::Controller.write("Start Upward Pass\n"); );

        typename FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();
        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
            // for each cells
            do{
                // We need the current cell and the child
                // child is an array (of 8 child) that may be null
                kernels->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild());
            } while(octreeIterator.moveRight());
            octreeIterator.moveUp();
            octreeIterator.gotoLeft();
        }

        FDEBUG( FDebug::Controller.write("Finished\n"); )
    }

    /** M2L L2L */
    void downardPass(){
        FDEBUG( FDebug::Controller.write("Start Downward Pass\n"); );

        { // first M2L
            typename FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>::Iterator octreeIterator(tree);
            octreeIterator.moveDown();
            CellClass* neighbors[208];
            // for each levels
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                // for each cells
                do{
                    const int counter = tree->getDistantNeighbors(neighbors, octreeIterator.getCurrentGlobalIndex(),idxLevel);
                    kernels->M2L( octreeIterator.getCurrentCell() , neighbors, counter);
                } while(octreeIterator.moveRight());
                octreeIterator.gotoLeft();
                octreeIterator.moveDown();
            }
        }
        { // second L2L
            typename FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>::Iterator octreeIterator(tree);
            octreeIterator.moveDown();
            const int heightMinusOne = OctreeHeight - 1;
            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                // for each cells
                do{
                    kernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild());
                } while(octreeIterator.moveRight());
                octreeIterator.gotoLeft();
                octreeIterator.moveDown();
            }
        }
        FDEBUG( FDebug::Controller.write("Finished\n"); )

    }

    /** P2P */
    void directPass(){
        FDEBUG( FDebug::Controller.write("Start Direct Pass\n"); );

        typename FOctree<ParticuleClass, CellClass, OctreeHeight, SubtreeHeight>::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        // There is a maximum of 26 neighbors
        FList<ParticuleClass*>* neighbors[26];
        // for each leafs
        do{
            kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentList());
            // need the current particules and neighbors particules
            const int counter = tree->getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalIndex(),OctreeHeight-1);
            kernels->P2P( octreeIterator.getCurrentList() , neighbors, counter);
        } while(octreeIterator.moveRight());

        FDEBUG( FDebug::Controller.write("Finished\n"); )

    }

};


#endif //FFMMALGORITHM_HPP

// [--LICENSE--]
