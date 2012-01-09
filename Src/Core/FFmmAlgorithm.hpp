#ifndef FFMMALGORITHM_HPP
#define FFMMALGORITHM_HPP
// [--License--]

#include "../Utils/FGlobal.hpp"
#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithm
* @brief
* Please read the license
*
* This class is a basic FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguements.
*/
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithm : protected FAssertable{

    OctreeClass* const tree;       //< The octree to work on
    KernelClass* const kernels;    //< The kernels

    const int OctreeHeight;

public:	
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithm(OctreeClass* const inTree, KernelClass* const inKernels)
                      : tree(inTree) , kernels(inKernels), OctreeHeight(tree->getHeight()) {

        fassert(tree, "tree cannot be null", __LINE__, __FILE__);
        fassert(kernels, "kernels cannot be null", __LINE__, __FILE__);

        FDEBUG(FDebug::Controller << "FFmmAlgorithm\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithm(){
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

        bottomPass();

        upwardPass();

        downardPass();

        directPass();         
    }

private:
    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M */
    void bottomPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic computationCounter);

        typename OctreeClass::Iterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            // We need the current cell that represent the leaf
            // and the list of particles
            FDEBUG(computationCounter.tic());
            kernels->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
            FDEBUG(computationCounter.tac());
        } while(octreeIterator.moveRight());

        FDEBUG( FDebug::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic computationCounter);

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
            // for each cells
            do{
                // We need the current cell and the child
                // child is an array (of 8 child) that may be null
                FDEBUG(computationCounter.tic());
                kernels->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                FDEBUG(computationCounter.tac());
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();
        }


        FDEBUG( FDebug::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downward
    /////////////////////////////////////////////////////////////////////////////

    /** M2L L2L */
    void downardPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

        { // first M2L
            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            FDEBUG(FTic computationCounter);

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.moveDown();

            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            const CellClass* neighbors[189];

            // for each levels
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                // for each cells
                do{
                    const int counter = tree->getDistantNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);
                    FDEBUG(computationCounter.tic());
                    if(counter) kernels->M2L( octreeIterator.getCurrentCell() , neighbors, counter, idxLevel);
                    FDEBUG(computationCounter.tac());
                } while(octreeIterator.moveRight());

                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;
            }
            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        }

        { // second L2L            
            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            FDEBUG(FTic computationCounter );

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.moveDown();

            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = OctreeHeight - 1;
            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                // for each cells
                do{
                    FDEBUG(computationCounter.tic());
                    kernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                    FDEBUG(computationCounter.tac());
                } while(octreeIterator.moveRight());

                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;
            }

            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        }


    }

    /////////////////////////////////////////////////////////////////////////////
    // Direct
    /////////////////////////////////////////////////////////////////////////////

    /** P2P */
    void directPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic computationCounterL2P);
        FDEBUG(FTic computationCounterP2P);

        const int heightMinusOne = OctreeHeight - 1;

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        // There is a maximum of 26 neighbors
        ContainerClass* neighbors[26];
        MortonIndex neighborsIndex[26];
        // for each leafs
        do{
            FDEBUG(computationCounterL2P.tic());
            ///kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
            FDEBUG(computationCounterL2P.tac());
            // need the current particles and neighbors particles
            const int counter = tree->getLeafsNeighborsWithIndex(neighbors, neighborsIndex, octreeIterator.getCurrentGlobalIndex(),heightMinusOne);
            FDEBUG(computationCounterP2P.tic());
            kernels->P2P(octreeIterator.getCurrentGlobalIndex(),octreeIterator.getCurrentListTargets(), octreeIterator.getCurrentListSrc() , neighbors, neighborsIndex, counter);
            FDEBUG(computationCounterP2P.tac());
        } while(octreeIterator.moveRight());


        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation L2P : " << computationCounterL2P.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation P2P : " << computationCounterP2P.cumulated() << " s\n" );

    }

};


#endif //FFMMALGORITHM_HPP

// [--END--]
