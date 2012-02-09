// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FFMMALGORITHMTSM_HPP
#define FFMMALGORITHMTSM_HPP


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
*
* The differences with FmmAlgorithm is that it used target source model.
*/
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmTsm : protected FAssertable{

    OctreeClass* const tree;                                                     //< The octree to work on
    KernelClass* const kernels;    //< The kernels

    const int OctreeHeight;

    FDEBUG(FTic counterTime);                                               //< In case of debug: to count the elapsed time
    FDEBUG(FTic computationCounter);                                        //< In case of debug: to  count computation time

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmTsm(OctreeClass* const inTree, KernelClass* const inKernels)
        : tree(inTree) , kernels(inKernels) , OctreeHeight(tree->getHeight()){

        fassert(tree, "tree cannot be null", __LINE__, __FILE__);
        fassert(kernels, "kernels cannot be null", __LINE__, __FILE__);

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
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

        bottomPass();

        upwardPass();

        transferPass();

        downardPass();

        directPass();

    }

    /** P2M */
    void bottomPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

        typename OctreeClass::Iterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            // We need the current cell that represent the leaf
            // and the list of particles
            FDEBUG(computationCounter.tic());
            ContainerClass* const sources = octreeIterator.getCurrentListSrc();
            if(sources->getSize()){
                octreeIterator.getCurrentCell()->setSrcChildTrue();
                kernels->P2M( octreeIterator.getCurrentCell() , sources);
            }
            if(octreeIterator.getCurrentListTargets()->getSize()){
                octreeIterator.getCurrentCell()->setTargetsChildTrue();
            }
            FDEBUG(computationCounter.tac());
            FDEBUG(totalComputation += computationCounter.elapsed());
        } while(octreeIterator.moveRight());

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );

    }

    /** M2M */
    void upwardPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Upward Pass\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

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

                CellClass* potentialChild[8];
                CellClass** const realChild = octreeIterator.getCurrentChild();
                CellClass* const currentCell = octreeIterator.getCurrentCell();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    potentialChild[idxChild] = 0;
                    if(realChild[idxChild]){
                        if(realChild[idxChild]->hasSrcChild()){
                            currentCell->setSrcChildTrue();
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
        FDEBUG( FDebug::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );

    }

    /** M2L */
    void transferPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.moveDown();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const CellClass* neighbors[189];

        // for each levels
        for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
            // for each cells
            do{
                FDEBUG(computationCounter.tic());
                CellClass* const currentCell = octreeIterator.getCurrentCell();
                if(currentCell->hasTargetsChild()){
                    const int counter = tree->getDistantNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),idxLevel);
                    int offsetTargetNeighbors = 0;
                    for(int idxRealNeighbors = 0 ; idxRealNeighbors < counter ; ++idxRealNeighbors, ++offsetTargetNeighbors){
                        if(neighbors[idxRealNeighbors]->hasSrcChild()){
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

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );
    }

    /** L2L */
    void downardPass(){
        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.moveDown();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

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

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << totalComputation << " s\n" );
    }



    /** P2P */
    void directPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG( counterTime.tic() );
        FDEBUG( double totalComputation = 0 );

        const int heightMinusOne = OctreeHeight - 1;

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        // There is a maximum of 26 neighbors
        ContainerClass* neighbors[26];
        // for each leafs
        do{
            FDEBUG(computationCounter.tic());
            kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
            // need the current particles and neighbors particles
            const int counter = tree->getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalIndex(),heightMinusOne);
            kernels->P2P( octreeIterator.getCurrentListTargets(), octreeIterator.getCurrentListSrc() , neighbors, counter);
            FDEBUG(computationCounter.tac());
            FDEBUG(totalComputation += computationCounter.elapsed());
        } while(octreeIterator.moveRight());

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation L2P + P2P : " << totalComputation << " s\n" );

    }

};


#endif //FFMMALGORITHMTSM_HPP


