// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================
#ifndef FFMMALGORITHMPERIODIC_HPP
#define FFMMALGORITHMPERIODIC_HPP


#include "../Utils/FGlobal.hpp"
#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmPeriodic
* @brief
* Please read the license
*
* This class is a basic FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguements.
*/
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmPeriodic : protected FAssertable{

    OctreeClass* const tree;       //< The octree to work on
    KernelClass* const kernels;    //< The kernels

    const int OctreeHeight;
    const int periodicLimit;

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmPeriodic(OctreeClass* const inTree, KernelClass* const inKernels, const int inPeriodicLimit = 0)
                      : tree(inTree) , kernels(inKernels), OctreeHeight(tree->getHeight()), periodicLimit(inPeriodicLimit) {

        fassert(tree, "tree cannot be null", __LINE__, __FILE__);
        fassert(kernels, "kernels cannot be null", __LINE__, __FILE__);

        FDEBUG(FDebug::Controller << "FFmmAlgorithmPeriodic\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmPeriodic(){
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
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 0 ; --idxLevel ){
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
            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            const CellClass* neighbors[189];
            FTreeCoordinate relativePosition[189];

            // for each levels
            for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
                // for each cells
                do{
                    const int counter = tree->getDistantNeighbors(neighbors, relativePosition, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);
                    FDEBUG(computationCounter.tic());
                    if(counter) kernels->M2L( octreeIterator.getCurrentCell() , neighbors, relativePosition, counter, idxLevel);
                    FDEBUG(computationCounter.tac());
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;
            }
            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        }

        processPeriodicLevels();

        { // second L2L
            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            FDEBUG(FTic computationCounter );

            typename OctreeClass::Iterator octreeIterator(tree);
            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = OctreeHeight - 1;
            // for each levels exepted leaf level
            for(int idxLevel = 1 ; idxLevel < heightMinusOne ; ++idxLevel ){
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
        FTreeCoordinate neighborsPosition[26];
        // for each leafs
        do{
            FDEBUG(computationCounterL2P.tic());
            kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
            FDEBUG(computationCounterL2P.tac());
            // need the current particles and neighbors particles
            const int counter = tree->getLeafsNeighborsWithIndex(neighbors, neighborsPosition, octreeIterator.getCurrentGlobalIndex(),heightMinusOne);
            FDEBUG(computationCounterP2P.tic());
            kernels->P2P(octreeIterator.getCurrentGlobalIndex(),octreeIterator.getCurrentListTargets(), octreeIterator.getCurrentListSrc() , neighbors, neighborsPosition, counter);
            FDEBUG(computationCounterP2P.tac());
        } while(octreeIterator.moveRight());


        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation L2P : " << computationCounterL2P.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation P2P : " << computationCounterP2P.cumulated() << " s\n" );

    }

    /////////////////////////////////////////////////////////////////////////////
    // Periodic levels = levels <= 0
    /////////////////////////////////////////////////////////////////////////////

    /** Periodicity */
    void processPeriodicLevels(){
        if( !periodicLimit ){
            return;
        }

        CellClass upperCells[periodicLimit];

        // First M2M from level 1 to level 0
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoLeft();
            kernels->M2M( &upperCells[0], octreeIterator.getCurrentBox(), 0);
        }
        // Then M2M from level 0 to level -LIMITE
        {
            CellClass* virtualChild[8];
            for(int idxLevel = 1 ; idxLevel < periodicLimit ; ++idxLevel){
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                    virtualChild[idxChild] = &upperCells[idxLevel-1];
                }
                kernels->M2M( &upperCells[idxLevel], virtualChild, -idxLevel);
            }
        }
        // Then M2L at all level
        {
            // We say that we are in the child index 0
            // So we can compute one time the relative indexes
            FTreeCoordinate relativePosition[189];
            {
                int counterPosition = 0;
                for(int idxX = -2 ; idxX <= 3 ; ++idxX){
                    for(int idxY = -2 ; idxY <= 3 ; ++idxY){
                        for(int idxZ = -2 ; idxZ <= 3 ; ++idxZ){
                            if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                                relativePosition[counterPosition++].setPosition( idxX, idxY, idxZ);
                            }
                        }
                    }
                }
            }

            const CellClass* neighbors[189];
            const int counter = 189;

            for(int idxLevel = 0 ; idxLevel < periodicLimit ; ++idxLevel ){
                for(int idxNeigh = 0 ; idxNeigh < 189 ; ++idxNeigh){
                    neighbors[idxNeigh] = &upperCells[idxLevel];
                }
                kernels->M2L( &upperCells[idxLevel] , neighbors, relativePosition, counter, -idxLevel);
            }

        }

        // Finally L2L until level 0
        {
            CellClass* virtualChild[8];
            memset(virtualChild, 0, sizeof(CellClass*) * 8);
            for(int idxLevel = periodicLimit - 1 ; idxLevel > 0  ; --idxLevel){
                virtualChild[0] = &upperCells[idxLevel-1];
                kernels->L2L( &upperCells[idxLevel], virtualChild, -idxLevel);
            }
        }

        // L2L from 0 to level 1
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoLeft();
            kernels->L2L( &upperCells[0], octreeIterator.getCurrentBox(), 0);
        }

    }

};


#endif // FFMMALGORITHMPERIODIC_HPP
