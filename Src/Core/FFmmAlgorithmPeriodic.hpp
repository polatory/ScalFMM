#ifndef FFMMALGORITHMPERIODIC_HPP
#define FFMMALGORITHMPERIODIC_HPP


#include "../Utils/FGlobal.hpp"
#include "../Utils/FGlobalPeriodic.hpp"
#include "../Utils/FAssert.hpp"
#include "../Utils/FLog.hpp"

#include "../Utils/FTic.hpp"
#include "../Utils/FMemUtils.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"

#include "FCoreCommon.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmPeriodic
* @brief
* Please read the license
*
* This class is a basic FMM algorithm with periodic behavior
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguments.
*/
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmPeriodic : public FAbstractAlgorithm{

    OctreeClass* const tree;        //< The octree to work on
    KernelClass* kernels;           //< The kernels

    const int OctreeHeight;         //< The height of the octree (real height)
    const int nbLevelsAboveRoot;    //< The nb of level the user ask to go above the tree (>= -1)
    const int offsetRealTree;       //< nbLevelsAboveRoot GetFackLevel


public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      * @param inUpperLevel this parameter defines the behavior of the periodicity refer to the main doc
      *
      */
    FFmmAlgorithmPeriodic(OctreeClass* const inTree, const int inUpperLevel = 0)
        : tree(inTree) , kernels(nullptr), OctreeHeight(tree->getHeight()),
          nbLevelsAboveRoot(inUpperLevel), offsetRealTree(inUpperLevel + 3) {

        FAssertLF(tree, "tree cannot be null");
        FAssertLF(-1 <= inUpperLevel, "inUpperLevel cannot be < -1");

        FLOG(FLog::Controller << "FFmmAlgorithmPeriodic\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmPeriodic(){
    }

    void setKernel(KernelClass*const inKernel){
        kernels = inKernel;
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(const unsigned operationsToProceed = FFmmNearAndFarFields){
        FAssertLF(kernels, "kernels cannot be null");

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L){
            transferPass();
            // before downward pass we have to perform the periodicity
            processPeriodicLevels();
        }

        if(operationsToProceed & FFmmL2L) downardPass();

        if((operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P)) directPass();
    }


    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M */
    void bottomPass(){
        FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        typename OctreeClass::Iterator octreeIterator(tree);

        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            // We need the current cell that represent the leaf
            // and the list of particles
            FLOG(computationCounter.tic());
            kernels->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
            FLOG(computationCounter.tac());
        } while(octreeIterator.moveRight());

        FLOG( FLog::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
        FLOG( FLog::Controller.write("\tStart Upward Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 0 ; --idxLevel ){
            FLOG(FTic counterTimeLevel);
            const int fackLevel = idxLevel + offsetRealTree;
            // for each cells
            do{
                // We need the current cell and the child
                // child is an array (of 8 child) that may be null
                FLOG(computationCounter.tic());
                kernels->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), fackLevel);
                FLOG(computationCounter.tac());
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();
            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << "(" << fackLevel << ") = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
        }


        FLOG( FLog::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Transfer
    /////////////////////////////////////////////////////////////////////////////

    /** M2L L2L */
    void transferPass(){
        FLOG( FLog::Controller.write("\tStart Downward Pass (M2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        typename OctreeClass::Iterator octreeIterator(tree);
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const CellClass* neighbors[343];

        // for each levels
        for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
            FLOG(FTic counterTimeLevel);
            const int fackLevel = idxLevel + offsetRealTree;
            // for each cells
            do{
                const int counter = tree->getPeriodicInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel, AllDirs);
                FLOG(computationCounter.tic());
                if(counter) kernels->M2L( octreeIterator.getCurrentCell() , neighbors, counter, fackLevel);
                FLOG(computationCounter.tac());
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;

            FLOG(computationCounter.tic());
            kernels->finishedLevelM2L(fackLevel);
            FLOG(computationCounter.tac());
            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << "(" << fackLevel << ") = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
        }
        FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downward
    /////////////////////////////////////////////////////////////////////////////


    void downardPass(){ // second L2L
        FLOG( FLog::Controller.write("\tStart Downward Pass (L2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter );

        typename OctreeClass::Iterator octreeIterator(tree);
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = OctreeHeight - 1;
        // for each levels exepted leaf level
        for(int idxLevel = 1 ; idxLevel < heightMinusOne ; ++idxLevel ){
            FLOG(FTic counterTimeLevel);
            const int fackLevel = idxLevel + offsetRealTree;

            // for each cells
            do{
                FLOG(computationCounter.tic());
                kernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), fackLevel);
                FLOG(computationCounter.tac());
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;
            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << "(" << fackLevel << ") = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
        }

        FLOG( FLog::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );


    }

    /////////////////////////////////////////////////////////////////////////////
    // Direct
    /////////////////////////////////////////////////////////////////////////////

    /** P2P */
    void directPass(){
        FLOG( FLog::Controller.write("\tStart Direct Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounterL2P);
        FLOG(FTic computationCounterP2P);

        const int heightMinusOne = OctreeHeight - 1;
        const FReal boxWidth = tree->getBoxWidth();

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        // There is a maximum of 26 neighbors
        ContainerClass* neighbors[27];
        FTreeCoordinate offsets[27];
        bool hasPeriodicLeaves;
        // for each leafs
        do{
            FLOG(computationCounterL2P.tic());
            kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
            FLOG(computationCounterL2P.tac());

            // need the current particles and neighbors particles
            const FTreeCoordinate centerOfLeaf = octreeIterator.getCurrentGlobalCoordinate();
            const int counter = tree->getPeriodicLeafsNeighbors( neighbors, offsets, &hasPeriodicLeaves, centerOfLeaf, heightMinusOne, AllDirs);
            int periodicNeighborsCounter = 0;

            if(hasPeriodicLeaves){
                ContainerClass* periodicNeighbors[27];
                memset(periodicNeighbors, 0, 27 * sizeof(ContainerClass*));

                for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
                    if( neighbors[idxNeig] && !offsets[idxNeig].equals(0,0,0) ){
                        // Put periodic neighbors into other array
                        periodicNeighbors[idxNeig] = neighbors[idxNeig];
                        neighbors[idxNeig] = nullptr;
                        ++periodicNeighborsCounter;

                        FReal*const positionsX = periodicNeighbors[idxNeig]->getWPositions()[0];
                        FReal*const positionsY = periodicNeighbors[idxNeig]->getWPositions()[1];
                        FReal*const positionsZ = periodicNeighbors[idxNeig]->getWPositions()[2];

                        for(int idxPart = 0; idxPart < periodicNeighbors[idxNeig]->getNbParticles() ; ++idxPart){
                            positionsX[idxPart] += boxWidth * FReal(offsets[idxNeig].getX());
                            positionsY[idxPart] += boxWidth * FReal(offsets[idxNeig].getY());
                            positionsZ[idxPart] += boxWidth * FReal(offsets[idxNeig].getZ());
                        }
                    }
                }

                FLOG(computationCounterP2P.tic());
                kernels->P2PRemote(octreeIterator.getCurrentGlobalCoordinate(),octreeIterator.getCurrentListTargets(),
                             octreeIterator.getCurrentListSrc(), periodicNeighbors, periodicNeighborsCounter);
                FLOG(computationCounterP2P.tac());

                for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
                    if( periodicNeighbors[idxNeig] ){
                        FReal*const positionsX = periodicNeighbors[idxNeig]->getWPositions()[0];
                        FReal*const positionsY = periodicNeighbors[idxNeig]->getWPositions()[1];
                        FReal*const positionsZ = periodicNeighbors[idxNeig]->getWPositions()[2];

                        for(int idxPart = 0; idxPart < periodicNeighbors[idxNeig]->getNbParticles() ; ++idxPart){
                            positionsX[idxPart] -= boxWidth * FReal(offsets[idxNeig].getX());
                            positionsY[idxPart] -= boxWidth * FReal(offsets[idxNeig].getY());
                            positionsZ[idxPart] -= boxWidth * FReal(offsets[idxNeig].getZ());
                        }
                    }
                }
            }

            FLOG(computationCounterP2P.tic());
            kernels->P2P(octreeIterator.getCurrentGlobalCoordinate(),octreeIterator.getCurrentListTargets(),
                         octreeIterator.getCurrentListSrc(), neighbors, counter - periodicNeighborsCounter);
            FLOG(computationCounterP2P.tac());


        } while(octreeIterator.moveRight());


        FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation L2P : " << computationCounterL2P.cumulated() << " s\n" );
        FLOG( FLog::Controller << "\t\t Computation P2P : " << computationCounterP2P.cumulated() << " s\n" );

    }

    /////////////////////////////////////////////////////////////////////////////
    // Periodic levels = levels <= 0
    /////////////////////////////////////////////////////////////////////////////


    /** This function process several M2M from level nbLevelsAboveRoot to level 0
      * and give the final result
      * @param result the cell at the last M2M
      * @param root the starting cell
      * @param startX the beginning of the index in x [0;endX]
      * @param endX the end of the index in x [startX;1]
      * @param startY the beginning of the index in y [0;endY]
      * @param endY the end of the index in y [startY;1]
      * @param startZ the beginning of the index in z [0;endZ]
      * @param endZ the end of the index in z [startZ;1]
      */
    void processTopM2MInIntervals( CellClass*const result, const CellClass& root, const int startX,
                              const int endX, const int startY, const int endY, const int startZ,
                              const int endZ){
        // allocate array
        CellClass*const cellsAtLevel = new CellClass[nbLevelsAboveRoot+2];
        // process by using other function
        processM2MInIntervals(cellsAtLevel,root,startX,endX,startY,endY,startZ,endZ);
        // copy result
        *result = cellsAtLevel[0];
        delete[] cellsAtLevel;
    }

    /** This function process several M2M from level nbLevelsAboveRoot to level 0
      * @param cellsAtLevel the intermediate results
      * @param root the starting cell
      * @param startX the beginning of the index in x [0;endX]
      * @param endX the end of the index in x [startX;1]
      * @param startY the beginning of the index in y [0;endY]
      * @param endY the end of the index in y [startY;1]
      * @param startZ the beginning of the index in z [0;endZ]
      * @param endZ the end of the index in z [startZ;1]
      */
    void  processM2MInIntervals( CellClass cellsAtLevel[], const CellClass& root, const int startX,
                              const int endX, const int startY, const int endY, const int startZ,
                              const int endZ){
        // start from the initial cell
        cellsAtLevel[nbLevelsAboveRoot+1] = root;
        // to create virtual children
        CellClass* virtualChild[8];
        // for all levels
        for(int idxLevel = nbLevelsAboveRoot ; idxLevel >= 0  ; --idxLevel){
            // reset children
            memset(virtualChild, 0, sizeof(CellClass*)*8);
            // fill the vector with previous result
            for(int idxX = startX ; idxX <= endX ; ++idxX){
                for(int idxY = startY ; idxY <= endY ; ++idxY){
                    for(int idxZ = startZ ; idxZ <= endZ ; ++idxZ){
                        virtualChild[childIndex(idxX,idxY,idxZ)] = &cellsAtLevel[idxLevel+1];
                    }
                }
            }
            // compute the M2M
            kernels->M2M( &cellsAtLevel[idxLevel], virtualChild, idxLevel + 2);
        }
    }

    /** Fill an interactions neighbors with some intervals
      * @param neighbors the vector to fill
      * @param source the source cell to fill the vector
      * @param startX the beginning of the index in x [-3;0]
      * @param endX the end of the index in x  [0;3]
      * @param startY the beginning of the index in y [-3;0]
      * @param endY the end of the index in y [0;3]
      * @param startZ the beginning of the index in z [-3;0]
      * @param endZ the end of the index in z [0;3]
      * @return the number of position filled
      */
    int  fillM2LVectorFromIntervals(const CellClass* neighbors[343], const CellClass& source,
                     const int startX, const int endX, const int startY, const int endY,
                     const int startZ, const int endZ){
        int counter = 0;
        // for all x in interval
        for(int idxX = startX ; idxX <= endX ; ++idxX){
            // for all y in interval
            for(int idxY = startY ; idxY <= endY ; ++idxY){
                // for all z in interval
                for(int idxZ = startZ ; idxZ <= endZ ; ++idxZ){
                    // do not fill close neigbors
                    if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1 ){
                        neighbors[neighIndex(idxX,idxY,idxZ)] = &source;
                        ++counter;
                    }
                }
            }
        }
        // return the number of position filled
        return counter;
    }

    /** Get the index of a child (for the M2M and the L2L)
      * @param x the x position in the children  (from 0 to +1)
      * @param y the y position in the children  (from 0 to +1)
      * @param z the z position in the children  (from 0 to +1)
      * @return the index (from 0 to 7)
      */
    int childIndex(const int x, const int y, const int z) const {
        return (x<<2) | (y<<1) | z;
    }

    /** Get the index of a interaction neighbors (for M2L)
      * @param x the x position in the interactions (from -3 to +3)
      * @param y the y position in the interactions (from -3 to +3)
      * @param z the z position in the interactions (from -3 to +3)
      * @return the index (from 0 to 342)
      */
    int neighIndex(const int x, const int y, const int z) const {
        return (((x+3)*7) + (y+3))*7 + (z + 3);
    }


    long long int theoricalRepetition() const {
        if( nbLevelsAboveRoot == -1 ){
            // we know it is 3 (-1;+1)
            return 3;
        }
        // Else we find the repetition in one dir and double it
        const long long int oneDirectionRepetition = (1<<(nbLevelsAboveRoot+2)); // 2^nbLevelsAboveRoot in each dim
        const long long int fullRepetition = 2 * oneDirectionRepetition;
        return fullRepetition;
    }


    void repetitionsIntervals(FTreeCoordinate*const min, FTreeCoordinate*const max) const {
        if( nbLevelsAboveRoot == -1 ){
            // We know it is (-1;1)
            min->setPosition(-1,-1,-1);
            max->setPosition(1,1,1);
        }
        else{
            const int halfRepeated = int(theoricalRepetition()/2);
            min->setPosition(-halfRepeated,-halfRepeated,-halfRepeated);
            // if we repeat the box 8 times, we go from [-4 to 3]
            max->setPosition(halfRepeated-1,halfRepeated-1,halfRepeated-1);
        }
    }


    FReal extendedBoxWidth() const {
        // The simulation box is repeated is repeated 4 times if nbLevelsAboveRoot==-1
        // And then it doubles by two
        return tree->getBoxWidth() * FReal(1<<(nbLevelsAboveRoot+3));
    }

    /** This function has to be used to init the kernel with correct args
      * it return the box cneter seen from a kernel point of view from the periodicity the user ask for
      * this is computed using the originalBoxWidth and originalBoxCenter given in parameter
      * @param originalBoxCenter the real system center
      * @param originalBoxWidth the real system size
      * @return the center the kernel should use
      */
    FPoint extendedBoxCenter() const {
        const FReal originalBoxWidth     = tree->getBoxWidth();
        const FReal originalBoxWidthDiv2 = originalBoxWidth/2.0;
        const FPoint originalBoxCenter   = tree->getBoxCenter();

        const FReal offset = extendedBoxWidth()/FReal(2.0);
        return FPoint( originalBoxCenter.getX() - originalBoxWidthDiv2 + offset,
                       originalBoxCenter.getY() - originalBoxWidthDiv2 + offset,
                       originalBoxCenter.getZ() - originalBoxWidthDiv2 + offset);
    }

    /** This function has to be used to init the kernel with correct args
      * it return the tree heigh seen from a kernel point of view from the periodicity the user ask for
      * this is computed using the originalTreeHeight given in parameter
      * @param originalTreeHeight the real tree heigh
      * @return the heigh the kernel should use
      */
    int extendedTreeHeight() const {
        // The real height
        return OctreeHeight + offsetRealTree;
    }

    /** Periodicity Core
      * This function is split in several part:
      * 1 - special case managment
      * There is nothing to do if nbLevelsAboveRoot == -1 and only
      * a M2L if nbLevelsAboveRoot == 0
      * 2 - if nbLevelsAboveRoot > 0
      * First we compute M2M and special M2M if needed for the border
      * Then the M2L by taking into account the periodicity directions
      * Then the border by using the precomputed M2M
      * Finally the L2L
      */
    void processPeriodicLevels(){
        FLOG( FLog::Controller.write("\tStart Periodic Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);

        if( nbLevelsAboveRoot != -1 ){
            // we will use offsetRealTree-1 cells but for simplicity allocate offsetRealTree
            // upperCells[offsetRealTree-1] is root cell
            CellClass*const upperCells = new CellClass[offsetRealTree];
            {
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoLeft();
                kernels->M2M( &upperCells[offsetRealTree-1], octreeIterator.getCurrentBox(), offsetRealTree);
            }
            {
                CellClass* virtualChild[8];
                for(int idxLevel = offsetRealTree-1 ; idxLevel > 1  ; --idxLevel){
                    FMemUtils::setall(virtualChild,&upperCells[idxLevel],8);
                    kernels->M2M( &upperCells[idxLevel-1], virtualChild, idxLevel);
                }
            }
            CellClass*const downerCells = new CellClass[offsetRealTree];

            {
                const int idxUpperLevel = 2;

                const CellClass* neighbors[343];
                memset(neighbors, 0, sizeof(CellClass*) * 343);
                int counter = 0;
                for(int idxX = -2 ; idxX <= 1 ; ++idxX){
                    for(int idxY = -2 ; idxY <= 1 ; ++idxY){
                        for(int idxZ = -2 ; idxZ <= 1 ; ++idxZ){
                            if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                                neighbors[neighIndex(idxX,idxY,idxZ)] = &upperCells[idxUpperLevel-1];
                                ++counter;
                            }
                        }
                    }
                }
                // compute M2L
                kernels->M2L( &downerCells[idxUpperLevel-1] , neighbors, counter, idxUpperLevel);
            }

            for(int idxUpperLevel = 3 ; idxUpperLevel <= offsetRealTree ; ++idxUpperLevel){
                const CellClass* neighbors[343];
                memset(neighbors, 0, sizeof(CellClass*) * 343);
                int counter = 0;
                for(int idxX = -2 ; idxX <= 3 ; ++idxX){
                    for(int idxY = -2 ; idxY <= 3 ; ++idxY){
                        for(int idxZ = -2 ; idxZ <= 3 ; ++idxZ){
                            if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                                neighbors[neighIndex(idxX,idxY,idxZ)] = &upperCells[idxUpperLevel-1];
                                ++counter;
                            }
                        }
                    }
                }

                // compute M2L
                kernels->M2L( &downerCells[idxUpperLevel-1] , neighbors, counter, idxUpperLevel);
            }

            {
                CellClass* virtualChild[8];
                memset(virtualChild, 0, sizeof(CellClass*) * 8);
                for(int idxLevel = 2 ; idxLevel <= offsetRealTree-1  ; ++idxLevel){
                    virtualChild[0] = &downerCells[idxLevel];
                    kernels->L2L( &downerCells[idxLevel-1], virtualChild, idxLevel);
                }
            }

            // L2L from 0 to level 1
            {
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoLeft();
                kernels->L2L( &downerCells[offsetRealTree-1], octreeIterator.getCurrentBox(), offsetRealTree);
            }

            delete[] upperCells;
            delete[] downerCells;
        }

        FLOG( FLog::Controller << "\tFinished (@Periodic = "  << counterTime.tacAndElapsed() << "s)\n" );
    }


};


#endif // FFMMALGORITHMPERIODIC_HPP
