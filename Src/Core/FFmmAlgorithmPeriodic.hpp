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
#ifndef FFMMALGORITHMPERIODIC_HPP
#define FFMMALGORITHMPERIODIC_HPP


#include "../Utils/FGlobal.hpp"
#include "../Utils/FGlobalPeriodic.hpp"
#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FMemUtils.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmPeriodic
* @brief
* Please read the license
*
* This class is a basic FMM algorithm with periodic behavior
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguements.
*/
template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmPeriodic : protected FAssertable{

    OctreeClass* const tree;        //< The octree to work on
    KernelClass* kernels;           //< The kernels

    const int OctreeHeight;         //< The heigh of the octree (real heigh)
    const int nbLevelsAboveRoot;    //< The nb of level the user ask to go above the tree (>= -1)
    const int offsetRealTree;       //< nbLevelsAboveRoot GetFackLevel
    const int periodicDirections;

    static int GetFackLevel(const int inLevelAboveRequiered){
        if( inLevelAboveRequiered == -1 ) return 1;
        if( inLevelAboveRequiered == 0  ) return 2;
        return inLevelAboveRequiered + 3;
    }

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      * @param inUpperLevel this parameter defins the behavior of the periodicity refer to the main doc
      *
      */
    FFmmAlgorithmPeriodic(OctreeClass* const inTree, const int inUpperLevel = 0, const int inPeriodicDirections = AllDirs)
        : tree(inTree) , kernels(0), OctreeHeight(tree->getHeight()),
          nbLevelsAboveRoot(inUpperLevel), offsetRealTree(GetFackLevel(inUpperLevel)),
          periodicDirections(inPeriodicDirections) {

        fassert(tree, "tree cannot be null", __LINE__, __FILE__);
        fassert(-1 <= inUpperLevel, "inUpperLevel cannot be < -1", __LINE__, __FILE__);

        FDEBUG(FDebug::Controller << "FFmmAlgorithmPeriodic\n");
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
    void execute(){
        fassert(kernels, "kernels cannot be null", __LINE__, __FILE__);
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

        bottomPass();

        upwardPass();

        transferPass();
        // before downward pass we have to perform the periodicity
        processPeriodicLevels();

        downardPass();

        directPass();
    }


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
            FDEBUG(FTic counterTimeLevel);
            const int fackLevel = idxLevel + offsetRealTree;
            // for each cells
            do{
                // We need the current cell and the child
                // child is an array (of 8 child) that may be null
                FDEBUG(computationCounter.tic());
                kernels->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), fackLevel);
                FDEBUG(computationCounter.tac());
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();
            FDEBUG( FDebug::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
        }


        FDEBUG( FDebug::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Transfer
    /////////////////////////////////////////////////////////////////////////////

    /** M2L L2L */
    void transferPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic computationCounter);

        typename OctreeClass::Iterator octreeIterator(tree);
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const CellClass* neighbors[343];

        // for each levels
        for(int idxLevel = 1 ; idxLevel < OctreeHeight ; ++idxLevel ){
            FDEBUG(FTic counterTimeLevel);
            const int fackLevel = idxLevel + offsetRealTree;
            // for each cells
            do{
                const int counter = tree->getPeriodicInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel, periodicDirections);
                FDEBUG(computationCounter.tic());
                if(counter) kernels->M2L( octreeIterator.getCurrentCell() , neighbors, counter, fackLevel);
                FDEBUG(computationCounter.tac());
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;

            FDEBUG(computationCounter.tic());
            kernels->finishedLevelM2L(fackLevel);
            FDEBUG(computationCounter.tac());
            FDEBUG( FDebug::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
        }
        FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downward
    /////////////////////////////////////////////////////////////////////////////


    void downardPass(){ // second L2L
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);
        FDEBUG(FTic computationCounter );

        typename OctreeClass::Iterator octreeIterator(tree);
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = OctreeHeight - 1;
        // for each levels exepted leaf level
        for(int idxLevel = 1 ; idxLevel < heightMinusOne ; ++idxLevel ){
            FDEBUG(FTic counterTimeLevel);
            const int fackLevel = idxLevel + offsetRealTree;

            // for each cells
            do{
                FDEBUG(computationCounter.tic());
                kernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), fackLevel);
                FDEBUG(computationCounter.tac());
            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;
            FDEBUG( FDebug::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
        }

        FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );


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
        const FReal boxWidth = tree->getBoxWidth();

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        // There is a maximum of 26 neighbors
        ContainerClass* neighbors[27];
        FTreeCoordinate offsets[27];
        bool hasPeriodicLeaves;
        // for each leafs
        do{
            FDEBUG(computationCounterL2P.tic());
            kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
            FDEBUG(computationCounterL2P.tac());

            // need the current particles and neighbors particles
            const FTreeCoordinate centerOfLeaf = octreeIterator.getCurrentGlobalCoordinate();
            const int counter = tree->getPeriodicLeafsNeighbors( neighbors, offsets, &hasPeriodicLeaves, centerOfLeaf, heightMinusOne, periodicDirections);
            int periodicNeighborsCounter = 0;

            if(hasPeriodicLeaves){
                ContainerClass* periodicNeighbors[27];
                memset(periodicNeighbors, 0, 27 * sizeof(ContainerClass*));

                for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
                    if( neighbors[idxNeig] && !offsets[idxNeig].equals(0,0,0) ){
                        // Put periodic neighbors into other array
                        periodicNeighbors[idxNeig] = neighbors[idxNeig];
                        neighbors[idxNeig] = 0;
                        ++periodicNeighborsCounter;
                        typename ContainerClass::BasicIterator iter(*periodicNeighbors[idxNeig]);
                        while( iter.hasNotFinished() ){
                            iter.data().incPosition(boxWidth * FReal(offsets[idxNeig].getX()),
                                                    boxWidth * FReal(offsets[idxNeig].getY()),
                                                    boxWidth * FReal(offsets[idxNeig].getZ()));
                            iter.gotoNext();
                        }
                    }
                }

                FDEBUG(computationCounterP2P.tic());
                kernels->P2PRemote(octreeIterator.getCurrentGlobalCoordinate(),octreeIterator.getCurrentListTargets(),
                             octreeIterator.getCurrentListSrc(), periodicNeighbors, periodicNeighborsCounter);
                FDEBUG(computationCounterP2P.tac());

                for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
                    if( periodicNeighbors[idxNeig] ){
                        typename ContainerClass::BasicIterator iter(*periodicNeighbors[idxNeig]);
                        while( iter.hasNotFinished() ){
                            iter.data().incPosition(-boxWidth * FReal(offsets[idxNeig].getX()),
                                                    -boxWidth * FReal(offsets[idxNeig].getY()),
                                                    -boxWidth * FReal(offsets[idxNeig].getZ()));
                            iter.gotoNext();
                        }
                    }
                }
            }

            FDEBUG(computationCounterP2P.tic());
            kernels->P2P(octreeIterator.getCurrentGlobalCoordinate(),octreeIterator.getCurrentListTargets(),
                         octreeIterator.getCurrentListSrc(), neighbors, counter - periodicNeighborsCounter);
            FDEBUG(computationCounterP2P.tac());


        } while(octreeIterator.moveRight());


        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation L2P : " << computationCounterL2P.cumulated() << " s\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation P2P : " << computationCounterP2P.cumulated() << " s\n" );

    }

    /////////////////////////////////////////////////////////////////////////////
    // Periodic levels = levels <= 0
    /////////////////////////////////////////////////////////////////////////////


    /** To know how many times the box is repeated in each direction
      * -x +x -y +y -z +z
      */
    int repeatedBox() const {
        return nbLevelsAboveRoot == -1 ? 3 : 3 * (1<<(nbLevelsAboveRoot+1)) + 1;
    }

    void repetitions(FTreeCoordinate*const min, FTreeCoordinate*const max) const {
        const int halfRepeated = (repeatedBox()-1) /2;
        min->setPosition(-ifDir(DirMinusX,halfRepeated,0),-ifDir(DirMinusY,halfRepeated,0),
                         -ifDir(DirMinusZ,halfRepeated,0));
        max->setPosition(ifDir(DirPlusX,halfRepeated,0),ifDir(DirPlusY,halfRepeated,0),
                         ifDir(DirPlusZ,halfRepeated,0));
    }

    FTreeCoordinate repetitions() const {
        const int halfRepeated = (repeatedBox()-1) /2;
        return FTreeCoordinate(ifDir(DirMinusX,halfRepeated,0) + ifDir(DirPlusX,halfRepeated,0) + 1,
                               ifDir(DirMinusY,halfRepeated,0) + ifDir(DirPlusY,halfRepeated,0) + 1,
                               ifDir(DirMinusZ,halfRepeated,0) + ifDir(DirPlusZ,halfRepeated,0) + 1);
    }

    /** This function has to be used to init the kernel with correct args
      * it return the box seen from a kernel point of view from the periodicity the user ask for
      * this is computed using the originalBoxWidth given in parameter
      * @param originalBoxWidth the real system size
      * @return the size the kernel should use
      */
    FReal boxwidthForKernel(const FReal originalBoxWidth) const {
        return originalBoxWidth * FReal(1<<offsetRealTree);
    }

    FPoint boxcenterForKernel(const FPoint originalBoxCenter, const FReal originalBoxWidth) const {
        const FReal offset = originalBoxWidth * FReal(1<<(offsetRealTree-1)) - originalBoxWidth/FReal(2.0);
        return FPoint( originalBoxCenter.getX() + offset,
                       originalBoxCenter.getY() + offset,
                       originalBoxCenter.getZ() + offset);
    }


    /** This function has to be used to init the kernel with correct args
      * it return the tree heigh seen from a kernel point of view from the periodicity the user ask for
      * this is computed using the originalTreeHeight given in parameter
      * @param originalTreeHeight the real tree heigh
      * @return the heigh the kernel should use
      */
    int treeHeightForKernel() const {
        // The real height
        return OctreeHeight + offsetRealTree;
    }

    bool usePerDir(const int testDir) const{
        return testPeriodicCondition(periodicDirections , PeriodicCondition(testDir));
    }

    template <class T>
    int ifDir(const PeriodicCondition testDir, const T& correctValue, const T& wrongValue) const {
        return (periodicDirections & testDir ? correctValue : wrongValue);
    }

    /** Periodicity */
    void processPeriodicLevels(){
        /////////////////////////////////////////////////////
        // If nb level == -1 nothing to do
        if( nbLevelsAboveRoot == -1 ){
            return;
        }
        /////////////////////////////////////////////////////
        // if nb level == 0 only M2L at real root level
        if( nbLevelsAboveRoot == 0 ){
            CellClass rootUp;
            // compute the root
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoLeft();
            kernels->M2M( &rootUp, octreeIterator.getCurrentBox(), 2);

            // build fack M2L vector from -3/+3 x/y/z
            const CellClass* neighbors[343];
            memset(neighbors, 0, sizeof(CellClass*) * 343);
            int counter = 0;
            for(int idxX = ifDir(DirMinusX,-3,0) ; idxX <= ifDir(DirPlusX,3,0) ; ++idxX){
                for(int idxY = ifDir(DirMinusY,-3,0) ; idxY <= ifDir(DirPlusY,3,0) ; ++idxY){
                    for(int idxZ = ifDir(DirMinusZ,-3,0) ; idxZ <= ifDir(DirPlusZ,3,0) ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            neighbors[neighIndex(idxX,idxY,idxZ)] = &rootUp;
                            ++counter;
                        }
                    }
                }
            }
            // compute M2L
            CellClass rootDown;
            kernels->M2L( &rootDown , neighbors, counter, 2);

            // put result in level 1
            kernels->L2L( &rootDown, octreeIterator.getCurrentBox(), 2);

            return;
        }
        /////////////////////////////////////////////////////
        // in other situation, we have to compute M2L from 0 to nbLevelsAboveRoot
        // but also at nbLevelsAboveRoot +1 for the rest
        CellClass*const upperCells = new CellClass[nbLevelsAboveRoot+2];

        CellClass*const cellsXAxis = new CellClass[nbLevelsAboveRoot+2];
        CellClass*const cellsYAxis = new CellClass[nbLevelsAboveRoot+2];
        CellClass*const cellsZAxis = new CellClass[nbLevelsAboveRoot+2];
        CellClass*const cellsXYAxis = new CellClass[nbLevelsAboveRoot+2];
        CellClass*const cellsYZAxis = new CellClass[nbLevelsAboveRoot+2];
        CellClass*const cellsXZAxis = new CellClass[nbLevelsAboveRoot+2];


        // First M2M from level 1 to level 0
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoLeft();
            kernels->M2M( &upperCells[nbLevelsAboveRoot+1], octreeIterator.getCurrentBox(), offsetRealTree);
        }

        // Then M2M from level 0 to level -LIMITE
        {
            CellClass* virtualChild[8];
            for(int idxLevel = nbLevelsAboveRoot ; idxLevel > 0  ; --idxLevel){
                FMemUtils::setall(virtualChild,&upperCells[idxLevel+1],8);
                kernels->M2M( &upperCells[idxLevel], virtualChild, idxLevel + 2);
            }

            // Cells on the axis of the center should be computed separatly.

            if( usePerDir(DirMinusX) && usePerDir(DirMinusY) ){
                FMemUtils::copyall(cellsZAxis,upperCells,nbLevelsAboveRoot+2);
            }
            else{
                getUpperSpecialCells(cellsZAxis,upperCells[nbLevelsAboveRoot+1],
                                    ifDir(DirMinusX,0,1),1,ifDir(DirMinusY,0,1),1,0,1);
            }
            if( usePerDir(DirMinusX) && usePerDir(DirMinusZ) ){
                FMemUtils::copyall(cellsYAxis,upperCells,nbLevelsAboveRoot+2);
            }
            else{
                getUpperSpecialCells(cellsYAxis,upperCells[nbLevelsAboveRoot+1],
                                    ifDir(DirMinusX,0,1),1,0,1,ifDir(DirMinusZ,0,1),1);
            }
            if( usePerDir(DirMinusY) && usePerDir(DirMinusZ) ){
                FMemUtils::copyall(cellsXAxis,upperCells,nbLevelsAboveRoot+2);
            }
            else{
                getUpperSpecialCells(cellsXAxis,upperCells[nbLevelsAboveRoot+1],
                                    0,1,ifDir(DirMinusY,0,1),1,ifDir(DirMinusZ,0,1),1);
            }

            // Then cells on the spaces should be computed separatly

            if( !usePerDir(DirMinusX) ){
                getUpperSpecialCells(cellsYZAxis,upperCells[nbLevelsAboveRoot+1],1,1,0,1,0,1);
            }
            else {
                FMemUtils::copyall(cellsYZAxis,upperCells,nbLevelsAboveRoot+2);
            }
            if( !usePerDir(DirMinusY) ){
                getUpperSpecialCells(cellsXZAxis,upperCells[nbLevelsAboveRoot+1],0,1,1,1,0,1);
            }
            else {
                FMemUtils::copyall(cellsXZAxis,upperCells,nbLevelsAboveRoot+2);
            }
            if( !usePerDir(DirMinusZ) ){
                getUpperSpecialCells(cellsXYAxis,upperCells[nbLevelsAboveRoot+1],0,1,0,1,1,1);
            }
            else {
                FMemUtils::copyall(cellsXYAxis,upperCells,nbLevelsAboveRoot+2);
            }

        }

        // Then M2L at all level
        {
            CellClass* positionedCells[343];
            memset(positionedCells, 0, 343 * sizeof(CellClass**));

            for(int idxX = ifDir(DirMinusX,-3,0) ; idxX <= ifDir(DirPlusX,3,0) ; ++idxX){
                for(int idxY = ifDir(DirMinusY,-3,0) ; idxY <= ifDir(DirPlusY,3,0) ; ++idxY){
                    for(int idxZ = ifDir(DirMinusZ,-3,0) ; idxZ <= ifDir(DirPlusZ,3,0) ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            if(idxX == 0 && idxY == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsZAxis;
                            }
                            else if(idxX == 0 && idxZ == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsYAxis;
                            }
                            else if(idxY == 0 && idxZ == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsXAxis;
                            }
                            else if(idxX == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsYZAxis;
                            }
                            else if(idxY == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsXZAxis;
                            }
                            else if(idxZ == 0){
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = cellsXYAxis;
                            }
                            else{
                                positionedCells[neighIndex(idxX,idxY,idxZ)] = upperCells;
                            }
                        }
                    }
                }
            }

            // We say that we are in the child index 0
            // So we can compute one time the relative indexes
            const CellClass* neighbors[343];
            memset(neighbors, 0, sizeof(CellClass*) * 343);
            int counter = 0;
            for(int idxX = ifDir(DirMinusX,-3,0) ; idxX <= ifDir(DirPlusX,2,0) ; ++idxX){
                for(int idxY = ifDir(DirMinusY,-3,0) ; idxY <= ifDir(DirPlusY,2,0) ; ++idxY){
                    for(int idxZ = ifDir(DirMinusZ,-3,0) ; idxZ <= ifDir(DirPlusZ,2,0) ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            neighbors[neighIndex(idxX,idxY,idxZ)] = reinterpret_cast<const CellClass*>(~0);
                            ++counter;
                        }
                    }
                }
            }

            for(int idxLevel = nbLevelsAboveRoot + 1 ; idxLevel > 1 ; --idxLevel ){
                for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
                    if(neighbors[idxNeigh]){
                        neighbors[idxNeigh] = &positionedCells[idxNeigh][idxLevel];
                    }
                }
                kernels->M2L( &upperCells[idxLevel] , neighbors, counter , idxLevel + 2);
            }

            memset(neighbors, 0, sizeof(CellClass*) * 343);
            counter = 0;
            for(int idxX = ifDir(DirMinusX,-2,0) ; idxX <= ifDir(DirPlusX,3,0) ; ++idxX){
                for(int idxY = ifDir(DirMinusY,-2,0) ; idxY <= ifDir(DirPlusY,3,0) ; ++idxY){
                    for(int idxZ = ifDir(DirMinusZ,-2,0) ; idxZ <= ifDir(DirPlusZ,3,0) ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            const int index = neighIndex(idxX,idxY,idxZ);
                            neighbors[index] = &positionedCells[index][1];
                            ++counter;
                        }
                    }
                }
            }
            kernels->M2L( &upperCells[1] , neighbors, 189, 3);
        }

        {
            if( usePerDir(AllDirs) ){
                CellClass leftborder, bottomborder, frontborder, angleborderlb,
                        angleborderfb, angleborderlf, angleborder;

                const CellClass* neighbors[343];
                memset(neighbors, 0, sizeof(CellClass*) * 343);
                int counter = 0;

                getUpperSpecialCell( &leftborder, upperCells[nbLevelsAboveRoot+1],     1,1 , 0,1 , 0,1);
                counter += setM2LVector(neighbors, leftborder,     -2,-2 , -1,1,  -1,1 );
                getUpperSpecialCell( &bottomborder, upperCells[nbLevelsAboveRoot+1],   0,1 , 0,1 , 1,1);
                counter += setM2LVector(neighbors, bottomborder,   -1,1  , -1,1,  -2,-2);
                getUpperSpecialCell( &frontborder, upperCells[nbLevelsAboveRoot+1],    0,1 , 1,1 , 0,1);
                counter += setM2LVector(neighbors, frontborder,    -1,1  , -2,-2, -1,1 );
                getUpperSpecialCell( &angleborderlb, upperCells[nbLevelsAboveRoot+1],  1,1 , 0,1 , 1,1);
                counter += setM2LVector(neighbors, angleborderlb,  -2,-2 , -1,1,  -2,-2);
                getUpperSpecialCell( &angleborderfb, upperCells[nbLevelsAboveRoot+1],  0,1 , 1,1 , 1,1);
                counter += setM2LVector(neighbors, angleborderfb,  -1,1 ,  -2,-2, -2,-2);
                getUpperSpecialCell( &angleborderlf, upperCells[nbLevelsAboveRoot+1],  1,1 , 1,1 , 0,1);
                counter += setM2LVector(neighbors, angleborderlf,  -2,-2 , -2,-2, -1,1 );
                getUpperSpecialCell( &angleborder, upperCells[nbLevelsAboveRoot+1],    1,1 , 1,1 , 1,1);
                counter += setM2LVector(neighbors, angleborder,    -2,-2 , -2,-2, -2,-2);


                kernels->M2L( &upperCells[0] , neighbors, counter, 2);

                CellClass* virtualChild[8];
                memset(virtualChild, 0, sizeof(CellClass*) * 8);
                virtualChild[childIndex(0,0,0)] = &upperCells[1];
                kernels->L2L( &upperCells[0], virtualChild, 2);
            }
            else {
                CellClass*const leftborder = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const bottomborder = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const frontborder = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const angleborderlb = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const angleborderfb = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const angleborderlf = new CellClass[nbLevelsAboveRoot+2];
                CellClass*const angleborder = new CellClass[nbLevelsAboveRoot+2];

                getUpperSpecialCells( leftborder,   upperCells[nbLevelsAboveRoot+1],     1,1 , 0,1 , 0,1);
                getUpperSpecialCells( bottomborder, upperCells[nbLevelsAboveRoot+1],   0,1 , 0,1 , 1,1);
                getUpperSpecialCells( frontborder,  upperCells[nbLevelsAboveRoot+1],    0,1 , 1,1 , 0,1);
                getUpperSpecialCells( angleborderlb,upperCells[nbLevelsAboveRoot+1],  1,1 , 0,1 , 1,1);
                getUpperSpecialCells( angleborderfb,upperCells[nbLevelsAboveRoot+1],  0,1 , 1,1 , 1,1);
                getUpperSpecialCells( angleborderlf,upperCells[nbLevelsAboveRoot+1],  1,1 , 1,1 , 0,1);
                getUpperSpecialCells( angleborder,  upperCells[nbLevelsAboveRoot+1],    1,1 , 1,1 , 1,1);

                const CellClass* neighbors[343];
                memset(neighbors, 0, sizeof(CellClass*) * 343);
                int counter = 0;

                if(usePerDir(DirMinusX) && usePerDir(DirMinusY) && usePerDir(DirMinusZ)){
                    neighbors[neighIndex(-2,-1,-1)] = &leftborder[0];
                    neighbors[neighIndex(-1,-2,-1)] = &frontborder[0];
                    neighbors[neighIndex(-1,-1,-2)] = &bottomborder[0];
                    neighbors[neighIndex(-2,-2,-1)] = &angleborderlf[0];
                    neighbors[neighIndex(-2,-1,-2)] = &angleborderlb[0];
                    neighbors[neighIndex(-1,-2,-2)] = &angleborderfb[0];
                    neighbors[neighIndex(-2,-2,-2)] = &angleborder[0];
                    counter += 7;
                }
                if(usePerDir(DirMinusX) && usePerDir(DirPlusY) && usePerDir(DirMinusZ)){
                    neighbors[neighIndex(-2,1,-1)] = &leftborder[0];
                    neighbors[neighIndex(-1,1,-2)] = &bottomborder[0];
                    neighbors[neighIndex(-2,1,-2)] = &angleborderlb[0];
                    counter += 3;
                }
                if(usePerDir(DirMinusX) && usePerDir(DirMinusY) && usePerDir(DirPlusZ)){
                    neighbors[neighIndex(-2,-1,1)] = &leftborder[0];
                    neighbors[neighIndex(-1,-2,1)] = &frontborder[0];
                    neighbors[neighIndex(-2,-2,1)] = &angleborderlf[0];
                    counter += 3;
                }
                if(usePerDir(DirMinusX) && usePerDir(DirPlusY) && usePerDir(DirPlusZ)){
                    neighbors[neighIndex(-2,1,1)] = &leftborder[0];
                    counter += 1;
                }
                if(usePerDir(DirPlusX) && usePerDir(DirMinusY) && usePerDir(DirMinusZ)){
                    neighbors[neighIndex(1,-2,-1)] = &frontborder[0];
                    neighbors[neighIndex(1,-1,-2)] = &bottomborder[0];
                    neighbors[neighIndex(1,-2,-2)] = &angleborderfb[0];
                    counter += 3;
                }
                if(usePerDir(DirPlusX) && usePerDir(DirMinusY) && usePerDir(DirPlusZ)){
                    neighbors[neighIndex(1,-2,1)] = &frontborder[0];
                    counter += 1;
                }
                if(usePerDir(DirPlusX) && usePerDir(DirPlusY) && usePerDir(DirMinusZ)){
                    neighbors[neighIndex(1,1,-2)] = &bottomborder[0];
                    counter += 1;
                }

                CellClass centerXFace;
                CellClass centerYFace;
                CellClass centerZFace;
                CellClass centerXZFace;
                CellClass centerXYFace;
                CellClass centerYXFace;
                CellClass centerYZFace;
                CellClass centerZXFace;
                CellClass centerZYFace;
                CellClass angleXZFace;
                CellClass angleXYFace;
                CellClass angleYZFace;

                if(usePerDir(DirMinusX)){
                    {
                        CellClass* virtualChild[8];
                        memset(virtualChild, 0, 8*sizeof(CellClass*));
                        if(usePerDir(DirPlusY) && usePerDir(DirPlusZ))  virtualChild[childIndex(1,1,1)] = &leftborder[1];
                        if(usePerDir(DirMinusY) && usePerDir(DirPlusZ)) virtualChild[childIndex(1,0,1)] = &leftborder[1];
                        else if( usePerDir(DirPlusZ))                   virtualChild[childIndex(1,0,1)] = &angleborderlf[1];
                        if(usePerDir(DirPlusY) && usePerDir(DirMinusZ)) virtualChild[childIndex(1,1,0)] = &leftborder[1];
                        else if(usePerDir(DirPlusY) )                   virtualChild[childIndex(1,1,0)] = &angleborderlb[1];

                        if(usePerDir(DirMinusY) && usePerDir(DirMinusZ)) virtualChild[childIndex(1,0,0)] = &leftborder[1];
                        else if(usePerDir(DirMinusZ))                    virtualChild[childIndex(1,0,0)] = &angleborderlf[1];
                        else if(usePerDir(DirMinusY))                    virtualChild[childIndex(1,0,0)] = &angleborderlb[1];
                        else                                             virtualChild[childIndex(1,0,0)] = &angleborder[1];

                        kernels->M2M( &centerXFace, virtualChild, 2);
                        neighbors[neighIndex(-2,0,0)] = &centerXFace;
                        counter += 1;
                    }
                    if(usePerDir(DirMinusZ) || usePerDir(DirPlusZ)){
                        if(usePerDir(DirY)){
                            centerXZFace = leftborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusY)){
                                virtualChild[childIndex(1,1,0)] = &leftborder[1];
                                virtualChild[childIndex(1,1,1)] = &leftborder[1];
                            }
                            if(usePerDir(DirMinusY)){
                                virtualChild[childIndex(1,0,0)] = &leftborder[1];
                                virtualChild[childIndex(1,0,1)] = &leftborder[1];
                            }
                            else{
                                virtualChild[childIndex(1,0,0)] = &angleborderlf[1];
                                virtualChild[childIndex(1,0,1)] = &angleborderlf[1];
                            }
                            kernels->M2M( &centerXZFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusZ) ){
                            neighbors[neighIndex(-2,0,1)] = &centerXZFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusZ) ){
                            neighbors[neighIndex(-2,0,-1)] = &centerXZFace;
                            counter += 1;

                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusY)){
                                virtualChild[childIndex(1,1,1)] = &angleborderlb[1];
                            }
                            if(usePerDir(DirMinusY)){
                                virtualChild[childIndex(1,0,1)] = &angleborderlb[1];
                            }
                            else{
                                virtualChild[childIndex(1,0,1)] = &angleborder[1];
                            }
                            kernels->M2M( &angleXZFace, virtualChild, 2);

                            neighbors[neighIndex(-2,0,-2)] = &angleXZFace;
                            counter += 1;
                        }
                    }
                    if(usePerDir(DirMinusY) || usePerDir(DirPlusY)){
                        if(usePerDir(DirZ)){
                            centerXYFace = leftborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusZ)){
                                virtualChild[childIndex(1,0,1)] = &leftborder[1];
                                virtualChild[childIndex(1,1,1)] = &leftborder[1];
                            }
                            if(usePerDir(DirMinusZ)){
                                virtualChild[childIndex(1,0,0)] = &leftborder[1];
                                virtualChild[childIndex(1,1,0)] = &leftborder[1];
                            }
                            else{
                                virtualChild[childIndex(1,0,0)] = &angleborderlb[1];
                                virtualChild[childIndex(1,1,0)] = &angleborderlb[1];
                            }
                            kernels->M2M( &centerXYFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusY) ){
                            neighbors[neighIndex(-2,1,0)] = &centerXYFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusY) ){
                            neighbors[neighIndex(-2,-1,0)] = &centerXYFace;
                            counter += 1;

                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusZ)){
                                virtualChild[childIndex(1,1,1)] = &angleborderlf[1];
                            }
                            if(usePerDir(DirMinusZ)){
                                virtualChild[childIndex(1,1,0)] = &angleborderlf[1];
                            }
                            else{
                                virtualChild[childIndex(1,1,0)] = &angleborder[1];
                            }
                            kernels->M2M( &angleXYFace, virtualChild, 2);

                            neighbors[neighIndex(-2,-2,0)] = &angleXYFace;
                            counter += 1;
                        }
                    }
                }
                if(usePerDir(DirMinusY)){
                    {
                        CellClass* virtualChild[8];
                        memset(virtualChild, 0, 8*sizeof(CellClass*));
                        if(usePerDir(DirPlusX) && usePerDir(DirPlusZ))  virtualChild[childIndex(1,1,1)] = &frontborder[1];
                        if(usePerDir(DirMinusX) && usePerDir(DirPlusZ)) virtualChild[childIndex(0,1,1)] = &frontborder[1];
                        else if(usePerDir(DirPlusZ))                    virtualChild[childIndex(0,1,1)] = &angleborderlf[1];
                        if(usePerDir(DirPlusX) && usePerDir(DirMinusZ)) virtualChild[childIndex(1,1,0)] = &frontborder[1];
                        else if(usePerDir(DirPlusX))                    virtualChild[childIndex(1,1,0)] = &angleborderfb[1];

                        if(usePerDir(DirMinusX) && usePerDir(DirMinusZ)) virtualChild[childIndex(0,1,0)] = &frontborder[1];
                        else if(usePerDir(DirMinusZ))                    virtualChild[childIndex(0,1,0)] = &angleborderlf[1];
                        else if(usePerDir(DirMinusX))                    virtualChild[childIndex(0,1,0)] = &angleborderfb[1];
                        else                                             virtualChild[childIndex(0,1,0)] = &angleborder[1];

                        kernels->M2M( &centerYFace, virtualChild, 2);
                        neighbors[neighIndex(0,-2,0)] = &centerYFace;
                        counter += 1;
                    }
                    if(usePerDir(DirMinusZ) || usePerDir(DirPlusZ)){
                        if(usePerDir(DirX)){
                            centerYZFace = frontborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusX)){
                                virtualChild[childIndex(1,1,0)] = &frontborder[1];
                                virtualChild[childIndex(1,1,1)] = &frontborder[1];
                            }
                            if(usePerDir(DirMinusX)){
                                virtualChild[childIndex(0,1,0)] = &frontborder[1];
                                virtualChild[childIndex(0,1,1)] = &frontborder[1];
                            }
                            else{
                                virtualChild[childIndex(0,1,0)] = &angleborderlf[1];
                                virtualChild[childIndex(0,1,1)] = &angleborderlf[1];
                            }
                            kernels->M2M( &centerYZFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusZ) ){
                            neighbors[neighIndex(0,-2,1)] = &centerYZFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusZ) ){
                            neighbors[neighIndex(0,-2,-1)] = &centerYZFace;
                            counter += 1;

                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusX)){
                                virtualChild[childIndex(1,1,1)] = &angleborderfb[1];
                            }
                            if(usePerDir(DirMinusX)){
                                virtualChild[childIndex(0,1,1)] = &angleborderfb[1];
                            }
                            else{
                                virtualChild[childIndex(0,1,1)] = &angleborder[1];
                            }
                            kernels->M2M( &angleYZFace, virtualChild, 2);

                            neighbors[neighIndex(0,-2,-2)] = &angleYZFace;
                            counter += 1;
                        }
                    }
                    if(usePerDir(DirMinusX) || usePerDir(DirPlusX)){
                        if(usePerDir(DirZ)){
                            centerYXFace = frontborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusZ)){
                                virtualChild[childIndex(0,1,1)] = &frontborder[1];
                                virtualChild[childIndex(1,1,1)] = &frontborder[1];
                            }
                            if(usePerDir(DirMinusZ)){
                                virtualChild[childIndex(0,1,0)] = &frontborder[1];
                                virtualChild[childIndex(1,1,0)] = &frontborder[1];
                            }
                            else{
                                virtualChild[childIndex(0,1,0)] = &angleborderfb[1];
                                virtualChild[childIndex(1,1,0)] = &angleborderfb[1];
                            }
                            kernels->M2M( &centerYXFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusX) ){
                            neighbors[neighIndex(1,-2,0)] = &centerYXFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusX) ){
                            neighbors[neighIndex(-1,-2,0)] = &centerYXFace;
                            counter += 1;
                        }
                    }
                }
                if(usePerDir(DirMinusZ)){
                    {
                        CellClass* virtualChild[8];
                        memset(virtualChild, 0, 8*sizeof(CellClass*));
                        if(usePerDir(DirPlusX) && usePerDir(DirPlusY))  virtualChild[childIndex(1,1,1)] = &bottomborder[1];
                        if(usePerDir(DirMinusX) && usePerDir(DirPlusY)) virtualChild[childIndex(0,1,1)] = &bottomborder[1];
                        else if( usePerDir(DirPlusY))                   virtualChild[childIndex(0,1,1)] = &angleborderlb[1];
                        if(usePerDir(DirPlusX) && usePerDir(DirMinusY)) virtualChild[childIndex(1,0,1)] = &bottomborder[1];
                        else if(usePerDir(DirPlusX) )                   virtualChild[childIndex(1,0,1)] = &angleborderfb[1];

                        if(usePerDir(DirMinusX) && usePerDir(DirMinusY)) virtualChild[childIndex(0,0,1)] = &bottomborder[1];
                        else if(usePerDir(DirMinusY))                    virtualChild[childIndex(0,0,1)] = &angleborderfb[1];
                        else if(usePerDir(DirMinusX))                    virtualChild[childIndex(0,0,1)] = &angleborderlb[1];
                        else                                             virtualChild[childIndex(0,0,1)] = &angleborder[1];

                        kernels->M2M( &centerZFace, virtualChild, 2);
                        neighbors[neighIndex(0,0,-2)] = &centerZFace;
                        counter += 1;
                    }
                    if(usePerDir(DirMinusY) || usePerDir(DirPlusY)){
                        if(usePerDir(DirX)){
                            centerZYFace = bottomborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusX)){
                                virtualChild[childIndex(1,0,1)] = &bottomborder[1];
                                virtualChild[childIndex(1,1,1)] = &bottomborder[1];
                            }
                            if(usePerDir(DirMinusX)){
                                virtualChild[childIndex(0,0,1)] = &bottomborder[1];
                                virtualChild[childIndex(0,1,1)] = &bottomborder[1];
                            }
                            else{
                                virtualChild[childIndex(0,0,1)] = &angleborderlb[1];
                                virtualChild[childIndex(0,1,1)] = &angleborderlb[1];
                            }
                            kernels->M2M( &centerZYFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusY) ){
                            neighbors[neighIndex(0,1,-2)] = &centerZYFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusY) ){
                            neighbors[neighIndex(0,-1,-2)] = &centerZYFace;
                            counter += 1;
                        }
                    }
                    if(usePerDir(DirMinusX) || usePerDir(DirPlusX)){
                        if(usePerDir(DirY)){
                            centerZXFace = bottomborder[0];
                        }
                        else{
                            CellClass* virtualChild[8];
                            memset(virtualChild, 0, 8*sizeof(CellClass*));
                            if(usePerDir(DirPlusY)){
                                virtualChild[childIndex(0,1,1)] = &bottomborder[1];
                                virtualChild[childIndex(1,1,1)] = &bottomborder[1];
                            }
                            if(usePerDir(DirMinusY)){
                                virtualChild[childIndex(0,0,1)] = &bottomborder[1];
                                virtualChild[childIndex(1,0,1)] = &bottomborder[1];
                            }
                            else{
                                virtualChild[childIndex(0,0,1)] = &angleborderlf[1];
                                virtualChild[childIndex(0,1,1)] = &angleborderlf[1];
                            }
                            kernels->M2M( &centerZXFace, virtualChild, 2);
                        }
                        if( usePerDir(DirPlusX) ){
                            neighbors[neighIndex(1,0,-2)] = &centerZXFace;
                            counter += 1;
                        }
                        if( usePerDir(DirMinusX) ){
                            neighbors[neighIndex(-1,0,-2)] = &centerZXFace;
                            counter += 1;
                        }
                    }
                }

                for(int idxX = -3 ; idxX <= 3 ; ++idxX){
                    for(int idxY = -3 ; idxY <= 3 ; ++idxY){
                        for(int idxZ = -3 ; idxZ <= 3 ; ++idxZ){
                            if(neighbors[neighIndex(idxX,idxY,idxZ)]){
                                printf("Used at %d %d %d\n", idxX, idxY, idxZ);
                            }
                        }
                    }
                }

                kernels->M2L( &upperCells[0] , neighbors, counter, 2);

                CellClass* virtualChild[8];
                memset(virtualChild, 0, sizeof(CellClass*) * 8);
                virtualChild[childIndex(0,0,0)] = &upperCells[1];
                kernels->L2L( &upperCells[0], virtualChild, 2);

                delete[] leftborder;
                delete[] bottomborder;
                delete[] frontborder;
                delete[] angleborderlb;
                delete[] angleborderfb;
                delete[] angleborderlf;
                delete[] angleborder;
            }

        }

        // Finally L2L until level 0
        {
            CellClass* virtualChild[8];
            memset(virtualChild, 0, sizeof(CellClass*) * 8);
            for(int idxLevel = 1 ; idxLevel <= nbLevelsAboveRoot  ; ++idxLevel){
                virtualChild[childIndex(1,1,1)] = &upperCells[idxLevel+1];
                kernels->L2L( &upperCells[idxLevel], virtualChild, idxLevel + 2);
            }
        }

        // L2L from 0 to level 1
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoLeft();
            kernels->L2L( &upperCells[nbLevelsAboveRoot+1], octreeIterator.getCurrentBox(), offsetRealTree);
        }
    }

    void getUpperSpecialCell( CellClass*const result, const CellClass& root, const int startX,
                              const int endX, const int startY, const int endY, const int startZ,
                              const int endZ){
        CellClass*const cellsAtLevel = new CellClass[nbLevelsAboveRoot+2];
        getUpperSpecialCells(cellsAtLevel,root,startX,endX,startY,endY,startZ,endZ);
        *result = cellsAtLevel[0];
        delete[] cellsAtLevel;
    }

    void getUpperSpecialCells( CellClass cellsAtLevel[], const CellClass& root, const int startX,
                              const int endX, const int startY, const int endY, const int startZ,
                              const int endZ){
        cellsAtLevel[nbLevelsAboveRoot+1] = root;
        CellClass* virtualChild[8];
        for(int idxLevel = nbLevelsAboveRoot ; idxLevel >= 0  ; --idxLevel){
            memset(virtualChild, 0, sizeof(CellClass*)*8);
            for(int idxX = startX ; idxX <= endX ; ++idxX)
                for(int idxY = startY ; idxY <= endY ; ++idxY)
                    for(int idxZ = startZ ; idxZ <= endZ ; ++idxZ)
                        virtualChild[childIndex(idxX,idxY,idxZ)] = &cellsAtLevel[idxLevel+1];

            kernels->M2M( &cellsAtLevel[idxLevel], virtualChild, idxLevel + 2);
        }
    }

    int setM2LVector(const CellClass* neighbors[343], const CellClass& source,
                     const int startX, const int endX, const int startY, const int endY,
                     const int startZ, const int endZ){
        int counter = 0;
        for(int idxX = startX ; idxX <= endX ; ++idxX)
            for(int idxY = startY ; idxY <= endY ; ++idxY)
                for(int idxZ = startZ ; idxZ <= endZ ; ++idxZ){
                    neighbors[neighIndex(idxX,idxY,idxZ)] = &source;
                    ++counter;
                }
        return counter;
    }

    int childIndex(const int x, const int y, const int z) const {
        return (x<<2) | (y<<1) | z;
    }

    int neighIndex(const int x, const int y, const int z) const {
        return (((x+3)*7) + (y+3))*7 + (z + 3);
    }

};


#endif // FFMMALGORITHMPERIODIC_HPP
