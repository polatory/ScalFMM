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
        bool isPeriodicLeaf;
        // for each leafs
        do{
            FDEBUG(computationCounterL2P.tic());
            kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
            FDEBUG(computationCounterL2P.tac());

            // need the current particles and neighbors particles
            const FTreeCoordinate centerOfLeaf = octreeIterator.getCurrentGlobalCoordinate();
            const int counter = tree->getPeriodicLeafsNeighbors( neighbors, offsets, &isPeriodicLeaf, centerOfLeaf, heightMinusOne, periodicDirections);

            if(isPeriodicLeaf){
                for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
                    if( neighbors[idxNeig] ){
                        typename ContainerClass::BasicIterator iter(*neighbors[idxNeig]);
                        while( iter.hasNotFinished() ){
                            iter.data().incPosition(boxWidth * FReal(offsets[idxNeig].getX()),
                                                    boxWidth * FReal(offsets[idxNeig].getY()),
                                                    boxWidth * FReal(offsets[idxNeig].getZ()));
                            iter.gotoNext();
                        }
                    }
                }
            }

            FDEBUG(computationCounterP2P.tic());
            kernels->P2P(octreeIterator.getCurrentGlobalCoordinate(),octreeIterator.getCurrentListTargets(),
                         octreeIterator.getCurrentListSrc(), neighbors, counter);
            FDEBUG(computationCounterP2P.tac());

            if(isPeriodicLeaf){
                for(int idxNeig = 0 ; idxNeig < 27 ; ++idxNeig){
                    if( neighbors[idxNeig] ){
                        typename ContainerClass::BasicIterator iter(*neighbors[idxNeig]);
                        while( iter.hasNotFinished() ){
                            iter.data().incPosition(-boxWidth * FReal(offsets[idxNeig].getX()),
                                                    -boxWidth * FReal(offsets[idxNeig].getY()),
                                                    -boxWidth * FReal(offsets[idxNeig].getZ()));
                            iter.gotoNext();
                        }
                    }
                }
            }

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
        min->setPosition(ifDir(DirMinusX,halfRepeated,0),ifDir(DirMinusY,halfRepeated,0),
                         ifDir(DirMinusZ,halfRepeated,0));
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
        // If nb level == -1 nothing to do
        if( nbLevelsAboveRoot == -1 ){
            return;
        }
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
            for(int idxX = -3 ; idxX <= 3 ; ++idxX){
                for(int idxY = -3 ; idxY <= 3 ; ++idxY){
                    for(int idxZ = -3 ; idxZ <= 3 ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            neighbors[(((idxX+3)*7) + (idxY+3))*7 + (idxZ + 3)] = &rootUp;
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

        // in other situation, we have to compute M2L from 0 to nbLevelsAboveRoot
        // but also at nbLevelsAboveRoot +1 for the rest
        CellClass upperCells[nbLevelsAboveRoot+2];

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
        }
        // Then M2L at all level
        {
            // We say that we are in the child index 0
            // So we can compute one time the relative indexes
            const CellClass* neighbors[343];
            memset(neighbors, 0, sizeof(CellClass*) * 343);
            int counter = 0;
            for(int idxX = -3 ; idxX <= 2 ; ++idxX){
                for(int idxY = -3 ; idxY <= 2 ; ++idxY){
                    for(int idxZ = -3 ; idxZ <= 2 ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            neighbors[(((idxX+3)*7) + (idxY+3))*7 + (idxZ + 3)] = reinterpret_cast<const CellClass*>(~0);
                            ++counter;
                        }
                    }
                }
            }

            for(int idxLevel = nbLevelsAboveRoot + 1 ; idxLevel > 1 ; --idxLevel ){
                for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
                    if(neighbors[idxNeigh]){
                        neighbors[idxNeigh] = &upperCells[idxLevel];
                    }
                }
                kernels->M2L( &upperCells[idxLevel] , neighbors, counter , idxLevel + 2);
            }

            memset(neighbors, 0, sizeof(CellClass*) * 343);
            for(int idxX = -2 ; idxX <= 3 ; ++idxX){
                for(int idxY = -2 ; idxY <= 3 ; ++idxY){
                    for(int idxZ = -2 ; idxZ <= 3 ; ++idxZ){
                        if( FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            neighbors[(((idxX+3)*7) + (idxY+3))*7 + (idxZ + 3)] = &upperCells[1];
                        }
                    }
                }
            }
            kernels->M2L( &upperCells[1] , neighbors, 189, 3);
        }

        {
            CellClass leftborder, bottomborder, frontborder, angleborderlb,
                    angleborderfb, angleborderlf, angleborder;
            const CellClass* neighbors[343];
            memset(neighbors, 0, sizeof(CellClass*) * 343);
            int counter = 0;

            //if( usePerDir(DirMinusX) && usePerDir(DirMinusY) && usePerDir(DirMinusZ)){
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
            //}
            /*else {
                if(usePerDir(DirMinusX)){
                    if(usePerDir(DirMinusY)){
                        getUpperSpecialCell( &leftborder, upperCells[nbLevelsAboveRoot+1],     1,1 , 0,1 , 1,1);
                        counter += setM2LVector(neighbors, leftborder,     -2,-2 , -1,1,  0,0 );
                        getUpperSpecialCell( &frontborder, upperCells[nbLevelsAboveRoot+1],    0,1 , 1,1 , 1,1);
                        counter += setM2LVector(neighbors, frontborder,    -1,1  , -2,-2, 0,0 );
                        getUpperSpecialCell( &angleborderlf, upperCells[nbLevelsAboveRoot+1],  1,1 , 1,1 , 1,1);
                        counter += setM2LVector(neighbors, angleborderlf,  -2,-2 , -2,-2, 0,0 );
                    }
                    else {
                        if(usePerDir(DirMinusZ)){
                            getUpperSpecialCell( &leftborder, upperCells[nbLevelsAboveRoot+1],     1,1 , 1,1 , 0,1);
                            counter += setM2LVector(neighbors, leftborder,     -2,-2 , 0,0,  -1,1 );
                            getUpperSpecialCell( &bottomborder, upperCells[nbLevelsAboveRoot+1],   0,1 , 1,1 , 1,1);
                            counter += setM2LVector(neighbors, bottomborder,   -1,1  , 0,0,  -2,-2);
                            getUpperSpecialCell( &angleborderlb, upperCells[nbLevelsAboveRoot+1],  1,1 , 1,1 , 1,1);
                            counter += setM2LVector(neighbors, angleborderlb,  -2,-2 , 0,0,  -2,-2);
                        }
                        else{
                            getUpperSpecialCell( &leftborder, upperCells[nbLevelsAboveRoot+1],     1,1 , 1,1 , 1,1);
                            counter += setM2LVector(neighbors, leftborder,     -2,-2 , 0,0,  0,0 );
                        }
                    }
                }
                else if(usePerDir(DirMinusY)){
                    if(usePerDir(DirMinusZ)){
                        getUpperSpecialCell( &bottomborder, upperCells[nbLevelsAboveRoot+1],   1,1 , 0,1 , 1,1);
                        counter += setM2LVector(neighbors, bottomborder,   0,0  , -1,1,  -2,-2);
                        getUpperSpecialCell( &frontborder, upperCells[nbLevelsAboveRoot+1],    1,1 , 1,1 , 0,1);
                        counter += setM2LVector(neighbors, frontborder,    0,0  , -2,-2, -1,1 );
                        getUpperSpecialCell( &angleborderfb, upperCells[nbLevelsAboveRoot+1],  1,1 , 1,1 , 1,1);
                        counter += setM2LVector(neighbors, angleborderfb,  0,0 ,  -2,-2, -2,-2);
                    }
                    else{
                        getUpperSpecialCell( &frontborder, upperCells[nbLevelsAboveRoot+1],    1,1 , 1,1 , 1,1);
                        counter += setM2LVector(neighbors, frontborder,    0,0  , -2,-2, 0,0 );
                    }
                }
                else if(usePerDir(DirMinusZ)){
                    getUpperSpecialCell( &bottomborder, upperCells[nbLevelsAboveRoot+1],   1,1 , 1,1 , 1,1);
                    counter += setM2LVector(neighbors, bottomborder,   0,0  , 0,0,  -2,-2);
                }
            }*/

            kernels->M2L( &upperCells[0] , neighbors, counter, 2);

            CellClass* virtualChild[8];
            memset(virtualChild, 0, sizeof(CellClass*) * 8);
            virtualChild[0] = &upperCells[1];
            kernels->L2L( &upperCells[0], virtualChild, 2);
        }

        // Finally L2L until level 0
        {
            CellClass* virtualChild[8];
            memset(virtualChild, 0, sizeof(CellClass*) * 8);
            for(int idxLevel = 1 ; idxLevel <= nbLevelsAboveRoot  ; ++idxLevel){
                virtualChild[7] = &upperCells[idxLevel+1];
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
        CellClass cellsAtLevel[nbLevelsAboveRoot+2];
        cellsAtLevel[nbLevelsAboveRoot+1] = root;
        CellClass* virtualChild[8];
        for(int idxLevel = nbLevelsAboveRoot ; idxLevel >= 0  ; --idxLevel){
            memset(virtualChild, 0, sizeof(CellClass*)*8);
            for(int idxX = startX ; idxX <= endX ; ++idxX)
                for(int idxY = startY ; idxY <= endY ; ++idxY)
                    for(int idxZ = startZ ; idxZ <= endZ ; ++idxZ)
                        virtualChild[(idxX<<2) | (idxY<<1) | idxZ] = &cellsAtLevel[idxLevel+1];

            kernels->M2M( &cellsAtLevel[idxLevel], virtualChild, idxLevel + 2);
        }
        *result = cellsAtLevel[0];
    }

    int setM2LVector(const CellClass* neighbors[343], const CellClass& source,
                     const int startX, const int endX, const int startY, const int endY,
                     const int startZ, const int endZ){
        int counter = 0;
        for(int idxX = startX ; idxX <= endX ; ++idxX)
            for(int idxY = startY ; idxY <= endY ; ++idxY)
                for(int idxZ = startZ ; idxZ <= endZ ; ++idxZ){
                    neighbors[(((idxX+3)*7) + (idxY+3))*7 + (idxZ + 3)] = &source;
                    ++counter;
                }
        return counter;
    }

};


#endif // FFMMALGORITHMPERIODIC_HPP
