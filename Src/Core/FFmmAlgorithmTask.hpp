// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FFMMALGORITHMTASK_HPP
#define FFMMALGORITHMTASK_HPP


#include "../Utils/FGlobal.hpp"
#include "../Utils/FAssert.hpp"
#include "../Utils/FLog.hpp"

#include "../Utils/FTic.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"

#include "FCoreCommon.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmTask
* @brief
* Please read the license
*
* This class is a basic FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
*
* Of course this class does not deallocate pointer given in arguements.
*/
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmTask : public FAbstractAlgorithm{

    OctreeClass* const tree;       //< The octree to work on
    KernelClass** kernels;    //< The kernels

    const int MaxThreads;

    const int OctreeHeight;

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmTask(OctreeClass* const inTree, KernelClass* const inKernels)
        : tree(inTree) , kernels(nullptr),
          MaxThreads(omp_get_max_threads()), OctreeHeight(tree->getHeight())
    {

        FAssertLF(tree, "tree cannot be null");
        FAssertLF(inKernels, "kernels cannot be null");

        this->kernels = new KernelClass*[MaxThreads];
        #pragma omp parallel for schedule(static)
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            #pragma omp critical (InitFFmmAlgorithmTask)
            {
                this->kernels[idxThread] = new KernelClass(*inKernels);
            }
        }

        FAbstractAlgorithm::setNbLevelsInTree(tree->getHeight());

        FLOG(FLog::Controller << "FFmmAlgorithmTask (Max Thread " << omp_get_max_threads() << ")\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmTask(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;
    }

protected:
    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void executeCore(const unsigned operationsToProceed) override {

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if((operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P)) directPass((operationsToProceed & FFmmP2P),(operationsToProceed & FFmmL2P));
    }

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M */
    void bottomPass(){
        FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
        FLOG(FTic counterTime);

        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                typename OctreeClass::Iterator octreeIterator(tree);

                // Iterate on leafs
                octreeIterator.gotoBottomLeft();
                do{
                    // We need the current cell that represent the leaf
                    // and the list of particles
                    #pragma omp task firstprivate(octreeIterator)
                    {
                        kernels[omp_get_thread_num()]->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
                    }
                } while(octreeIterator.moveRight());

                #pragma omp taskwait
            }
        }

        FLOG( FLog::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////

    /** M2M */
    void upwardPass(){
        FLOG( FLog::Controller.write("\tStart Upward Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);

        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                // Start from leal level - 1
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoBottomLeft();
                octreeIterator.moveUp();

                for(int idxLevel = OctreeHeight - 2 ; idxLevel > FAbstractAlgorithm::lowerWorkingLevel-2 ; --idxLevel){
                    octreeIterator.moveUp();
                }

                typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

                // for each levels
                for(int idxLevel = FAbstractAlgorithm::lowerWorkingLevel - 2 ; idxLevel >= FAbstractAlgorithm::upperWorkingLevel ; --idxLevel ){
                    FLOG(FTic counterTimeLevel);
                    // for each cells
                    do{
                        // We need the current cell and the child
                        // child is an array (of 8 child) that may be null
                        #pragma omp task firstprivate(octreeIterator) shared(idxLevel)
                        {
                            kernels[omp_get_thread_num()]->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                        }
                    } while(octreeIterator.moveRight());

                    avoidGotoLeftIterator.moveUp();
                    octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();

                    #pragma omp taskwait
                    FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
                }
            }
        }


        FLOG( FLog::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Transfer
    /////////////////////////////////////////////////////////////////////////////

    /** M2L L2L */
    void transferPass(){

        FLOG( FLog::Controller.write("\tStart Downward Pass (M2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                const CellClass* neighbors[343];

                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.moveDown();

                for(int idxLevel = 2 ; idxLevel < FAbstractAlgorithm::upperWorkingLevel ; --idxLevel){
                    octreeIterator.moveDown();
                }

                typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

                // for each levels
                for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < FAbstractAlgorithm::lowerWorkingLevel ; ++idxLevel ){
                    FLOG(FTic counterTimeLevel);
                    // for each cells
                    do{
                        int counter = tree->getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);
                        if(counter){
                            #pragma omp task firstprivate(octreeIterator, neighbors, counter) shared(idxLevel)
                            {
                                kernels[omp_get_thread_num()]->M2L( octreeIterator.getCurrentCell() , neighbors, counter, idxLevel);
                            }
                        }

                    } while(octreeIterator.moveRight());

                    avoidGotoLeftIterator.moveDown();
                    octreeIterator = avoidGotoLeftIterator;

                    #pragma omp taskwait

                    for( int idxThread = 0 ; idxThread < omp_get_num_threads() ; ++idxThread){
                        #pragma omp task
                        {
                            kernels[idxThread]->finishedLevelM2L(idxLevel);
                        }
                    }

                    #pragma omp taskwait
                    FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
                }
            }
        }
        FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downward
    /////////////////////////////////////////////////////////////////////////////

    void downardPass(){ // second L2L
        FLOG( FLog::Controller.write("\tStart Downward Pass (L2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);

        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.moveDown();

                for(int idxLevel = 2 ; idxLevel < FAbstractAlgorithm::upperWorkingLevel ; --idxLevel){
                    octreeIterator.moveDown();
                }

                typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

                const int heightMinusOne = FAbstractAlgorithm::lowerWorkingLevel - 1;
                // for each levels exepted leaf level
                for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; idxLevel < heightMinusOne ; ++idxLevel ){
                    FLOG(FTic counterTimeLevel);
                    // for each cells
                    do{
                        #pragma omp task firstprivate(octreeIterator) shared(idxLevel)
                        {
                            kernels[omp_get_thread_num()]->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
                        }

                    } while(octreeIterator.moveRight());

                    avoidGotoLeftIterator.moveDown();
                    octreeIterator = avoidGotoLeftIterator;

                    #pragma omp taskwait
                    FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
                }
            }
        }

        FLOG( FLog::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
    }


    /////////////////////////////////////////////////////////////////////////////
    // Direct
    /////////////////////////////////////////////////////////////////////////////

    /** P2P */
    void directPass(const bool p2pEnabled, const bool l2pEnabled){
        FLOG( FLog::Controller.write("\tStart Direct Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        const int heightMinusOne = OctreeHeight - 1;

        #pragma omp parallel
        {

            #pragma omp single nowait
            {
                // There is a maximum of 26 neighbors
                ContainerClass* neighbors[27];

                const int SizeShape = 3*3*3;
                FVector<typename OctreeClass::Iterator> shapes[SizeShape];

                typename OctreeClass::Iterator octreeIterator(tree);
                octreeIterator.gotoBottomLeft();

                // for each leafs
                do{
                    const FTreeCoordinate& coord = octreeIterator.getCurrentGlobalCoordinate();
                    const int shapePosition = (coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3);

                    shapes[shapePosition].push(octreeIterator);

                } while(octreeIterator.moveRight());

                FLOG( computationCounter.tic() );

                for( int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
                    const int nbLeaf = shapes[idxShape].getSize();
                    for(int iterLeaf = 0 ; iterLeaf < nbLeaf ; ++iterLeaf ){
                        typename OctreeClass::Iterator toWork = shapes[idxShape][iterLeaf];
                        #pragma omp task firstprivate(neighbors, toWork, l2pEnabled, p2pEnabled)
                        {
                            if(l2pEnabled){
                                kernels[omp_get_thread_num()]->L2P(toWork.getCurrentCell(), toWork.getCurrentListTargets());
                            }
                            if(p2pEnabled){
                                const int counter = tree->getLeafsNeighbors(neighbors, toWork.getCurrentGlobalCoordinate(),heightMinusOne);
                                kernels[omp_get_thread_num()]->P2P(toWork.getCurrentGlobalCoordinate(), toWork.getCurrentListTargets(),
                                                 toWork.getCurrentListSrc(), neighbors, counter);
                            }
                        }
                    }

                    #pragma omp taskwait
                }

                FLOG( computationCounter.tac() );
            }
        }


        FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation L2P + P2P : " << computationCounter.cumulated() << " s\n" );
    }

};


#endif //FFMMALGORITHMTASK_HPP


