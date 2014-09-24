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
#ifndef FFMMALGORITHMTHREAD_HPP
#define FFMMALGORITHMTHREAD_HPP


#include "../Utils/FAssert.hpp"
#include "../Utils/FLog.hpp"

#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"

#include "../Containers/FOctree.hpp"

#include "FCoreCommon.hpp"

#include <omp.h>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmThread
* @brief
* Please read the license
*
* This class is a threaded FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
* It used the inspector-executor model :
* iterates on the tree and builds an array to work in parallel on this array
*
* Of course this class does not deallocate pointer given in arguments.
*
* When using this algorithm the P2P is thread safe.
*/
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmThread : public FAbstractAlgorithm{
    OctreeClass* const tree;                  //< The octree to work on
    KernelClass** kernels;                    //< The kernels

    typename OctreeClass::Iterator* iterArray;
    int leafsNumber;

    static const int SizeShape = 3*3*3;
    int shapeLeaf[SizeShape];    

    const int MaxThreads;

    const int OctreeHeight;

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThread(OctreeClass* const inTree, KernelClass* const inKernels)
        : tree(inTree) , kernels(nullptr), iterArray(nullptr), leafsNumber(0),
          MaxThreads(omp_get_max_threads()), OctreeHeight(tree->getHeight()) {

        FAssertLF(tree, "tree cannot be null");

        this->kernels = new KernelClass*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass(*inKernels);
        }

        FLOG(FLog::Controller << "FFmmAlgorithmThread (Max Thread " << omp_get_max_threads() << ")\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmThread(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(const unsigned operationsToProceed = FFmmNearAndFarFields){

        for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
            this->shapeLeaf[idxShape] = 0;
        }

        // Count leaf
        leafsNumber = 0;
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            ++leafsNumber;
            const FTreeCoordinate& coord = octreeIterator.getCurrentCell()->getCoordinate();
            ++this->shapeLeaf[(coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3)];

        } while(octreeIterator.moveRight());
        iterArray = new typename OctreeClass::Iterator[leafsNumber];
        FAssertLF(iterArray, "iterArray bad alloc");

        if(operationsToProceed & FFmmP2M) bottomPass();

        if(operationsToProceed & FFmmM2M) upwardPass();

        if(operationsToProceed & FFmmM2L) transferPass();

        if(operationsToProceed & FFmmL2L) downardPass();

        if((operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P)) directPass();

        delete [] iterArray;
        iterArray = nullptr;
    }

private:
    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** P2M */
    void bottomPass(){
        FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
        FLOG(FTic counterTime);

        typename OctreeClass::Iterator octreeIterator(tree);
        int leafs = 0;
        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            iterArray[leafs] = octreeIterator;
            ++leafs;
        } while(octreeIterator.moveRight());

        const int chunkSize = FMath::Max(1 , leafs/(omp_get_max_threads()*omp_get_max_threads()));

        FLOG(FTic computationCounter);
        #pragma omp parallel
        {
            KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for nowait schedule(dynamic, chunkSize)
            for(int idxLeafs = 0 ; idxLeafs < leafs ; ++idxLeafs){
                // We need the current cell that represent the leaf
                // and the list of particles
                myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , iterArray[idxLeafs].getCurrentListSrc());
            }
        }
        FLOG(computationCounter.tac() );

        FLOG( FLog::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );

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
        for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
            FLOG(FTic counterTimeLevel);
            int numberOfCells = 0;
            // for each cells
            do{
                iterArray[numberOfCells] = octreeIterator;
                ++numberOfCells;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();

            const int chunkSize = FMath::Max(1 , numberOfCells/(omp_get_max_threads()*omp_get_max_threads()));

            FLOG(computationCounter.tic());
            #pragma omp parallel
            {
                KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                #pragma omp for nowait  schedule(dynamic, chunkSize)
                for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                    // We need the current cell and the child
                    // child is an array (of 8 child) that may be null
                    myThreadkernels->M2M( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
                }
            }

            FLOG(computationCounter.tac());
            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
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
        octreeIterator.moveDown();
        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
            FLOG(FTic counterTimeLevel);
            int numberOfCells = 0;
            // for each cells
            do{
                iterArray[numberOfCells] = octreeIterator;
                ++numberOfCells;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;

            const int chunkSize = FMath::Max(1 , numberOfCells/(omp_get_max_threads()*omp_get_max_threads()));

            FLOG(computationCounter.tic());
            #pragma omp parallel
            {
                KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                const CellClass* neighbors[343];

                #pragma omp for  schedule(dynamic, chunkSize) nowait
                for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                    const int counter = tree->getInteractionNeighbors(neighbors,  iterArray[idxCell].getCurrentGlobalCoordinate(),idxLevel);
                    if(counter) myThreadkernels->M2L( iterArray[idxCell].getCurrentCell() , neighbors, counter, idxLevel);
                }

                myThreadkernels->finishedLevelM2L(idxLevel);
            }
            FLOG(computationCounter.tac());
            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
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
        FLOG(FTic computationCounter);

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.moveDown();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = OctreeHeight - 1;
        // for each levels excepted leaf level
        for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
            FLOG(FTic counterTimeLevel);
            int numberOfCells = 0;
            // for each cells
            do{
                iterArray[numberOfCells] = octreeIterator;
                ++numberOfCells;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;

            FLOG(computationCounter.tic());
            const int chunkSize = FMath::Max(1 , numberOfCells/(omp_get_max_threads()*omp_get_max_threads()));
            #pragma omp parallel
            {
                KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                #pragma omp for nowait schedule(dynamic, chunkSize)
                for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                    myThreadkernels->L2L( iterArray[idxCell].getCurrentCell() , iterArray[idxCell].getCurrentChild(), idxLevel);
                }
            }
            FLOG(computationCounter.tac());
            FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
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
        FLOG(FTic computationCounter);
        FLOG(FTic computationCounterP2P);

        omp_lock_t lockShape[SizeShape];
        for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
            omp_init_lock(&lockShape[idxShape]);
        }

        struct LeafData{
            MortonIndex index;
            CellClass* cell;
            ContainerClass* targets;
            ContainerClass* sources;
        };
        LeafData* const leafsDataArray = new LeafData[this->leafsNumber];

        const int LeafIndex = OctreeHeight - 1;

        int startPosAtShape[SizeShape];
        startPosAtShape[0] = 0;
        for(int idxShape = 1 ; idxShape < SizeShape ; ++idxShape){
            startPosAtShape[idxShape] = startPosAtShape[idxShape-1] + this->shapeLeaf[idxShape-1];
        }

        #pragma omp parallel
        {

            const float step = float(this->leafsNumber) / float(omp_get_num_threads());
            const int start = int(FMath::Ceil(step * float(omp_get_thread_num())));
            const int tempEnd = int(FMath::Ceil(step * float(omp_get_thread_num()+1)));
            const int end = (tempEnd > this->leafsNumber ? this->leafsNumber : tempEnd);

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();

            for(int idxPreLeaf = 0 ; idxPreLeaf < start ; ++idxPreLeaf){
                octreeIterator.moveRight();
            }

            // for each leafs
            for(int idxMyLeafs = start ; idxMyLeafs < end ; ++idxMyLeafs){
                //iterArray[leafs] = octreeIterator;
                //++leafs;
                const FTreeCoordinate& coord = octreeIterator.getCurrentGlobalCoordinate();
                const int shapePosition = (coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3);

                omp_set_lock(&lockShape[shapePosition]);
                const int positionToWork = startPosAtShape[shapePosition]++;
                omp_unset_lock(&lockShape[shapePosition]);

                leafsDataArray[positionToWork].index   = octreeIterator.getCurrentGlobalIndex();
                leafsDataArray[positionToWork].cell    = octreeIterator.getCurrentCell();
                leafsDataArray[positionToWork].targets = octreeIterator.getCurrentListTargets();
                leafsDataArray[positionToWork].sources = octreeIterator.getCurrentListSrc();

                octreeIterator.moveRight();
            }

            #pragma omp barrier

            FLOG(if(!omp_get_thread_num()) computationCounter.tic());

            KernelClass& myThreadkernels = (*kernels[omp_get_thread_num()]);
            // There is a maximum of 26 neighbors
            ContainerClass* neighbors[27];
            int previous = 0;

            for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
                const int endAtThisShape = this->shapeLeaf[idxShape] + previous;
                const int chunkSize = FMath::Max(1 , endAtThisShape/(omp_get_num_threads()*omp_get_num_threads()));
                #pragma omp for schedule(dynamic, chunkSize)
                for(int idxLeafs = previous ; idxLeafs < endAtThisShape ; ++idxLeafs){
                    LeafData& currentIter = leafsDataArray[idxLeafs];
                    myThreadkernels.L2P(currentIter.cell, currentIter.targets);
                    // need the current particles and neighbors particles
                    FLOG(if(!omp_get_thread_num()) computationCounterP2P.tic());
                    const int counter = tree->getLeafsNeighbors(neighbors, currentIter.cell->getCoordinate(), LeafIndex);
                    myThreadkernels.P2P(currentIter.cell->getCoordinate(), currentIter.targets,
                                        currentIter.sources, neighbors, counter);
                    FLOG(if(!omp_get_thread_num()) computationCounterP2P.tac());
                }

                previous = endAtThisShape;
            }
        }

        FLOG(computationCounter.tac());

        delete [] leafsDataArray;
        for(int idxShape = 0 ; idxShape < SizeShape ; ++idxShape){
            omp_destroy_lock(&lockShape[idxShape]);
        }


        FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation L2P + P2P : " << computationCounter.cumulated()    << " s\n" );
        FLOG( FLog::Controller << "\t\t Computation P2P :       " << computationCounterP2P.cumulated() << " s\n" );

    }

};


#endif //FFMMALGORITHMTHREAD_HPP


