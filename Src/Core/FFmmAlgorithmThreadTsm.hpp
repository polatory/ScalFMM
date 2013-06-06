// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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
#ifndef FFMMALGORITHMTHREADTSM_HPP
#define FFMMALGORITHMTHREADTSM_HPP


#include "../Utils/FAssertable.hpp"
#include "../Utils/FDebug.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"
#include "../Utils/FGlobal.hpp"

#include "../Containers/FOctree.hpp"


#include <omp.h>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmmAlgorithmThreadTsm
* @brief
* Please read the license
*
* This class is a threaded FMM algorithm
* It just iterates on a tree and call the kernels with good arguments.
* It used the inspector-executor model :
* iterates on the tree and builds an array to work in parallel on this array
*
* Of course this class does not deallocate pointer given in arguements.
*
* Because this is a Target source model you do not need the P2P to be safe.
* You should not write on sources in the P2P method!
*/
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmThreadTsm : protected FAssertable{
    OctreeClass* const tree;                  //< The octree to work on
    KernelClass** kernels;                    //< The kernels

    typename OctreeClass::Iterator* iterArray;

    const int MaxThreads;

    const int OctreeHeight;

public:
    /** The constructor need the octree and the kernels used for computation
      * @param inTree the octree to work on
      * @param inKernels the kernels to call
      * An assert is launched if one of the arguments is null
      */
    FFmmAlgorithmThreadTsm(OctreeClass* const inTree, KernelClass* const inKernels)
                      : tree(inTree) , kernels(0), iterArray(0),
                      MaxThreads(omp_get_max_threads()) , OctreeHeight(tree->getHeight()) {

        fassert(tree, "tree cannot be null", __LINE__, __FILE__);

        this->kernels = new KernelClass*[MaxThreads];
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            this->kernels[idxThread] = new KernelClass(*inKernels);
        }

        FDEBUG(FDebug::Controller << "FFmmAlgorithmThreadTsm\n");
    }

    /** Default destructor */
    virtual ~FFmmAlgorithmThreadTsm(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;
    }

    /**
      * To execute the fmm algorithm
      * Call this function to run the complete algorithm
      */
    void execute(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

        // Count leaf
        int numberOfLeafs = 0;
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            ++numberOfLeafs;
        } while(octreeIterator.moveRight());
        iterArray = new typename OctreeClass::Iterator[numberOfLeafs];
        fassert(iterArray, "iterArray bad alloc", __LINE__, __FILE__);

        bottomPass();

        upwardPass();

        transferPass();

        downardPass();

        directPass();

        delete [] iterArray;
        iterArray = 0;


    }

    /** P2M */
    void bottomPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Bottom Pass\n").write(FDebug::Flush) );
        FDEBUG( FTic counterTime );

        typename OctreeClass::Iterator octreeIterator(tree);
        int numberOfLeafs = 0;
        // Iterate on leafs
        octreeIterator.gotoBottomLeft();
        do{
            iterArray[numberOfLeafs] = octreeIterator;
            ++numberOfLeafs;
        } while(octreeIterator.moveRight());

        FDEBUG(FTic computationCounter);
        #pragma omp parallel
        {
            KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
            #pragma omp for nowait
            for(int idxLeafs = 0 ; idxLeafs < numberOfLeafs ; ++idxLeafs){
                // We need the current cell that represent the leaf
                // and the list of particles
                ContainerClass* const sources = iterArray[idxLeafs].getCurrentListSrc();
                if(sources->getNbParticles()){
                    iterArray[idxLeafs].getCurrentCell()->setSrcChildTrue();
                    myThreadkernels->P2M( iterArray[idxLeafs].getCurrentCell() , sources);
                }
                if(iterArray[idxLeafs].getCurrentListTargets()->getNbParticles()){
                    iterArray[idxLeafs].getCurrentCell()->setTargetsChildTrue();
                }
            }
        }
        FDEBUG(computationCounter.tac());

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.elapsed() << " s\n" );

    }

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
            FDEBUG(FTic counterTimeLevel);
            int numberOfCells = 0;
            // for each cells
            do{
                iterArray[numberOfCells] = octreeIterator;
                ++numberOfCells;
            } while(octreeIterator.moveRight());
            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();

            FDEBUG(computationCounter.tic());
            #pragma omp parallel
            {
                KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                #pragma omp for nowait
                for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                    // We need the current cell and the child
                    // child is an array (of 8 child) that may be null
                    CellClass* potentialChild[8];
                    CellClass** const realChild = iterArray[idxCell].getCurrentChild();
                    CellClass* const currentCell = iterArray[idxCell].getCurrentCell();
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
                    myThreadkernels->M2M( currentCell , potentialChild, idxLevel);
                }
            }
            FDEBUG(computationCounter.tac());
            FDEBUG( FDebug::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
        }

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.elapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );

    }

    /** M2L */
    void transferPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (M2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            FDEBUG(FTic computationCounter);

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.moveDown();
            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            // for each levels
            for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
                FDEBUG(FTic counterTimeLevel);
                int numberOfCells = 0;
                // for each cells
                do{
                    iterArray[numberOfCells] = octreeIterator;
                    ++numberOfCells;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                FDEBUG(computationCounter.tic());
                #pragma omp parallel
                {
                    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                    const CellClass* neighbors[343];

                    #pragma omp for  schedule(dynamic) nowait
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        CellClass* const currentCell = iterArray[idxCell].getCurrentCell();
                        if(currentCell->hasTargetsChild()){
                            const int counter = tree->getInteractionNeighbors(neighbors, iterArray[idxCell].getCurrentGlobalCoordinate(),idxLevel);
                            if( counter ){
                                int counterWithSrc = 0;
                                for(int idxRealNeighbors = 0 ; idxRealNeighbors < 343 ; ++idxRealNeighbors ){
                                    if(neighbors[idxRealNeighbors] && neighbors[idxRealNeighbors]->hasSrcChild()){
                                        ++counterWithSrc;
                                    }
                                    else{
                                        neighbors[idxRealNeighbors] = 0;
                                    }
                                }
                                if(counterWithSrc){
                                    myThreadkernels->M2L( currentCell , neighbors, counterWithSrc, idxLevel);
                                }
                            }
                        }
                    }

                    FDEBUG(computationCounter.tic());
                    myThreadkernels->finishedLevelM2L(idxLevel);
                    FDEBUG(computationCounter.tac());
                }
                FDEBUG(computationCounter.tac());
                FDEBUG( FDebug::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
            }
            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        }

        /* L2L */
        void downardPass(){
            FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

            FDEBUG( FDebug::Controller.write("\tStart Downward Pass (L2L)\n").write(FDebug::Flush); );
            FDEBUG(FTic counterTime);
            FDEBUG(FTic computationCounter);

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.moveDown();

            typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

            const int heightMinusOne = OctreeHeight - 1;
            // for each levels exepted leaf level
            for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
                FDEBUG(FTic counterTimeLevel);
                int numberOfCells = 0;
                // for each cells
                do{
                    iterArray[numberOfCells] = octreeIterator;
                    ++numberOfCells;
                } while(octreeIterator.moveRight());
                avoidGotoLeftIterator.moveDown();
                octreeIterator = avoidGotoLeftIterator;

                FDEBUG(computationCounter.tic());
                #pragma omp parallel
                {
                    KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
                    #pragma omp for nowait
                    for(int idxCell = 0 ; idxCell < numberOfCells ; ++idxCell){
                        CellClass* potentialChild[8];
                        CellClass** const realChild = iterArray[idxCell].getCurrentChild();
                        CellClass* const currentCell = iterArray[idxCell].getCurrentCell();
                        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                            if(realChild[idxChild] && realChild[idxChild]->hasTargetsChild()){
                                potentialChild[idxChild] = realChild[idxChild];
                            }
                            else{
                                potentialChild[idxChild] = 0;
                            }
                        }
                        myThreadkernels->L2L( currentCell , potentialChild, idxLevel);
                    }
                }
                FDEBUG(computationCounter.tac());
                FDEBUG( FDebug::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
            }
            FDEBUG( FDebug::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
            FDEBUG( FDebug::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
        }


    /** P2P */
    void directPass(){
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
        FDEBUG( FDebug::Controller.write("\tStart Direct Pass\n").write(FDebug::Flush); );
        FDEBUG(FTic counterTime);

        int numberOfLeafs = 0;
        {
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            // for each leaf
            do{
                iterArray[numberOfLeafs] = octreeIterator;
                ++numberOfLeafs;
            } while(octreeIterator.moveRight());
        }

        const int heightMinusOne = OctreeHeight - 1;
        FDEBUG(FTic computationCounter);
        #pragma omp parallel
        {
            KernelClass * const myThreadkernels = kernels[omp_get_thread_num()];
            // There is a maximum of 26 neighbors
            ContainerClass* neighbors[27];

            #pragma omp for schedule(dynamic) nowait
            for(int idxLeafs = 0 ; idxLeafs < numberOfLeafs ; ++idxLeafs){
                myThreadkernels->L2P(iterArray[idxLeafs].getCurrentCell(), iterArray[idxLeafs].getCurrentListTargets());
                // need the current particles and neighbors particles
                const int counter = tree->getLeafsNeighbors(neighbors, iterArray[idxLeafs].getCurrentGlobalCoordinate(),heightMinusOne);
                myThreadkernels->P2P( iterArray[idxLeafs].getCurrentGlobalCoordinate(), iterArray[idxLeafs].getCurrentListTargets(),
                                      iterArray[idxLeafs].getCurrentListSrc() , neighbors, counter);
            }
        }
        FDEBUG(computationCounter.tac());

        FDEBUG( counterTime.tac() );
        FDEBUG( FDebug::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FDEBUG( FDebug::Controller << "\t\t Computation L2P + P2P : " << computationCounter.elapsed() << " s\n" );

    }

};


#endif //FFMMALGORITHMTHREADTSM_HPP


