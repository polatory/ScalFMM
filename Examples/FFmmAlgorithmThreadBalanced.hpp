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
#ifndef FFMMALGORITHMTHREADBALANCED_HPP
#define FFMMALGORITHMTHREADBALANCED_HPP


#include "../Src/Utils/FAssert.hpp"
#include "../Src/Utils/FLog.hpp"

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FGlobal.hpp"
#include "../Src/Utils/FAlgorithmTimers.hpp"

#include "../Src/Containers/FOctree.hpp"

#include "../Src/Core/FCoreCommon.hpp"

#include <vector>

#include <omp.h>

/**
 * \brief Implements a threaded FMM algorithm using OpenMP.
 *
 * \author Quentin Khan, original file: Berenger Bramas <berenger.bramas@inria.fr>
 * \copyright Please read the license.
 *
 * This class runs a threaded FMM algorithm. It just iterates on a tree and call
 * the kernels with good arguments.  The inspector-executor model is used : the
 * class iterates on the tree and builds an array and works in parallel on this
 * array.
 *
 * This algorithm uses the P2P in a thread safe manner, even if the kernel does
 * not initially take care of it. When working on a leaf, a kernel may want to
 * write to the leaf direct neighbours. To avoid any concurrent write, we use 27
 * colours (which is the maximum number of neighbours for a point in a 3D grid).
 * All leaves of a given colour are separated by a least 2 leaves. This means
 * that all threads can work on the same colour at the same time.
 *
 * For example, in 2D, one would have a grid looking like the following, where
 * each number represents a coloured cell. The grid has been cut to work on
 * cells which have a colour value of 4.
 *
 *     0 1 2 | 0 1 2
 *     3 4 5 | 3 4 5
 *     6 7 8 | 6 7 8
 *     ------+------
 *     0 1 2 | 0 1 2
 *     3 4 5 | 3 4 5
 *     6 7 8 | 6 7 8
 *
 * \note Upon destruction, this class does not free pointers given to its constructor.
 */
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FFmmAlgorithmThreadBalanced : public FAbstractAlgorithm, public FAlgorithmTimers{

    /// Shortened tree iterator class.
    using TreeIterator   = typename OctreeClass::Iterator;
    /// Factorisation of the class holding the zone bounds.
    using ZoneBoundClass = std::pair<MortonIndex, int>;

    OctreeClass* const tree;  ///< The octree to work on.
    KernelClass** kernels;    ///< The kernels.

    /// An array of iterators that points to specific cells in the tree. 
    typename OctreeClass::Iterator* iterArray; 
    int leafCount;            ///< The number of leaves.

    static const int SizeColor = 3*3*3; ///< Leaf colors count, see.
    int colorsLeafCount[SizeColor];     ///< Leaf count for each color.

    const int MaxThreads;     ///< The maximum number of threads.
    const int OctreeHeight;   ///< The height of the given tree.

    /// The vector containing the costzones
    const std::vector<std::vector<ZoneBoundClass>>& costzones;
    
public:
    /** 
     * \brief Class constructor
     * 
     * The constructor needs the octree and the kernel used for computation.
     *
     * \warning Internally, one kernel is built for each thread, and each works
     * on its own copy. This means the kernel cannot assume anything about the
     * parts of the tree it will be executed on.
     *
     * \param inTree      The octree to work on.
     * \param inKernels   The kernel to call.
     * \param inCostzones The cost zones the threads will follow.
     *
     * \except An exception is thrown if one of the arguments is NULL.
     */
    FFmmAlgorithmThreadBalanced(OctreeClass* const inTree,
                                KernelClass* const inKernel,
                                const std::vector<std::vector<ZoneBoundClass>>&
                                inCostzones) :
        tree(inTree) , 
        kernels(nullptr),
        iterArray(nullptr),
        leafCount(0),
        MaxThreads(omp_get_max_threads()),
        OctreeHeight(tree->getHeight()),
        costzones(inCostzones) {
        
        FAssertLF(tree, "Tree cannot be null.");
        FAssertLF(inCostzones.size() <= static_cast<unsigned int>(MaxThreads),
                  std::string("Thread count is inferior to cost zone count.") +
                  std::to_string(MaxThreads) +
                  std::string(" / ") +
                  std::to_string(inCostzones.size()));

        this->kernels = new KernelClass*[MaxThreads];

        #pragma omp parallel for schedule(static)
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread) {
            #pragma omp critical (InitFFmmAlgorithmThreadBalanced)
            {
                this->kernels[idxThread] = new KernelClass(*inKernel);
            }
        }

        FAbstractAlgorithm::setNbLevelsInTree(OctreeHeight);

        FLOG(FLog::Controller << "FFmmAlgorithmThreadBalanced (Max Thread " << MaxThreads << ")\n");
    }

    /** \brief Default destructor */
    virtual ~FFmmAlgorithmThreadBalanced(){
        for(int idxThread = 0 ; idxThread < MaxThreads ; ++idxThread){
            delete this->kernels[idxThread];
        }
        delete [] this->kernels;
    }

protected:
    /**
     * \brief Runs the complete algorithm.
     *
     * \param operationsToProceed A flag combinaison to specifiy the operators to use. See FFmmOperations in FCoreCommon.hpp.
     */
    void executeCore(const unsigned operationsToProceed) override {

        for(int idxColor = 0 ; idxColor < SizeColor ; ++idxColor){
            this->colorsLeafCount[idxColor] = 0;
        }

        // Count leaves and color them.
        leafCount = 0;
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            ++leafCount;
            const FTreeCoordinate& coord =
                octreeIterator.getCurrentCell()->getCoordinate();
            ++this->colorsLeafCount[(coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3)];

        } while(octreeIterator.moveRight());

        iterArray = new typename OctreeClass::Iterator[leafCount];
        FAssertLF(iterArray, "iterArray bad alloc");

        Timers[P2MTimer].tic();
        if(operationsToProceed & FFmmP2M)
            bottomPass();
        Timers[P2MTimer].tac();

        Timers[M2MTimer].tic();
        if(operationsToProceed & FFmmM2M)
            upwardPass();
        Timers[M2MTimer].tac();

        Timers[M2LTimer].tic();
        if(operationsToProceed & FFmmM2L)
            transferPass();
        Timers[M2LTimer].tac();

        Timers[L2LTimer].tic();
        if(operationsToProceed & FFmmL2L)
            downardPass();
        Timers[L2LTimer].tac();

        Timers[NearTimer].tic();
        if( (operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P) )
            directPass((operationsToProceed & FFmmP2P),(operationsToProceed & FFmmL2P));
        Timers[NearTimer].tac();

        delete [] iterArray;
        iterArray = nullptr;
    }

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** \brief Runs the P2M kernel. */
    void bottomPass(){
        FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
        FLOG( FTic counterTime );

        TreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();

        // One pair per zone.
        std::vector< std::pair<TreeIterator, int> > iterVector(0);

        std::cerr << "Defining iterators" << std::endl;
        // Find iterators to leaf portion of each zone.
        for( std::vector<ZoneBoundClass> zone : costzones ) {
            iterVector.push_back(
                std::pair<TreeIterator,int>(
                    octreeIterator,       // Iterator to the current cell
                    zone.back().second)); // Cell count in zone

            // Move iterator to end of zone (which is the first of the next zone)
            for( int idx = 0; idx < zone.back().second; idx++) {
                octreeIterator.moveRight();
            }
        }

        FLOG( FTic computationCounter );

        #pragma omp parallel
        {
            const int threadIdx = omp_get_thread_num();
            KernelClass * const myThreadkernels = kernels[threadIdx];
            TreeIterator zoneIterator = iterVector.at(threadIdx).first;
            int zoneCellCount = iterVector[threadIdx].second;

            // Call P2M on cells
            while ( zoneCellCount-- > 0 ) {
                myThreadkernels->P2M(zoneIterator.getCurrentCell(),     // Cell
                                     zoneIterator.getCurrentListSrc()); // Particles
                zoneIterator.moveRight();
            }
        }

        FLOG( computationCounter.tac() );

        FLOG( FLog::Controller << "\tFinished (@Bottom Pass (P2M) = "
              << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : "
              << computationCounter.elapsed() << " s\n" );

    }
    
    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////
    
    /** \brief Runs the M2M kernel. */
    void upwardPass() {
        FLOG( FLog::Controller.write("\tStart Upward Pass\n").write(FLog::Flush); );
        FLOG( FTic counterTime );
        FLOG( FTic computationCounter );

        // Start from leaf level - 1
        TreeIterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        octreeIterator.moveUp(); // Avoid leaf level

        for(int idxLevel = OctreeHeight - 2 ;
            idxLevel > FAbstractAlgorithm::lowerWorkingLevel-1 ;
            --idxLevel)
        {
            octreeIterator.moveUp();
        }

        // Avoids jumping to the left of the tree when changing levels
        TreeIterator avoidGotoLeftIterator(octreeIterator);


        // Stores the iterators to the beginning of zones *per level* in REVERSE order!
        // ie. level in REVERSE ORDER > zones > (iterator,cell_count)
        std::vector< std::vector< std::pair<TreeIterator, int> > >
            reverseLevelIterVector(0);

        const int startIdx =
            FMath::Min(OctreeHeight - 2,
                       FAbstractAlgorithm::lowerWorkingLevel - 1);
        for( int idxLevel = startIdx ;
             idxLevel >= FAbstractAlgorithm::upperWorkingLevel ;
             --idxLevel )
        {
            std::vector< std::pair<TreeIterator, int> > tempVect;
            // Find iterators to leaf portion of each zone.
            for( std::vector<ZoneBoundClass> zone : costzones ) {
                tempVect.push_back(
                    std::pair<TreeIterator,int>(
                        octreeIterator,          // Iterator to the current cell
                        zone[idxLevel].second)); // Cell count in zone
                // Get iterator to end of zone (which is the first of the next zone)
                for( int idx = 0; idx < zone[idxLevel].second; idx++) {
                    octreeIterator.moveRight();
                }

            }
            reverseLevelIterVector.emplace_back(tempVect);

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();
        }

        // for each level from bottom to top
        for( std::vector< std::pair<TreeIterator, int> > levelIterVector :
                 reverseLevelIterVector ) {

            int idxLevel = -1;
            FLOG(FTic counterTimeLevel);

            FLOG(computationCounter.tic());
            #pragma omp parallel
            {
                const int threadNum = omp_get_thread_num();
                KernelClass * const myThreadkernels = kernels[threadNum];
                TreeIterator zoneIterator = levelIterVector[threadNum].first;
                int zoneCellCount = levelIterVector[threadNum].second;

                while(zoneCellCount-- > 0) {
                    // We need the current cell and the child
                    // child is an array (of 8 child) that may be null
                    myThreadkernels->M2M( zoneIterator.getCurrentCell(),
                                          zoneIterator.getCurrentChild(),
                                          zoneIterator.level());
                    zoneIterator.moveRight();
                }

                idxLevel = zoneIterator.level();
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

    /** \brief Runs the M2L kernel. */
    void transferPass(){

        FLOG( FLog::Controller.write("\tStart Downward Pass (M2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        TreeIterator octreeIterator(tree);
        octreeIterator.moveDown();

        for(int idxLevel = 2 ; idxLevel < FAbstractAlgorithm::upperWorkingLevel ; ++idxLevel){
            octreeIterator.moveDown();
        }

        TreeIterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ;
            idxLevel < FAbstractAlgorithm::lowerWorkingLevel ;
            ++idxLevel )
        {
            FLOG(FTic counterTimeLevel);

            std::vector< std::pair<TreeIterator, int> > iterVector;
            // Find iterators to leaf portion of each zone.
            for( std::vector<ZoneBoundClass> zone : costzones ) {
                iterVector.push_back(
                    std::pair<TreeIterator,int>(
                        octreeIterator,          // Iterator to the current cell
                        zone[idxLevel].second)); // Cell count in zone
                // Get iterator to end of zone (which is the first of the next zone)
                for( int idx = 0; idx < zone[idxLevel].second; idx++) {
                    octreeIterator.moveRight();
                }

            }

            octreeIterator.moveDown();
            octreeIterator.gotoLeft();

            FLOG(computationCounter.tic());

            #pragma omp parallel
            {
                const int threadNum = omp_get_thread_num();
                KernelClass * const myThreadkernels = kernels[threadNum];
                const CellClass* neighbours[343];
                TreeIterator zoneIterator = iterVector[threadNum].first;
                int zoneCellCount = iterVector[threadNum].second;
                
                while(zoneCellCount-- > 0) {
                    const int counter =
                        tree->getInteractionNeighbors(
                            neighbours,
                            zoneIterator.getCurrentGlobalCoordinate(),
                            idxLevel);
                    if(counter)
                        myThreadkernels->M2L(
                            zoneIterator.getCurrentCell(),
                            neighbours,
                            counter,
                            idxLevel);
                    zoneIterator.moveRight();
                }

                myThreadkernels->finishedLevelM2L(idxLevel);

            }

            FLOG( computationCounter.tac() );
            FLOG( FLog::Controller << "\t\t>> Level " 
                                   << idxLevel << " = "  
                                   << counterTimeLevel.tacAndElapsed()
                                   << "s\n" );
        }
                
        FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "
              << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation : "
              << computationCounter.cumulated() << " s\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downward
    /////////////////////////////////////////////////////////////////////////////

    /** \brief Runs the L2L kernel. */
    void downardPass(){

        FLOG( FLog::Controller.write("\tStart Downward Pass (L2L)\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);

        TreeIterator octreeIterator(tree);
        octreeIterator.moveDown();

        for(int idxLevel = 2 ; idxLevel < FAbstractAlgorithm::upperWorkingLevel ; ++idxLevel){
            octreeIterator.moveDown();
        }

        TreeIterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = FAbstractAlgorithm::lowerWorkingLevel - 1;
        // for each levels excepted leaf level
        for(int idxLevel = FAbstractAlgorithm::upperWorkingLevel ; 
            idxLevel < heightMinusOne ;
            ++idxLevel ) 
        {
            FLOG(FTic counterTimeLevel);

            std::vector< std::pair<TreeIterator, int> > iterVector;
            // Find iterators to leaf portion of each zone.
            for( std::vector<ZoneBoundClass> zone : costzones ) {
                iterVector.push_back(
                    std::pair<TreeIterator,int>(
                        octreeIterator,          // Iterator to the current cell
                        zone[idxLevel].second)); // Cell count in zone
                // Get iterator to end of zone (which is the first of the next zone)
                for( int idx = 0; idx < zone[idxLevel].second; idx++) {
                    octreeIterator.moveRight();
                }
            }
            octreeIterator.gotoLeft();
            octreeIterator.moveDown();

            FLOG(computationCounter.tic());

            #pragma omp parallel
            {
                const int threadNum = omp_get_thread_num();
                KernelClass * const myThreadkernels = kernels[threadNum];
                TreeIterator zoneIterator = iterVector[threadNum].first;
                int zoneCellCount = iterVector[threadNum].second;
                
                while( zoneCellCount-- > 0 ) {
                    myThreadkernels->L2L(
                        zoneIterator.getCurrentCell(),
                        zoneIterator.getCurrentChild(),
                        idxLevel);
                    zoneIterator.moveRight();
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

    /**
     * \brief Runs the P2P & L2P kernels.
     *
     * \param p2pEnabled Run the P2P kernel.
     * \param l2pEnabled Run the L2P kernel.
     */
    void directPass(const bool p2pEnabled, const bool l2pEnabled){
        FLOG( FLog::Controller.write("\tStart Direct Pass\n").write(FLog::Flush); );
        FLOG(FTic counterTime);
        FLOG(FTic computationCounter);
        FLOG(FTic computationCounterP2P);

        omp_lock_t lockColor[SizeColor];
        for(int idxColor = 0 ; idxColor < SizeColor ; ++idxColor){
            omp_init_lock(&lockColor[idxColor]);
        }

        struct LeafData{
            MortonIndex index;
            CellClass* cell;
            ContainerClass* targets;
            ContainerClass* sources;
        };
        LeafData* const leafsDataArray = new LeafData[this->leafCount];

        const int LeafIndex = OctreeHeight - 1;

        int startPosAtColor[SizeColor];
        startPosAtColor[0] = 0;
        for(int idxColor = 1 ; idxColor < SizeColor ; ++idxColor){
            startPosAtColor[idxColor] = startPosAtColor[idxColor-1] + this->colorsLeafCount[idxColor-1];
        }

        #pragma omp parallel
        {

            const float step = float(this->leafCount) / float(omp_get_num_threads());
            const int start = int(FMath::Ceil(step * float(omp_get_thread_num())));
            const int tempEnd = int(FMath::Ceil(step * float(omp_get_thread_num()+1)));
            const int end = (tempEnd > this->leafCount ? this->leafCount : tempEnd);

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
                const int colorPosition = (coord.getX()%3)*9 + (coord.getY()%3)*3 + (coord.getZ()%3);

                omp_set_lock(&lockColor[colorPosition]);
                const int positionToWork = startPosAtColor[colorPosition]++;
                omp_unset_lock(&lockColor[colorPosition]);

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

            for( int idxColor = 0 ; idxColor < SizeColor ; ++idxColor ) {
                const int endAtThisColor = this->colorsLeafCount[idxColor] + previous;
                const int chunkSize = FMath::Max(1 , endAtThisColor / (omp_get_num_threads() * omp_get_num_threads()));
                #pragma omp for schedule(dynamic, chunkSize)
                for(int idxLeafs = previous ; idxLeafs < endAtThisColor ; ++idxLeafs){
                    LeafData& currentIter = leafsDataArray[idxLeafs];
                    if(l2pEnabled) {
                        myThreadkernels.L2P(currentIter.cell, currentIter.targets);
                    }
                    if(p2pEnabled){
                        // need the current particles and neighbors particles
                        FLOG(if(!omp_get_thread_num()) computationCounterP2P.tic());
                        const int counter = tree->getLeafsNeighbors(neighbors, currentIter.cell->getCoordinate(), LeafIndex);
                        myThreadkernels.P2P(currentIter.cell->getCoordinate(), currentIter.targets,
                                            currentIter.sources, neighbors, counter);
                        FLOG(if(!omp_get_thread_num()) computationCounterP2P.tac());
                    }
                }

                previous = endAtThisColor;
            }
        }

        FLOG(computationCounter.tac());

        delete [] leafsDataArray;
        for(int idxColor = 0 ; idxColor < SizeColor ; ++idxColor){
            omp_destroy_lock(&lockColor[idxColor]);
        }


        FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
        FLOG( FLog::Controller << "\t\t Computation L2P + P2P : " << computationCounter.cumulated()    << " s\n" );
        FLOG( FLog::Controller << "\t\t Computation P2P :       " << computationCounterP2P.cumulated() << " s\n" );

    }

};


#endif //FFMMALGORITHMTHREADBALANCED_HPP
