// ===================================================================================
// Copyright ScalFmm 2016 INRIA, Olivier Coulaud, Bérenger Bramas,
// Matthias Messner olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the
// FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
// An extension to the license is given to allow static linking of scalfmm
// inside a proprietary application (no matter its license).
// See the main license file for more details.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================

// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE



#ifndef FFMMALGORITHMTHREADBALANCED_HPP
#define FFMMALGORITHMTHREADBALANCED_HPP


#include "../Src/Utils/FAssert.hpp"
#include "../Src/Utils/FLog.hpp"

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FGlobal.hpp"
#include "../Src/Utils/FAlgorithmTimers.hpp"

#include "../Src/Containers/FOctree.hpp"

#include "../Src/Core/FCoreCommon.hpp"

#include "../Src/BalanceTree/FCoordColour.hpp"
#include "../Src/BalanceTree/FCostZones.hpp"
#include "../Utils/FEnv.hpp"


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
    using ZoneBoundClass = typename FCostZones<OctreeClass, CellClass>::BoundClass;

    OctreeClass* const tree;  ///< The octree to work on.
    KernelClass** kernels;    ///< The kernels.

    const int MaxThreads;     ///< The maximum number of threads.
    const int OctreeHeight;   ///< The height of the given tree.

    /// The vector containing the internal costzones
    const std::vector<std::vector<ZoneBoundClass>>& costzones;

    /// The vector containing the leaf level costzones
    const std::vector<std::vector<ZoneBoundClass>>& leafcostzones;
    
public:
    /** 
     * \brief Class constructor
     * 
     * The constructor needs the octree, the kernel used for computation and the
     * cost zones.
     *
     * \warning Internally, one kernel is built for each thread, and each works
     * on its own copy. This means the kernel cannot assume anything about the
     * parts of the tree it will be executed on.
     *
     * \param inTree      The octree to work on.
     * \param inKernels   The kernel to call.
     * \param internalCostZones The far-field cost zones for each thread.
     * \param leafCostZones The near-field cost zones for each thread.
     *
     * \except An exception is thrown if one of the arguments is NULL.
     * \except An assertion checks that the number of threads is the same as the
     * number of zones.
     */
    FFmmAlgorithmThreadBalanced(
        OctreeClass* const inTree,
        KernelClass* const inKernel,
        const std::vector<std::vector<ZoneBoundClass>>& internalCostzones,
        const std::vector<std::vector<ZoneBoundClass>>& leafCostzones) :
        tree(inTree) , 
        kernels(nullptr),
        MaxThreads(FEnv::GetValue("SCALFMM_ALGO_NUM_THREADS",omp_get_max_threads())),
        OctreeHeight(tree->getHeight()),
        costzones(internalCostzones),
        leafcostzones(leafCostzones) {
        
        FAssertLF(tree, "Tree cannot be null.");
        FAssertLF(internalCostzones.size() == static_cast<unsigned int>(MaxThreads),
                  std::string("Thread count is different from cost zone count (") +
                  std::to_string(MaxThreads) + std::string(" : ") +
                  std::to_string(internalCostzones.size()) + ")." );

        // Allocation of one kernel per thread (in case of impossible concurent use)
        this->kernels = new KernelClass*[MaxThreads];
        // Allocation is done so that kernels are closest to their core.
        #pragma omp parallel num_threads(MaxThreads)
        {
            #pragma omp critical (InitFFmmAlgorithmThreadBalanced)
            {
                this->kernels[omp_get_thread_num()] = new KernelClass(*inKernel);
            }
        }

        FAbstractAlgorithm::setNbLevelsInTree(OctreeHeight);

        FLOG(std::cerr << "FFmmAlgorithmThreadBalanced (Max Thread " << MaxThreads << ")\n");
    }

    /** \brief Destructor */
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
     * \param operationsToProceed A flag combinaison to specify the operators
     * to use. See FFmmOperations in FCoreCommon.hpp.
     */
    void executeCore(const unsigned operationsToProceed) override {

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

    }

    /////////////////////////////////////////////////////////////////////////////
    // P2M
    /////////////////////////////////////////////////////////////////////////////

    /** \brief Runs the P2M kernel.
     *
     * We retrieve the far-field cost zones for the tree leaf level. Each thread
     * runs the P2M kernel on the cells within its zone.
     *
     * \note The lower level of the far-field zones is used.
     */
    void bottomPass(){
        FLOG( std::cerr << "\tStart Bottom Pass" << std::endl );
        FLOG( FTic timer );

        #pragma omp parallel num_threads(MaxThreads)
        {
            // Thread index ( = zone index )
            const int threadIdx = omp_get_thread_num();
            // Thread kernels
            KernelClass * const threadKernels = kernels[threadIdx];
            // Iterator to the zone's first cell
            TreeIterator zoneIterator = costzones.at(threadIdx).at(OctreeHeight-1).first;
            // Cell count in zone
            int zoneCellCount         = costzones.at(threadIdx).at(OctreeHeight-1).second;
            
            // Call P2M on cells
            while ( zoneCellCount-- > 0 ) {
                threadKernels->P2M(zoneIterator.getCurrentCell(),     // Cell
                                   zoneIterator.getCurrentListSrc()); // Particles
                zoneIterator.moveRight();
            }
        } // sync barrier

        FLOG( timer.tac() );
        FLOG( std::cerr << "\tFinished (@Bottom Pass (P2M) = "
                        << timer.elapsed() << " s)" << std::endl );
    }
    
    /////////////////////////////////////////////////////////////////////////////
    // Upward
    /////////////////////////////////////////////////////////////////////////////
    
    /** \brief Runs the M2M kernel.
     *
     * The M2M kernel is run on every node of the tree starting from the
     * lowest level up to the topmost. There is a synchronisation barrier
     * between each level.
     *
     */
    void upwardPass() {
        FLOG( std::cerr << "\tStart Upward Pass" << std::endl );
        FLOG( FTic timer );

        // Running the kernel from the lowest up to the topmost
        for(int level = FAbstractAlgorithm::lowerWorkingLevel-2;
            level >= FAbstractAlgorithm::upperWorkingLevel;
            level--) {

            FLOG( FTic levelTimer );
            FLOG( std::cerr << "\t\t>> Level " << level << std::flush);

            #pragma omp parallel num_threads(MaxThreads)
            {
                // Thread index ( = zone index)
                const int threadNum = omp_get_thread_num(); 
                // Thread kernels
                KernelClass * const myThreadkernels = kernels[threadNum];
                // Iterator to the zone's first cell (for this level)
                TreeIterator zoneIterator = costzones.at(threadNum).at(level).first;
                // Cell count in zone (for this level)
                int zoneCellCount         = costzones.at(threadNum).at(level).second;
               
                // Call M2M on cells
                while(zoneCellCount-- > 0) {
                    // We need the current cell and the child
                    // child is an array (of 8 child) that may be null
                    myThreadkernels->M2M( zoneIterator.getCurrentCell(),
                                          zoneIterator.getCurrentChild(),
                                          zoneIterator.level());
                    zoneIterator.moveRight();
                }
            } // barrier between levels

            FLOG( levelTimer.tac() );
            FLOG( std::cerr << " = "  << levelTimer.elapsed()
                                   << "s" << std::endl );
        }

        FLOG( std::cerr << "\tFinished (@Upward Pass (M2M) = "
                               << timer.tacAndElapsed() << " s)" << std::endl );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Transfer
    /////////////////////////////////////////////////////////////////////////////

    /** \brief Runs the M2L kernel.
     *
     * The M2L kernel is run on every node of the tree starting from the
     * topmost level up to the lowest. There is a synchronisation barrier
     * between each level.
     */
    void transferPass() {
        FLOG( std::cerr << "\tStart Downward Pass (M2L)" << std::endl );
        FLOG( FTic timer );

        // Running the kernel from the topmost to the lowest
        for( int level = FAbstractAlgorithm::upperWorkingLevel ;
             level < FAbstractAlgorithm::lowerWorkingLevel ;
             ++level )
        {
            FLOG( std::cerr << "\t\t>> Level " << level << std::flush );
            FLOG( FTic levelTimer );

            #pragma omp parallel num_threads(MaxThreads)
            {

                // Thread index ( = zone index)
                const int threadNum = omp_get_thread_num(); 
                // Thread kernels
                KernelClass * const myThreadkernels = kernels[threadNum];
                // Iterator to the zone's first cell (for this level)
                TreeIterator zoneIterator = costzones.at(threadNum).at(level).first;
                // Cell count in zone (for this level)
                int zoneCellCount         = costzones.at(threadNum).at(level).second;
                // Array for cell neighbours
                const CellClass* neighbours[342];
                int neighborPositions[342];
                
                // Call M2L kernel on cells
                while(zoneCellCount-- > 0) {
                    const int counter =
                        tree->getInteractionNeighbors(
                            neighbours, neighborPositions,
                            zoneIterator.getCurrentGlobalCoordinate(),
                            level);
                    if(counter)
                        myThreadkernels->M2L(
                            zoneIterator.getCurrentCell(),
                            neighbours,
                            neighborPositions,
                            counter,
                            level);
                    zoneIterator.moveRight();
                }

                myThreadkernels->finishedLevelM2L(level);
            }

            FLOG( levelTimer.tac() );
            FLOG( std::cerr << " = " << levelTimer.elapsed()
                                   << " s"  << std::endl );
        }
                
        FLOG( std::cerr << "\tFinished (@Downward Pass (M2L) = "
              << timer.tacAndElapsed() << " s)\n" );
    }

    /////////////////////////////////////////////////////////////////////////////
    // Downward
    /////////////////////////////////////////////////////////////////////////////

    /** \brief Runs the L2L kernel. 
     *
     * The L2L kernel is run on every node of the tree starting from the
     * topmost level up to the lowest. There is a synchronisation barrier
     * between each level.
     */
    void downardPass(){
        FLOG( std::cerr << "\tStart Downward Pass (L2L)" << std::endl );
        FLOG( FTic timer );

        // for each levels excepted leaf level
        for(int level = FAbstractAlgorithm::upperWorkingLevel ; 
            level < FAbstractAlgorithm::lowerWorkingLevel - 1;
            ++level ) 
        {
            FLOG( std::cerr << "\t\t>> Level " << level << std::flush);
            FLOG( FTic levelTimer );

            #pragma omp parallel num_threads(MaxThreads)
            {
                // Thread index ( = zone index)
                const int threadNum = omp_get_thread_num(); 
                // Thread kernels
                KernelClass * const myThreadkernels = kernels[threadNum];
                // Iterator to the zone's first cell (for this level)
                TreeIterator zoneIterator = costzones.at(threadNum).at(level).first;
                // Cell count in zone (for this level)
                int zoneCellCount         = costzones.at(threadNum).at(level).second;
                
                // Call L2L kernel on cells
                while( zoneCellCount-- > 0 ) {
                    myThreadkernels->L2L(
                        zoneIterator.getCurrentCell(),
                        zoneIterator.getCurrentChild(),
                        level);
                    zoneIterator.moveRight();
                }
            }

            FLOG(levelTimer.tac());
            FLOG( std::cerr << " = "  << levelTimer.elapsed() << " s\n" );
        }

        FLOG( std::cerr << "\tFinished (@Downward Pass (L2L) = " 
                               << timer.tacAndElapsed() << " s)\n" );
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
        FLOG( std::cerr << "\tStart Direct Pass" << std::endl );
        FLOG( FTic timer );

        /** \brief Data handling structure
         *
         * This structure is used to store data related to a leaf.
         */
        struct LeafData {
            MortonIndex index;       ///< Leaf Morton index
            CellClass* cell;         ///< Pointer to the cell in the tree
            ContainerClass* targets; ///< Targets for L2P & P2P kernel
            ContainerClass* sources; ///< Sources for P2P kernel
        };

        #pragma omp parallel num_threads(MaxThreads)
        {
            // Thread index ( = zone index)
            const int threadNum = omp_get_thread_num(); 
            // Thread kernels
            KernelClass * const threadKernels = kernels[threadNum];
            // Iterator to the zone's first cell (for this level)
            TreeIterator zoneIterator;
            // Cell count in zone (for this level)
            int zoneCellCount;
            // Cell neighbours
            ContainerClass* neighbours[26];
            int neighborPositions[26];
            // Cell data
            LeafData leafdata;
            
            // The cells are coloured in order to make concurent work
            // possible. Cells in a colour can all be computed at the same time.
            for( int colourIdx = 0; colourIdx < FCoordColour::range; colourIdx++ ) {
                
                // Get the iterator of the first cell and cell count for this
                // thread ( = zone)
                zoneIterator  = leafcostzones.at(threadNum).at(colourIdx).first;
                zoneCellCount = leafcostzones.at(threadNum).at(colourIdx).second;

                // Skip through all the leaves, work only on those with the
                // right colour
                while( zoneCellCount > 0) {
                    if( FCoordColour::coord2colour(zoneIterator.getCurrentCell()
                                                   ->getCoordinate())
                        == colourIdx) {

                        // Store leaf data
                        leafdata.index   = zoneIterator.getCurrentGlobalIndex();
                        leafdata.cell    = zoneIterator.getCurrentCell();
                        leafdata.targets = zoneIterator.getCurrentListTargets();
                        leafdata.sources = zoneIterator.getCurrentListSrc();
                        
                        // Call L2P kernel
                        if( l2pEnabled ) {
                            threadKernels->L2P(leafdata.cell, leafdata.targets);
                        }

                        // Call P2P kernel
                        if( p2pEnabled ) {
                            // needs the current particles and neighbours particles
                            const int counter =
                                tree->getLeafsNeighbors(neighbours,neighborPositions,
                                                        leafdata.cell->getCoordinate(),
                                                        OctreeHeight - 1); // Leaf level

                            threadKernels->P2P(leafdata.cell->getCoordinate(),
                                                leafdata.targets,
                                                leafdata.sources,
                                                neighbours,
                                                neighborPositions,
                                                counter);
                        }

                        // Decrease cell count only if a cell of the right colour was encountered.
                        zoneCellCount--;
                    }

                    zoneIterator.moveRight();
                }
                // Barrier between computation on two different colours
                #pragma omp barrier
            }
            
        }

        FLOG( timer.tac() );
        FLOG( std::cerr << "\tFinished (@Direct Pass (L2P + P2P) = "
                               << timer.elapsed() << " s)" << std::endl );
    }

};


#endif //FFMMALGORITHMTHREADBALANCED_HPP
