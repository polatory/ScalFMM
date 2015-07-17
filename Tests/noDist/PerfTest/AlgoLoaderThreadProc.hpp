// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _ALGOLOADERTHREADPROC_HPP_
#define _ALGOLOADERTHREADPROC_HPP_

#include <memory>
#include <sstream>

#include "PerfTestUtils.hpp"

#include "Core/FFmmAlgorithmThreadProc.hpp"
#include "Utils/FMpi.hpp"

/**
 * \brief Algorithm loader for FFmmAlgorithmThread
 *
 * See FAlgoLoader.
 */
template <class _TreeLoader, template<typename> class _KernelLoader>
class AlgoLoaderThreadProc : public FAlgoLoader<_TreeLoader, _KernelLoader> {
public:

    // Type definitions, allows them to be reused by other classes
    using TreeLoader     = _TreeLoader;
    using KernelLoader   = _KernelLoader<TreeLoader>;

    using FReal          = typename TreeLoader::FReal;
    using CellClass      = typename TreeLoader::CellClass;
    using ContainerClass = typename TreeLoader::ContainerClass;
    using LeafClass      = typename TreeLoader::LeafClass;
    using OctreeClass    = typename TreeLoader::OctreeClass;
    using KernelClass    = typename KernelLoader::KernelClass;

    /// FMM algorithm class
    using FMMClass = FFmmAlgorithmThreadProc<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>;
    
    FMpi* _mpiContext;

    /// The tree loader (FTreeLoader) that was used
    TreeLoader& _treeLoader;

    /// The kernel loader (FKernelLoader) that was used
    KernelLoader& _kernelLoader;

    /// The #FMMClass algorithm instance
    std::unique_ptr<FMMClass> _algo;
    
    /// Array of MPI gathered cumulated times
    double timers[FAlgorithmTimers::nbTimers] {0};


    AlgoLoaderThreadProc(FPerfTestParams& params,
                     TreeLoader& treeLoader,
                     KernelLoader& kernelLoader) :
        _mpiContext(params.mpiContext),
        _treeLoader(treeLoader),
        _kernelLoader(kernelLoader),
        _algo(nullptr) {
        
    }
    

    void run() {
        _algo = std::unique_ptr<FMMClass>(
            new FMMClass(_mpiContext->global(), &(_treeLoader._tree), &(_kernelLoader._kernel)));
        _algo->execute();
        
        for( int idxTimer = 0; idxTimer < FAlgorithmTimers::nbTimers; ++idxTimer ) {
            timers[idxTimer] = _algo->getCumulatedTime(FAlgorithmTimers::FTimers(idxTimer));
        }

        if( _mpiContext->global().processId() == 0) {
            MPI_Reduce(MPI_IN_PLACE, timers, FAlgorithmTimers::nbTimers, MPI_DOUBLE, MPI_MAX, 0, _mpiContext->global().getComm());
        } else {
            MPI_Reduce(timers, NULL, FAlgorithmTimers::nbTimers, MPI_DOUBLE, MPI_MAX, 0, _mpiContext->global().getComm());
        }
    }

    double getCumulatedTime(FAlgorithmTimers::FTimers timerName) const {
        return timers[timerName];
    }

};

#endif
