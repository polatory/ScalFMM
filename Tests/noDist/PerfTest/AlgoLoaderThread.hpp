// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _ALGOLOADERTHREAD_HPP_
#define _ALGOLOADERTHREAD_HPP_

#include <memory>
#include <sstream>

#include "PerfTestUtils.hpp"

#include "Core/FFmmAlgorithmThread.hpp"

/**
 * \brief Algorithm loader for FFmmAlgorithmThread
 *
 * See FAlgoLoader.
 */
template <class _TreeLoader, template<typename> class _KernelLoader>
class AlgoLoaderThread : public FAlgoLoader<_TreeLoader, _KernelLoader> {
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
    using FMMClass = FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>;
    
    /// The tree loader (FTreeLoader) that was used
    TreeLoader& _treeLoader;

    /// The kernel loader (FKernelLoader) that was used
    KernelLoader& _kernelLoader;

    unsigned int _omp_chunk_size; ///< Chunk size for OpenMP

    /// The #FMMClass algorithm instance
    std::unique_ptr<FMMClass> _algo;

    AlgoLoaderThread(FPerfTestParams& params,
                     TreeLoader& treeLoader,
                     KernelLoader& kernelLoader) :
        _treeLoader(treeLoader),
        _kernelLoader(kernelLoader),
        _omp_chunk_size(params.omp_chunk_size),
        _algo(nullptr) {
        
    }

    void run() {
        _algo = std::unique_ptr<FMMClass>(
            new FMMClass(&(_treeLoader._tree), &(_kernelLoader._kernel)));
        _algo->setChunkSize(_omp_chunk_size);

        _algo->execute();
    }


    virtual std::string getRunInfoString() const {
        std::stringstream sstr;
        sstr << "chunksize:" << _omp_chunk_size << " ";
        return sstr.str();
    }

    double getCumulatedTime(FAlgorithmTimers::FTimers timerName) const {
        return _algo->getCumulatedTime(timerName);
    }


};

#endif
