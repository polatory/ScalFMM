#ifndef _ALGOLOADERTHREAD_HPP_
#define _ALGOLOADERTHREAD_HPP_

#include "PerfTestUtils.hpp"

#include "Core/FFmmAlgorithmThread.hpp"


template <class _TreeLoader, template<typename> class _KernelLoader>
class AlgoLoaderThread : public FAlgoLoader<_TreeLoader, _KernelLoader> {
public:
    using TreeLoader     = _TreeLoader;
    using KernelLoader   = _KernelLoader<TreeLoader>;

    using FReal          = typename TreeLoader::FReal;
    using CellClass      = typename TreeLoader::CellClass;
    using ContainerClass = typename TreeLoader::ContainerClass;
    using LeafClass      = typename TreeLoader::LeafClass;
    using OctreeClass    = typename TreeLoader::OctreeClass;
    using KernelClass    = typename KernelLoader::KernelClass;

    using FMMClass = FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>;
    
    TreeLoader& _treeLoader;
    KernelLoader& _kernelLoader;

    bool _omp_static_schedule;

    AlgoLoaderThread(FPerfTestParams& params,
                     TreeLoader& treeLoader,
                     KernelLoader& kernelLoader) :
        _treeLoader(treeLoader),
        _kernelLoader(kernelLoader),
        _omp_static_schedule(params.omp_static_schedule) {
        
    }


    void run() {
        FMMClass algo(&(_treeLoader._tree), &(_kernelLoader._kernel),
                      _omp_static_schedule);
        algo.execute();
    }
};

#endif
