// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _ALGOLOADERSECTIONTASK_HPP_
#define _ALGOLOADERSECTIONTASK_HPP_

#include <memory>

#include "PerfTestUtils.hpp"

#include "Core/FFmmAlgorithmSectionTask.hpp"


template <class _TreeLoader, template<typename> class _KernelLoader>
class AlgoLoaderSectionTask : public FAlgoLoader<_TreeLoader, _KernelLoader> {
public:
    using TreeLoader     = _TreeLoader;
    using KernelLoader   = _KernelLoader<TreeLoader>;


    using FReal          = typename TreeLoader::FReal;
    using CellClass      = typename TreeLoader::CellClass;
    using ContainerClass = typename TreeLoader::ContainerClass;
    using LeafClass      = typename TreeLoader::LeafClass;
    using OctreeClass    = typename TreeLoader::OctreeClass;
    using KernelClass    = typename KernelLoader::KernelClass;

    using FMMClass = FFmmAlgorithmSectionTask<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>;
    
    TreeLoader& _treeLoader;
    KernelLoader& _kernelLoader;

    std::unique_ptr<FMMClass> _algo;

    AlgoLoaderSectionTask(FPerfTestParams& /*params*/,
                   TreeLoader& treeLoader,
                   KernelLoader& kernelLoader) :
        _treeLoader(treeLoader),
        _kernelLoader(kernelLoader),
        _algo(nullptr) {
        
    }


    void run() {
        _algo = std::unique_ptr<FMMClass>(
            new FMMClass(&(_treeLoader._tree), &(_kernelLoader._kernel)));
        _algo->execute();
    }

    double getCumulatedTime(FAlgorithmTimers::FTimers timerName) const {
        return _algo->getCumulatedTime(timerName);
    }
};



#endif
