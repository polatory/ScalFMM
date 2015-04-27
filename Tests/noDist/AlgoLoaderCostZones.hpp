#ifndef _ALGOLOADERCOSTZONES_HPP_
#define _ALGOLOADERCOSTZONES_HPP_

#include "PerfTestUtils.hpp"

#include "Core/FFmmAlgorithm.hpp"

#include "BalanceTree/FFmmAlgorithmThreadBalanced.hpp"
#include "BalanceTree/FCostCell.hpp"
#include "BalanceTree/FCostZones.hpp"

/**
 * \brief Algorithm loader for the CostZones algorithm.
 *
 * \warning : This loader requires that the KernelLoader supply a type definition
 * for a `CostKernelClass`
 */
template <class _TreeLoader, template<typename> class _KernelLoader>
class AlgoLoaderCostZones : public FAlgoLoader<_TreeLoader, _KernelLoader> {
public:
    using TreeLoader     = _TreeLoader;
    using KernelLoader   = _KernelLoader<TreeLoader>;

    using FReal          = typename TreeLoader::FReal;
    using CellClass      = typename TreeLoader::CellClass;
    using ContainerClass = typename TreeLoader::ContainerClass;
    using LeafClass      = typename TreeLoader::LeafClass;
    using OctreeClass    = typename TreeLoader::OctreeClass;
    using KernelClass    = typename KernelLoader::KernelClass;
    using CostKernelClass= typename KernelLoader::CostKernelClass;

    static_assert(std::is_base_of<FCostCellTypeTrait, CellClass>::value,
        "The tree cells must derive from FCostCell.");

    using FMMClass = FFmmAlgorithmThreadBalanced
        <OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>;    
    using CostFmmClass = FFmmAlgorithm
        <OctreeClass, CellClass, ContainerClass, CostKernelClass, LeafClass>;

    TreeLoader& _treeLoader;
    KernelLoader& _kernelLoader;

    /// Builds the loader
    AlgoLoaderCostZones(FPerfTestParams& /*params*/,
                        TreeLoader& treeLoader,
                        KernelLoader& kernelLoader) :
        _treeLoader(treeLoader),
        _kernelLoader(kernelLoader) {
        
    }

    /// Computes the tree cells costs then runs the costzones and FMM algorithms.
    void run() {

        OctreeClass* p_tree = &(_treeLoader._tree);

        // Compute tree cells costs
        CostFmmClass costAlgo(p_tree, &(_kernelLoader._costKernel));

        this->time.tic();
        costAlgo.execute();
        this->time.tac();
        std::cout << "Generating tree cost: " << this->time.elapsed() << "s.\n";

        FCostZones<OctreeClass, CellClass> costzones(p_tree, omp_get_max_threads());

        this->time.tic();
        costzones.run();
        this->time.tac();
        std::cout << "Generating cost zones: " << this->time.elapsed() << "s.\n";


        this->time.tic();
        FMMClass algo(p_tree, &(_kernelLoader._kernel), costzones.getZoneBounds(), costzones.getLeafZoneBounds());
        algo.execute();
        this->time.tac();
    }
};



#endif
