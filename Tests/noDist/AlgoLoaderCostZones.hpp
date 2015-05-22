// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _ALGOLOADERCOSTZONES_HPP_
#define _ALGOLOADERCOSTZONES_HPP_

#include <memory>
#include <sstream>

#include "PerfTestUtils.hpp"

#include "Core/FFmmAlgorithmThread.hpp"

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
    using CostFmmClass = FFmmAlgorithmThread
        <OctreeClass, CellClass, ContainerClass, CostKernelClass, LeafClass>;

    std::stringstream _infostring;
    TreeLoader& _treeLoader;
    KernelLoader& _kernelLoader;

    std::unique_ptr<FMMClass> _algo;
    

    /// Builds the loader
    AlgoLoaderCostZones(FPerfTestParams& /*params*/,
                        TreeLoader& treeLoader,
                        KernelLoader& kernelLoader) :
        _treeLoader(treeLoader),
        _kernelLoader(kernelLoader),
        _algo(nullptr) {
        
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
        _infostring << "costgen:" << this->time.elapsed() << " ";

        FCostZones<OctreeClass, CellClass> costzones(p_tree, omp_get_max_threads());

        this->time.tic();
        costzones.run();
        this->time.tac();
        std::cout << "Generating cost zones: " << this->time.elapsed() << "s.\n";
        _infostring << "zonegen:" << this->time.elapsed() << " ";

        this->time.tic();
        _algo = std::unique_ptr<FMMClass>(
            new FMMClass(p_tree, &(_kernelLoader._kernel),
                         costzones.getZoneBounds(), costzones.getLeafZoneBounds()));
        _algo->execute();
        this->time.tac();
    }

    std::string getRunInfoString() const {
        return _infostring.str();
    }
};



#endif
