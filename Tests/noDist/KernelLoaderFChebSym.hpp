#ifndef _KERNELLOADERFCHEBSYM_HPP_
#define _KERNELLOADERFCHEBSYM_HPP_

#include "PerfTestUtils.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "BalanceTree/FChebSymCostKernel.hpp"

/**
 * \brief Kernel loader for the symetric Chebyshev kernel.
 *
 * \warning This loader requires that TreeLoader::CellClass inherits from
 * FChebCell.
 *
 * \note This loader also provides the typedef CostKernelClass and a member
 * _costKernel that cam be used by the AlgoLoaderCostZones.
 */
template <typename _TreeLoader>
class KernelLoaderFChebSym : public FKernelLoader<_TreeLoader> {

public:
    // Required type definitions
    using TreeLoader     = _TreeLoader;
    using FReal          = typename TreeLoader::FReal;
    /// Must derive from FChebCell
    using CellClass      = typename TreeLoader::CellClass;
    using ContainerClass = typename TreeLoader::ContainerClass;
    using OctreeClass    = typename TreeLoader::OctreeClass;

    using MatrixKernelClass = FInterpMatrixKernelR<FReal>;
    using KernelClass       = FChebSymKernel <FReal, CellClass, ContainerClass,
                                              MatrixKernelClass, TreeLoader::ORDER>;
    /// Kernel class used to compute the tree cell costs.
    using CostKernelClass = FChebSymCostKernel<FReal, CellClass, ContainerClass, 
                                               MatrixKernelClass, TreeLoader::ORDER,
                                               OctreeClass>;

    // Meaningfull (?) error message.
    static_assert(std::is_base_of<FChebCell<FReal,TreeLoader::ORDER>,CellClass>::value,
                  "TreeLoader::CellClass must derive from FChebCell");    


    const FReal epsilon = 1e-4;
    const MatrixKernelClass _matrixKernel;
    /// kernel used to compute the tree cells interactions.
    KernelClass _kernel;
    /// Kernel used to compute the tree cells costs.
    CostKernelClass _costKernel;

    /// Builds the loader and loads the kernel.
    KernelLoaderFChebSym(FPerfTestParams& /*params*/, TreeLoader& treeLoader) :
        _matrixKernel(),
        _kernel(treeLoader._tree.getHeight(),
                treeLoader._tree.getBoxWidth(),
                treeLoader._tree.getBoxCenter(),
                &_matrixKernel),
        _costKernel(&(treeLoader._tree), epsilon){

    }


};


#endif
