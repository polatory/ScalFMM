// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _TREELOADERFCHEB_HPP_
#define _TREELOADERFCHEB_HPP_

#include "PerfTestUtils.hpp"

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Containers/FOctree.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "BalanceTree/FCostCell.hpp"

/**
 * \brief Tree loader for a Chebyshev cell type tree.
 *
 * See FTreeLoader documentation.
 */
template <typename _FReal = double>
class TreeLoaderFCheb : public FTreeLoader {
public:
    using FReal = _FReal;
    enum {ORDER = 7}; ///< Chebyshev interpolation order.

    // Required type definitions.
    using CellClass          = FCostCell<FChebCell<FReal, ORDER>>;
    using ContainerClass     = FP2PParticleContainerIndexed<FReal>;
    using LeafClass          = FSimpleLeaf<FReal, ContainerClass >;
    using OctreeClass        = FOctree<FReal, CellClass, ContainerClass, LeafClass>;

    /// File loader.
    FFmaGenericLoader<FReal> _loader;
    /// Required tree member.
    OctreeClass _tree;

    /// Constructs the loader and loads the tree.
    TreeLoaderFCheb(FPerfTestParams& params):
        _loader(params.filename),
        _tree(params.treeHeight,
              params.subTreeHeight,
              _loader.getBoxWidth(),
              _loader.getCenterOfBox()) {
        this->loadTree(_loader, _tree);
    }


};

#endif
