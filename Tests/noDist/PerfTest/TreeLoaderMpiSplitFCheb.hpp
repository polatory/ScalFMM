// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _TREELOADERMPISPLITFCHEB_HPP_
#define _TREELOADERMPISPLITFCHEB_HPP_

#include "Kernels/Chebyshev/FChebCell.hpp"

#include "TreeLoaderMpiSplit.hpp"


/**
 * \brief Tree loader for a Chebyshev cell type tree.
 *
 * See FTreeLoader and TreeLoaderBasic documentation.
 */
template <typename _FReal, int _ORDER>
class TreeLoaderMpiSplitFCheb : public TreeLoaderMpiSplit<_FReal, FChebCell<_FReal, _ORDER> > {
public:

    enum {ORDER=_ORDER};

    /// Constructs the loader and loads the tree.
    TreeLoaderMpiSplitFCheb(FPerfTestParams& params):
        TreeLoaderMpiSplit<_FReal, FChebCell<_FReal, _ORDER>>(params)
    {}

};

#endif
