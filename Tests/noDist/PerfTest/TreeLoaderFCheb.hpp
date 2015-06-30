// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _TREELOADERFCHEB_HPP_
#define _TREELOADERFCHEB_HPP_

#include "Kernels/Chebyshev/FChebCell.hpp"

#include "TreeLoaderBasic.hpp"


/**
 * \brief Tree loader for a Chebyshev cell type tree.
 *
 * See FTreeLoader and TreeLoaderBasic documentation.
 */
template <typename _FReal, int _ORDER>
class TreeLoaderFCheb : public TreeLoaderBasic<_FReal, FChebCell<_FReal, _ORDER> > {
public:

    enum {ORDER=_ORDER};

    /// Constructs the loader and loads the tree.
    TreeLoaderFCheb(FPerfTestParams& params):
        TreeLoaderBasic<_FReal, FChebCell<_FReal, _ORDER>>(params)
    {}

};

#endif
