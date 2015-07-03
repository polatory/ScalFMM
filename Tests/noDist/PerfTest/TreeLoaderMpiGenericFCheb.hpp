// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _TREELOADERMPIGENERICFCHEB_HPP_
#define _TREELOADERMPIGENERICFCHEB_HPP_

#include "Kernels/Chebyshev/FChebCell.hpp"

#include "TreeLoaderMpiGeneric.hpp"


/**
 * \brief Tree loader for a Chebyshev cell type tree.
 *
 * See FTreeLoader and TreeLoaderBasic documentation.
 */
template <typename _FReal, int _ORDER>
class TreeLoaderMpiGenericFCheb : public TreeLoaderMpiGeneric<_FReal, FChebCell<_FReal, _ORDER> > {
public:

    enum {ORDER=_ORDER};

    /// Constructs the loader and loads the tree.
    TreeLoaderMpiGenericFCheb(FPerfTestParams& params):
        TreeLoaderMpiGeneric<_FReal, FChebCell<_FReal, _ORDER>>(params)
    {}

};

#endif
