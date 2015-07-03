// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _TREELOADERBASIC_HPP_
#define _TREELOADERBASIC_HPP_

#include "PerfTestUtils.hpp"

#include "Containers/FOctree.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "BalanceTree/FCostCell.hpp"

/**
 * \brief Basic tree loader.
 *
 * See FTreeLoader documentation.
 */
template <typename _FReal, typename _BaseClass>
class TreeLoaderBasic : public FTreeLoader {
public:
    using FReal     = _FReal;
    using BaseClass = _BaseClass;

    // Required type definitions.
    using CellClass          = FCostCell<BaseClass>;
    using ContainerClass     = FP2PParticleContainerIndexed<FReal>;
    using LeafClass          = FSimpleLeaf<FReal, ContainerClass >;
    using OctreeClass        = FOctree<FReal, CellClass, ContainerClass, LeafClass>;

    /// File loader.
    FFmaGenericLoader<FReal> _loader;
    /// Required tree member.
    OctreeClass _tree;

    /// Constructs the loader and loads the tree.
    TreeLoaderBasic(FPerfTestParams& params):
        _loader(params.filename),
        _tree(params.treeHeight,
              params.subTreeHeight,
              _loader.getBoxWidth(),
              _loader.getCenterOfBox()) {
        this->loadTree();
    }

    virtual void loadTree() {
        std::cout << "Creating & inserting particles" << std::flush;
        
        time.tic();
        
        FPoint<FReal> position;
        FReal physicalValue = 0.0;
        for(FSize idxPart = 0 ; idxPart < _loader.getNumberOfParticles() ; ++idxPart) {
            // Read particle per particle from file
            _loader.fillParticle(&position,&physicalValue);
            // put particle in octree
            _tree.insert(position, idxPart, physicalValue);
        }

        time.tac();
        std::cout << " Done  (" << time.elapsed() << " s)." << std::endl;
    }

};

#endif
