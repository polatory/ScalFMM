// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _TREELOADERMPISPLIT_HPP_
#define _TREELOADERMPISPLIT_HPP_

#include "PerfTestUtils.hpp"
#include "Utils/FMpi.hpp"

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "BalanceTree/FCostCell.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Containers/FOctree.hpp"

#include "Files/FMpiSplitFmaLoader.hpp"


/**
 * \brief Splitted FMA file tree loader.
 *
 * See FTreeLoader documentation.
 */
template <typename _FReal, class _BaseClass>
class TreeLoaderMpiSplit : public FTreeLoader {
public:
    using FReal = _FReal;

    // Required type definitions.
    using BaseClass          = _BaseClass;
    using CellClass          = FCostCell<BaseClass>;
    using ContainerClass     = FP2PParticleContainerIndexed<FReal>;
    using LeafClass          = FSimpleLeaf<FReal, ContainerClass >;
    using OctreeClass        = FOctree<FReal, CellClass, ContainerClass, LeafClass>;

    /// Mpi application context
    FMpi* _mpiContext;
    /// File loader.
    FMpiSplitFmaLoader<FReal> _loader;
    /// Required tree member.
    OctreeClass _tree;

    /// Constructs the loader and loads the tree.
    TreeLoaderMpiSplit(FPerfTestParams& params):
        _mpiContext(params.mpiContext),
        _loader(params.filename,_mpiContext->global().processId()),
        _tree(params.treeHeight,
              params.subTreeHeight,
              _loader.getBoxWidth(),
              _loader.getCenterOfBox()) {
        if( nullptr == _mpiContext ) {
            std::cerr << "No MPI context available" << std::endl;
            exit(-1);
        }

        this->loadTree();
    }

    void loadTree() {
        std::cout << "Creating & inserting particles" << std::flush;
        
        time.tic();
        
        FPoint<FReal> position;
        FReal physicalValue = 0.0;
        for(FSize idxPart = 0 ; idxPart < _loader.getMyNumberOfParticles() ; ++idxPart) {
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
