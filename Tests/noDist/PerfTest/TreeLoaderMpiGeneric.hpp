// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _TREELOADERMPIGENERIC_HPP_
#define _TREELOADERMPIGENERIC_HPP_

#include "PerfTestUtils.hpp"
#include "Utils/FMpi.hpp"

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "BalanceTree/FCostCell.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Containers/FOctree.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Files/FMpiFmaGenericLoader.hpp"
#include "Files/FMpiTreeBuilder.hpp"


/**
 * \brief Genericted FMA file tree loader.
 *
 * See FTreeLoader documentation.
 */
template <typename _FReal, class _BaseClass>
class TreeLoaderMpiGeneric : public FTreeLoader {
public:
    using FReal = _FReal;

    // Required type definitions.
    using BaseClass          = _BaseClass;
    using CellClass          = FCostCell<BaseClass>;
    using ContainerClass     = FP2PParticleContainerIndexed<FReal>;
    using LeafClass          = FSimpleLeaf<FReal, ContainerClass >;
    using OctreeClass        = FOctree<FReal, CellClass, ContainerClass, LeafClass>;

    /// MPI applcation context.
    FMpi* _mpiContext;
    /// File loader.
    FMpiFmaGenericLoader<FReal> _loader;
    /// Required tree member.
    OctreeClass _tree;

    /* Mock particle structure to balance the tree over the processes. */
    struct TestParticle{
        FSize index;             // Index of the particle in the original file.
        FPoint<FReal> position;  // Spatial position of the particle.
        FReal physicalValue;     // Physical value of the particle.
        /* Returns the particle position. */
        const FPoint<FReal>& getPosition(){
            return position;
        }
    };


    /// Constructs the loader and loads the tree.
    TreeLoaderMpiGeneric(FPerfTestParams& params):
        _mpiContext(params.mpiContext),
        _loader(params.filename, _mpiContext->global()),
        _tree(params.treeHeight,
              params.subTreeHeight,
              _loader.getBoxWidth(),
              _loader.getCenterOfBox()) {
        this->loadTree();
    }

    void loadTree() {
        if( 0 == _mpiContext->global().processId())
            std::cout << "Creating & inserting particles" << std::flush;
        
        time.tic();

        // Temporary array of particles read by this process.
        TestParticle* particles = new TestParticle[_loader.getMyNumberOfParticles()];
        memset(particles, 0, (sizeof(TestParticle) * _loader.getMyNumberOfParticles()));

        // Index (in file) of the first particle that will be read by this process.
        FSize idxStart = _loader.getStart();

        // Read particles from parts.
        for(FSize idxPart = 0 ; idxPart < _loader.getMyNumberOfParticles() ; ++idxPart){
            // Store the index (in the original file) the particle.
            particles[idxPart].index = idxPart + idxStart;
            // Read particle from file
            _loader.fillParticle(&particles[idxPart].position,
                                &particles[idxPart].physicalValue);
        }

        // Final vector of particles
        FVector<TestParticle> finalParticles;
        FLeafBalance balancer;

        // Redistribute particules between processes
        FMpiTreeBuilder< FReal, TestParticle >::
            DistributeArrayToContainer(_mpiContext->global(),
                                       particles,
                                       _loader.getMyNumberOfParticles(),
                                       _tree.getBoxCenter(),
                                       _tree.getBoxWidth(),
                                       _tree.getHeight(),
                                       &finalParticles,
                                       &balancer);
        
        // Free temporary array memory.
        delete[] particles;

        // Insert final particles into tree.
        for(FSize idx = 0 ; idx < finalParticles.getSize(); ++idx){
            _tree.insert(finalParticles[idx].position,
                        finalParticles[idx].index,
                        finalParticles[idx].physicalValue);
        }

        time.tac();
        double elapsedTime = time.elapsed(), minTime, maxTime;

        MPI_Reduce(&elapsedTime,&minTime,1,MPI_DOUBLE,MPI_MIN,0,_mpiContext->global().getComm());
        MPI_Reduce(&elapsedTime,&maxTime,1,MPI_DOUBLE,MPI_MAX,0,_mpiContext->global().getComm());

        if( 0 == _mpiContext->global().processId()) {
            std::cout << " Done  ( min-time:" << minTime
                      << " max-time:" << maxTime
                      << " )"
                      << std::endl;

        }
    }

};

#endif
