// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================

// ==== CMAKE =====
//
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include "ScalFmmConfig.h"
#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FFmaScanfLoader.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

//
// taylor kernel
#include "../../Src/Kernels/Taylor/FTaylorCell.hpp"
#include "../../Src/Kernels/Taylor/FTaylorKernel_Opt.hpp"
//
#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

/**
 * This program compares two different kernels, eg., the Chebyshev kernel with
 * the SphericalBlas kernel.
 */


// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    // get info from commandline
    const char* const filename       = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);
    const unsigned int NbThreads     = FParameters::getValue(argc, argv, "-t", omp_get_max_threads());
    const int DevP                   = FParameters::getValue(argc, argv, "-p", 7);
#ifdef _OPENMP
    omp_set_num_threads(NbThreads);
    std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
#else
    std::cout << "\n>> Sequential version.\n" << std::
#endif
    // init timer
    FTic time;

    struct TestParticle{
        FPoint position;
        FReal forces[3];
        FReal physicalValue;
        FReal potential;
    };
    // open particle file
    FFmaScanfLoader loader(filename);
    if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

    TestParticle* const particles = new TestParticle[loader.getNumberOfParticles()];
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint position;
        FReal physicalValue = 0.0;
        loader.fillParticle(&position,&physicalValue);
        // get copy
        particles[idxPart].position = position;
        particles[idxPart].physicalValue = physicalValue;
        particles[idxPart].potential = 0.0;
        particles[idxPart].forces[0] = 0.0;
        particles[idxPart].forces[1] = 0.0;
        particles[idxPart].forces[2] = 0.0;
    }
    {
        for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
            for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                FP2P::MutualParticles(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                      particles[idxTarget].position.getZ(),particles[idxTarget].physicalValue,
                                      &particles[idxTarget].forces[0],&particles[idxTarget].forces[1],
                                      &particles[idxTarget].forces[2],&particles[idxTarget].potential,
                                particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                particles[idxOther].position.getZ(),particles[idxOther].physicalValue,
                                &particles[idxOther].forces[0],&particles[idxOther].forces[1],
                                &particles[idxOther].forces[2],&particles[idxOther].potential);
            }
        }
    }
    //
    ////////////////////////////////////////////////////////////////////
    //
    {	// begin Taylor kernel

        // accuracy
        const unsigned int ORDER = 7;

        // typedefs
	typedef FTaylorCell<ORDER,1>                                 CellClass;

        typedef FP2PParticleContainerIndexed                          ContainerClass;
        typedef FSimpleLeaf< ContainerClass >                         LeafClass;
        typedef FOctree< CellClass, ContainerClass , LeafClass >      OctreeClass;
	typedef FTaylorKernel<CellClass,ContainerClass,ORDER,1>       KernelClass;
        typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
	//  typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        // init cell class and oct-tree
        OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

        { // -----------------------------------------------------
            std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                      << " particles ..." << std::endl;
            std::cout << "\tHeight : " << TreeHeight << " \t sub-height : "
                      << SubTreeHeight << std::endl;
            time.tic();

            for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                // put in tree
                tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue);
            }

            time.tac();
            std::cout << "Done  " << "(@Creating and Inserting Particles = "
                      << time.elapsed() << "s)." << std::endl;
        } // -----------------------------------------------------

        // -----------------------------------------------------
        std::cout << "\nFFmaTaylor FMM ..." << std::endl;
        time.tic();
        KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClass algorithm(&tree, &kernels);
        algorithm.execute();
        time.tac();
        std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;
        { // Check that each particle has been summed with all other

            tree.forEachLeaf([&](LeafClass* leaf){
                const FReal*const potentials = leaf->getTargets()->getPotentials();
                const FReal*const forcesX = leaf->getTargets()->getForcesX();
                const FReal*const forcesY = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
                const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                const FVector<int>& indexes = leaf->getTargets()->getIndexes();
                for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const int indexPartOrig = indexes[idxPart];
                    potentialDiff.add(particles[indexPartOrig].potential,potentials[idxPart]);
                    fx.add(particles[indexPartOrig].forces[0],forcesX[idxPart]);
                    fy.add(particles[indexPartOrig].forces[1],forcesY[idxPart]);
                    fz.add(particles[indexPartOrig].forces[2],forcesZ[idxPart]);
               }
            });
        }

        // Print for information
        std::cout << "Potential " << potentialDiff << std::endl;
        std::cout << "Fx " << fx << std::endl;
        std::cout << "Fy " << fy << std::endl;
        std::cout << "Fz " << fz << std::endl;
    } // end FFTaylor kernel
    delete[] particles;

    return 0;
}
