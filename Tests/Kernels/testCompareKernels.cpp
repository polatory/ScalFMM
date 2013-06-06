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
// @FUSE_BLAS
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FFmaScanfLoader.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

// chebyshev kernel
#include "../../Src/Kernels/Chebyshev/FChebLeaf.hpp"
#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Chebyshev/FChebMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymKernel.hpp"

// spherical kernel
#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalBlasKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalBlockBlasKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalRotationKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

/**
 * This program compares two different kernels, eg., the Chebyshev kernel with
 * the SphericalBlas kernel.
 */



FReal computeL2norm(const unsigned int N, const FReal *const u, const FReal *const v)
{
    FReal      dot = FReal(0.);
    FReal diff_dot = FReal(0.);
    for (unsigned int i=0; i<N; ++i) {
        FReal w = v[i] - u[i];
        diff_dot += w    * w;
        dot      += u[i] * u[i];
    }
    return FMath::Sqrt(diff_dot / dot);
}



FReal computeINFnorm(const unsigned int N, const FReal *const u, const FReal *const v)
{
    FReal      max = FReal(0.);
    FReal diff_max = FReal(0.);
    for (unsigned int n=0; n<N; ++n) {
        if (     max<std::abs(u[n]))           max = std::abs(u[n]);
        if (diff_max<std::abs(u[n]-v[n])) diff_max = std::abs(u[n]-v[n]);
    }
    return diff_max / max;
}


// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    // get info from commandline
    const char* const filename       = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);
    const unsigned int NbThreads     = FParameters::getValue(argc, argv, "-t", omp_get_max_threads());

    omp_set_num_threads(NbThreads);

    std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;

    // init timer
    FTic time;

    // only for direct computation of nt1 target particles
    unsigned int nt1 = 0;
    FReal* p10; p10 = NULL;
    FReal* f10; p10 = NULL;

    ////////////////////////////////////////////////////////////////////
    {	// begin Chebyshef kernel

        // accuracy
        const unsigned int ORDER = 7;
        const FReal epsilon = FReal(1e-7);

        FReal* p1;  p1  = NULL;
        FReal* f1;  p1  = NULL;

        // typedefs
        typedef FP2PParticleContainer ContainerClass;
        typedef FSimpleLeaf<ContainerClass> LeafClass;
        typedef FChebMatrixKernelR MatrixKernelClass;
        typedef FChebCell<ORDER> CellClass;
        typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;

        //typedef FChebKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FChebSymKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        //typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
        typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;


        // open particle file
        FFmaScanfLoader loader(filename);
        if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

        // init oct-tree
        OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

        { // -----------------------------------------------------
            std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                      << " particles ..." << std::endl;
            std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
            time.tic();

            FPoint particlePosition;
            FReal physicalValue = 0.0;
            for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                loader.fillParticle(&particlePosition,&physicalValue);
                tree.insert(particlePosition, physicalValue);
            }

            time.tac();
            std::cout << "Done  " << "(@Creating and Inserting Particles = "
                      << time.elapsed() << "s)." << std::endl;
        } // -----------------------------------------------------

        { // -----------------------------------------------------
            std::cout << "\nChebyshev FMM ... " << std::endl;
            time.tic();
            KernelClass kernels(TreeHeight, loader.getCenterOfBox(), loader.getBoxWidth(), epsilon);
            FmmClass algorithm(&tree, &kernels);
            algorithm.execute();
            time.tac();
            std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
        } // -----------------------------------------------------

        /*{ // -----------------------------------------------------
            // read potential from particles and write to array p1
            p1 = new FReal [loader.getNumberOfParticles()];
            f1 = new FReal [loader.getNumberOfParticles() * 3];
            OctreeClass::Iterator iLeafs(&tree);
            iLeafs.gotoBottomLeft();
            unsigned int counter = 0;
            do {
                ContainerClass *const Targets = iLeafs.getCurrentListSrc();
                ContainerClass::BasicIterator iTarget(*Targets);
                while(iTarget.hasNotFinished()) {
                    p1[counter] = iTarget.data().getPotential();
                    f1[counter*3 + 0] = iTarget.data().getForces().getX();
                    f1[counter*3 + 1] = iTarget.data().getForces().getY();
                    f1[counter*3 + 2] = iTarget.data().getForces().getZ();
                    counter++;
                    iTarget.gotoNext();
                }
            } while(iLeafs.moveRight());
        } // -----------------------------------------------------*/

        // compute direct interaction
        const unsigned int NumTargetCells = 3;
        //DirectInteractionComputer<OctreeClass,MatrixKernelClass,ContainerClass> direct;
        //nt1 = direct(tree, NumTargetCells, p10, f10);

        std::cout << "\nPotential error:" << std::endl;
        std::cout << "Relative L2 error  = " << computeL2norm( nt1, p10, p1) << std::endl;
        std::cout << "Relative Lmax error = "  << computeINFnorm(nt1, p10, p1) << std::endl;

        std::cout << "\nForce error:" << std::endl;
        std::cout << "Relative L2 error  = " << computeL2norm( nt1*3, f10, f1) << std::endl;
        std::cout << "Relative Lmax error = "  << computeINFnorm(nt1*3, f10, f1) << std::endl << std::endl;


        if (p1 !=NULL) delete [] p1;
        if (f1 !=NULL) delete [] f1;

    } // end Chebyshev kernel


    ////////////////////////////////////////////////////////////////////
    {	// begin FFmaBlas kernel

        // accuracy
        const int DevP = FParameters::getValue(argc, argv, "-p", 5);

        // typedefs
        typedef FSphericalCell                 CellClass;
        typedef FP2PParticleContainer         ContainerClass;
        typedef FSimpleLeaf< ContainerClass >                     LeafClass;
        typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
        //typedef FSphericalBlasKernel< CellClass, ContainerClass > KernelClass;
        //typedef FSphericalBlockBlasKernel< CellClass, ContainerClass > KernelClass;
        //typedef FSphericalKernel< CellClass, ContainerClass > KernelClass;
        typedef FSphericalBlockBlasKernel< CellClass, ContainerClass > KernelClass;
        //typedef FSphericalRotationKernel< CellClass, ContainerClass > KernelClass;
        //typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
        typedef FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        // open particle file
        FFmaScanfLoader loader(filename);
        if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

        // init cell class and oct-tree
        CellClass::Init(DevP, true); // only for blas
        CellClass::Init(DevP, false);
        OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

        { // -----------------------------------------------------
            std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                      << " particles ..." << std::endl;
            std::cout << "\tHeight : " << TreeHeight << " \t sub-height : "
                      << SubTreeHeight << std::endl;
            time.tic();

            FPoint particlePosition;
            FReal physicalValue = 0.0;
            for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                loader.fillParticle(&particlePosition,&physicalValue);
                tree.insert(particlePosition, physicalValue);
            }

            time.tac();
            std::cout << "Done  " << "(@Creating and Inserting Particles = "
                      << time.elapsed() << "s)." << std::endl;
        } // -----------------------------------------------------

        // -----------------------------------------------------
        std::cout << "\nFFmaBlas FMM ..." << std::endl;
        time.tic();
        KernelClass kernels(DevP, TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClass algorithm(&tree, &kernels);
        algorithm.execute();
        time.tac();
        std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

        FReal*const fmmParticlesPotential = new FReal[loader.getNumberOfParticles()];
        FReal*const fmmParticlesForces = new FReal [loader.getNumberOfParticles() * 3];
        /*{ // -----------------------------------------------------
            // read potential from particles and write to array p2
            OctreeClass::Iterator leavesIterator(&tree);
            leavesIterator.gotoBottomLeft();
            unsigned int particlePosition = 0;
            do {
                ContainerClass *const Targets = leavesIterator.getCurrentListSrc();
                ContainerClass::BasicIterator particlesIterator(*Targets);
                while(particlesIterator.hasNotFinished()) {
                    // Copy current data into array
                    fmmParticlesPotential[particlePosition]    = particlesIterator.data().getPotential();
                    fmmParticlesForces[particlePosition*3 + 0] = particlesIterator.data().getForces().getX();
                    fmmParticlesForces[particlePosition*3 + 1] = particlesIterator.data().getForces().getY();
                    fmmParticlesForces[particlePosition*3 + 2] = particlesIterator.data().getForces().getZ();
                    ++particlePosition;
                    // Reset particle
                    particlesIterator.data().setPotential(FReal(0.));
                    particlesIterator.data().setForces(FReal(0.), FReal(0.), FReal(0.));
                    particlesIterator.gotoNext();
                }
            } while(leavesIterator.moveRight());
        } // -----------------------------------------------------*/

        std::cout << "\nPotential error:" << std::endl;
        std::cout << "Relative L2 error   = "
                  << computeL2norm( nt1, p10, fmmParticlesPotential) << std::endl;
        std::cout << "Relative Lmax error = "
                  << computeINFnorm(nt1, p10, fmmParticlesPotential) << std::endl;

        std::cout << "\nForces error:" << std::endl;
        std::cout << "Relative L2 error   = "
                  << computeL2norm( nt1*3, f10, fmmParticlesForces) << std::endl;
        std::cout << "Relative Lmax error = "
                  << computeINFnorm(nt1*3, f10, fmmParticlesForces) << std::endl;



        delete [] fmmParticlesPotential;
        delete [] fmmParticlesForces;

    } // end FFmaBlas kernel


    // free memory
    if (p10!=NULL) delete [] p10;
    if (f10!=NULL) delete [] f10;

    return 0;
}
