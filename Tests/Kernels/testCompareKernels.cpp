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
#include "../../Src/Kernels/Chebyshev/FChebParticle.hpp"
#include "../../Src/Kernels/Chebyshev/FChebLeaf.hpp"
#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Chebyshev/FChebMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymKernel.hpp"

// spherical kernel
#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalBlasKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalBlockBlasKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalRotationKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"


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


template <class OctreeClass, class MatrixKernelClass, class ContainerClass>
class DirectInteractionComputer
{

	const MatrixKernelClass MatrixKernel;

public:
	DirectInteractionComputer()
		: MatrixKernel()
	{}
	
	unsigned int operator()(OctreeClass& tree, const unsigned int NumTargetCells,
													FReal* &p, FReal* &f) const
	{
		// begin direct computation of first leaf cell
		typename OctreeClass::Iterator iLeafs(&tree);
		iLeafs.gotoBottomLeft();

		const ContainerClass* *const TargetSets = new const ContainerClass* [NumTargetCells];

		unsigned int n = 0;
		for (unsigned int t=0; t<NumTargetCells; ++t) {
			TargetSets[t] = iLeafs.getCurrentListTargets();
			n += TargetSets[t]->getSize();
			iLeafs.moveRight();
		}
		p = new FReal [n];
		FBlas::setzero(n, p);
		f = new FReal [n * 3];
		FBlas::setzero(n * 3, f);

		std::cout << "\nDirect computation of " << n << " target particles ..." << std::endl;

		unsigned int start = 0;
		for (unsigned int t=0; t<NumTargetCells; ++t) {
			iLeafs.gotoBottomLeft();

			// retrieve targets
			const ContainerClass *const Targets = TargetSets[t];
			
			do {
				const ContainerClass *const Sources = iLeafs.getCurrentListSrc();
				unsigned int counter = start;
				typename ContainerClass::ConstBasicIterator iTarget(*Targets);
				while(iTarget.hasNotFinished()) {
					const FReal wt = iTarget.data().getPhysicalValue();
					typename ContainerClass::ConstBasicIterator iSource(*Sources);
					while(iSource.hasNotFinished()) {
						if (&iTarget.data() != &iSource.data()) {
							const FReal ws = iSource.data().getPhysicalValue();
							const FReal one_over_r = MatrixKernel.evaluate(iTarget.data().getPosition(),
																														 iSource.data().getPosition());
							// potential
							p[counter] += one_over_r * ws;
							// force
							FPoint force(iSource.data().getPosition() - iTarget.data().getPosition());
							force *= ((ws*wt) * (one_over_r*one_over_r*one_over_r));
							f[counter*3 + 0] += force.getX();
							f[counter*3 + 1] += force.getY();
							f[counter*3 + 2] += force.getZ();
						}
						iSource.gotoNext();
					}
					counter++;
					iTarget.gotoNext();
				}
			} while(iLeafs.moveRight());

			start += Targets->getSize();
		}

		delete [] TargetSets;

		return n;
	}

};



// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    // get info from commandline
    const char* const filename       = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);
		const unsigned int NbThreads     = FParameters::getValue(argc, argv, "-t", 1);
		
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
        typedef FChebParticle ParticleClass;
        typedef FVector<FChebParticle> ContainerClass;
        typedef FChebLeaf<ParticleClass,ContainerClass> LeafClass;
        typedef FChebMatrixKernelR MatrixKernelClass;
        typedef FChebCell<ORDER> CellClass;
        typedef FOctree<ParticleClass,CellClass,ContainerClass,LeafClass> OctreeClass;

				//typedef FChebKernel<ParticleClass,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
				typedef FChebSymKernel<ParticleClass,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
				//typedef FFmmAlgorithm<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
				typedef FFmmAlgorithmThread<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;


        // open particle file
        FFmaScanfLoader<ParticleClass> loader(filename);
        if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

        // init oct-tree
        OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

        { // -----------------------------------------------------
            std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
											<< " particles ..." << std::endl;
            std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
            time.tic();
            loader.fillTree(tree);
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

        { // -----------------------------------------------------
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
        } // -----------------------------------------------------

				// compute direct interaction
				const unsigned int NumTargetCells = 3;
				DirectInteractionComputer<OctreeClass,MatrixKernelClass,ContainerClass> direct;
				nt1 = direct(tree, NumTargetCells, p10, f10);


        //		for (unsigned int n=0; n<nt1; ++n)
        //			std::cout << p10[n] << " - " << p1[n] << " = " << p10[n]-p1[n] << std::endl;

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
        typedef FSphericalParticle             ParticleClass;
        typedef FSphericalCell                 CellClass;
        typedef FVector<ParticleClass>         ContainerClass;
        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
        //typedef FSphericalBlasKernel<ParticleClass, CellClass, ContainerClass > KernelClass;
        //typedef FSphericalBlockBlasKernel<ParticleClass, CellClass, ContainerClass > KernelClass;
        //typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass > KernelClass;
        typedef FSphericalBlockBlasKernel<ParticleClass, CellClass, ContainerClass > KernelClass;
        //typedef FSphericalRotationKernel<ParticleClass, CellClass, ContainerClass > KernelClass;
        //typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
        typedef FFmmAlgorithmThread<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        // open particle file
        FFmaScanfLoader<ParticleClass> loader(filename);
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
            loader.fillTree(tree);
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
        { // -----------------------------------------------------
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
        } // -----------------------------------------------------

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







/*
// Check if particles are strictly within its containing cells
const FReal BoxWidthLeaf = BoxWidth / FReal(FMath::pow(2, TreeHeight-1));
OctreeClass::Iterator octreeIterator(&tree);
octreeIterator.gotoBottomLeft();
do{
const CellClass *const LeafCell = octreeIterator.getCurrentCell();
const FPoint& LeafCellCenter = LeafCell -> getPosition();
const ContainerClass *const Particles = octreeIterator.getCurrentListSrc();
ContainerClass::ConstBasicIterator particleIterator(*Particles);
while(particleIterator.hasNotFinished()) {
const FPoint distance(LeafCellCenter-particleIterator.data().getPosition());
std::cout << "center - particle = " << distance << " < " << BoxWidthLeaf/FReal(2.) << std::endl;
if (std::abs(distance.getX())>BoxWidthLeaf/FReal(2.) ||
std::abs(distance.getY())>BoxWidthLeaf/FReal(2.) ||
std::abs(distance.getZ())>BoxWidthLeaf/FReal(2.)) {
std::cout << "stop" << std::endl;
exit(-1);
}
particleIterator.gotoNext();
}
} while(octreeIterator.moveRight());
*/







        //unsigned int sizeFirstLeaf = 0;
        //FReal* firstLeafPotential = 0;
        //FReal* firstLeafForces = 0;
        //{ // -----------------------------------------------------
        //    // begin direct computation of first leaf cell
        //    OctreeClass::Iterator leavesIterator(&tree);
        //    leavesIterator.gotoBottomLeft();
        //    ContainerClass *const firstLeaf = leavesIterator.getCurrentListTargets();
				//
        //    { // Compute the first leaf with it self
        //        ContainerClass::BasicIterator targetParticles(*firstLeaf);
        //        while(targetParticles.hasNotFinished()) {
        //            ContainerClass::ConstBasicIterator sourceParticles(*firstLeaf);
        //            while(sourceParticles.hasNotFinished()) {
        //                if (&targetParticles.data() != &sourceParticles.data()){
        //                    kernels.directInteraction( &targetParticles.data(), sourceParticles.data());
        //                }
        //                sourceParticles.gotoNext();
        //            }
        //            targetParticles.gotoNext();
        //        }
        //    }
				//
        //    while( leavesIterator.moveRight() ) {
        //        ContainerClass::ConstBasicIterator sourceParticles(*leavesIterator.getCurrentListSrc());
        //        while(sourceParticles.hasNotFinished()) {
        //            ContainerClass::BasicIterator targetParticles(*firstLeaf);
        //            while(targetParticles.hasNotFinished()) {
        //                kernels.directInteraction( &targetParticles.data(), sourceParticles.data());
        //                targetParticles.gotoNext();
        //            }
        //            sourceParticles.gotoNext();
        //        }
        //    }
				//
        //    // Copy data
        //    sizeFirstLeaf = firstLeaf->getSize();
				//
        //    firstLeafPotential = new FReal[sizeFirstLeaf];
        //    FBlas::setzero(sizeFirstLeaf, firstLeafPotential);
				//
        //    firstLeafForces = new FReal[sizeFirstLeaf * 3];
        //    FBlas::setzero(sizeFirstLeaf * 3, firstLeafForces);
				//
        //    unsigned int particlePosition = 0;
        //    ContainerClass::ConstBasicIterator targetParticles(*firstLeaf);
        //    while(targetParticles.hasNotFinished()) {
        //        firstLeafPotential[particlePosition] = targetParticles.data().getPotential();                
        //        firstLeafForces[particlePosition*3 + 0] = targetParticles.data().getForces().getX();
        //        firstLeafForces[particlePosition*3 + 1] = targetParticles.data().getForces().getY();
        //        firstLeafForces[particlePosition*3 + 2] = targetParticles.data().getForces().getZ();
        //        targetParticles.gotoNext();
        //        ++particlePosition;
        //    }
        //} // -----------------------------------------------------

        //std::cout << "\nPotential error:" << std::endl;
        //std::cout << "Relative L2 error   = " << computeL2norm( sizeFirstLeaf, firstLeafPotential, fmmParticlesPotential) << std::endl;
        //std::cout << "Relative Lmax error = " << computeINFnorm(sizeFirstLeaf, firstLeafPotential, fmmParticlesPotential) << std::endl;
				//
        //std::cout << "\nForces error:" << std::endl;
        //std::cout << "Relative L2 error   = " << computeL2norm( sizeFirstLeaf*3, firstLeafForces, fmmParticlesForces) << std::endl;
        //std::cout << "Relative Lmax error = " << computeINFnorm(sizeFirstLeaf*3, firstLeafForces, fmmParticlesForces) << std::endl;
