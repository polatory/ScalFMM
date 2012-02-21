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

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Files/FFmaScanfLoader.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

// chebyshev kernel
#include "../Src/Chebyshev/FChebParticle.hpp"
#include "../Src/Chebyshev/FChebLeaf.hpp"
#include "../Src/Chebyshev/FChebCell.hpp"
#include "../Src/Chebyshev/FChebMatrixKernel.hpp"
#include "../Src/Chebyshev/FChebKernels.hpp"

// spherical kernel
#include "../Src/Components/FSimpleLeaf.hpp"
#include "../Src/Kernels/FSphericalKernel.hpp"
#include "../Src/Kernels/FSphericalBlasKernel.hpp"
#include "../Src/Kernels/FSphericalCell.hpp"
#include "../Src/Kernels/FSphericalParticle.hpp"


/**
 * This program compares two different kernels, eg., the Chebyshev kernel with
 * the FmbBlas kernel.
 */



const FReal computeL2norm(unsigned int N, FReal *const u, FReal *const v)
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



const FReal computeINFnorm(unsigned int N, FReal *const u, FReal *const v)
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
	const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
	const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
	const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);

	// init timer
	FTic time;

//	// potentials
//	FReal* p1; p1 = NULL;
//	FReal* p2; p2 = NULL;

//	FReal* p10; p10 = NULL; unsigned int nt1 = 0;
//	FReal* p20; p20 = NULL; unsigned int nt2 = 0;

//	// number of particles
//	unsigned int NumParticles = 0;
	

	////////////////////////////////////////////////////////////////////
	{	// begin Chebyshef kernel

		// accuracy
		const unsigned int ORDER = 3;
		const FReal epsilon = 1e-3;

		unsigned int nt1 = 0;
		FReal* p1;  p1  = NULL;
		FReal* p10; p10 = NULL;
		FReal* f1;  p1  = NULL;
		FReal* f10; p10 = NULL;
		
		// typedefs
		typedef FChebParticle ParticleClass;
		typedef FVector<FChebParticle> ContainerClass;
		typedef FChebLeaf<ParticleClass,ContainerClass> LeafClass;
		typedef FChebMatrixKernelR MatrixKernelClass;
		typedef FChebCell<ORDER> CellClass;
		typedef FOctree<ParticleClass,CellClass,ContainerClass,LeafClass> OctreeClass;
		typedef FChebKernels<ParticleClass,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
		typedef FFmmAlgorithm<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
		//typedef FFmmAlgorithmThread<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
		
		// open particle file
		FFmaScanfLoader<ParticleClass> loader(filename);
		if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");

//		// loader set number of particles
//		NumParticles = loader.getNumberOfParticles();

		// init oct-tree
		OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
		
		{ // -----------------------------------------------------
			std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
			std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
			time.tic();
			loader.fillTree(tree);
			time.tac();
			std::cout << "Done  " << "(@Creating and Inserting Particles = " << time.elapsed() << "s)." << std::endl;
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
		
		
		
		{ // -----------------------------------------------------
			// begin direct computation of first leaf cell
			OctreeClass::Iterator iLeafs(&tree);
			iLeafs.gotoBottomLeft();
			
			// retrieve targets
			const ContainerClass *const Targets = iLeafs.getCurrentListTargets();
			nt1 = Targets->getSize();
			p10 = new FReal [nt1];
			FBlas::setzero(nt1, p10);
			f10 = new FReal [nt1 * 3];
			FBlas::setzero(nt1 * 3, f10);
			
			std::cout << "\nDirect computation of " << nt1 << " target particles ..." << std::endl;
			const MatrixKernelClass MatrixKernel;
			do {
				const ContainerClass *const Sources = iLeafs.getCurrentListSrc();
				unsigned int counter = 0;
				ContainerClass::ConstBasicIterator iTarget(*Targets);
				while(iTarget.hasNotFinished()) {
					const FReal wt = iTarget.data().getPhysicalValue();
					ContainerClass::ConstBasicIterator iSource(*Sources);
					while(iSource.hasNotFinished()) {
						if (&iTarget.data() != &iSource.data()) {
							const FReal ws = iSource.data().getPhysicalValue();
							const FReal one_over_r = MatrixKernel.evaluate(iTarget.data().getPosition(),
																														 iSource.data().getPosition());
							// potential
							p10[counter] += one_over_r * ws;
							// force
							F3DPosition force(iTarget.data().getPosition() - iSource.data().getPosition());
							force *= ((ws*wt) * (one_over_r*one_over_r*one_over_r));
							f10[counter*3 + 0] += force.getX();
							f10[counter*3 + 1] += force.getY();
							f10[counter*3 + 2] += force.getZ();
						}
						iSource.gotoNext();
					}
					counter++;
					iTarget.gotoNext();
				}
			} while(iLeafs.moveRight());
		} // -----------------------------------------------------
		
//		for (unsigned int n=0; n<nt1; ++n)
//			std::cout << p10[n] << " - " << p1[n] << " = " << p10[n]-p1[n] << std::endl;

		std::cout << "\nPotential error:" << std::endl;
		std::cout << "Relative L2 error  = " << computeL2norm( nt1, p10, p1) << std::endl;
		std::cout << "Relative Lmax error = "  << computeINFnorm(nt1, p10, p1) << std::endl;

		std::cout << "\nForce error:" << std::endl;
		std::cout << "Relative L2 error  = " << computeL2norm( nt1*3, f10, f1) << std::endl;
		std::cout << "Relative Lmax error = "  << computeINFnorm(nt1*3, f10, f1) << std::endl << std::endl;


		if (p10!=NULL) delete [] p10;
		if (p1 !=NULL) delete [] p1;
		if (f10!=NULL) delete [] f10;
		if (f1 !=NULL) delete [] f1;

	} // end Chebyshev kernel


	////////////////////////////////////////////////////////////////////
	{	// begin FFmaBlas kernel

		// accuracy
    const int DevP = 8;

		unsigned int nt2 = 0;
		FReal* p2;  p2  = NULL;
		FReal* p20; p20 = NULL;
		FReal* f2;  f2  = NULL;
		FReal* f20; f20 = NULL;

		// typedefs
    typedef FSphericalParticle             ParticleClass;
    typedef FSphericalCell                 CellClass;
    typedef FVector<ParticleClass>         ContainerClass;
    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    //typedef FSphericalBlasKernel<ParticleClass, CellClass, ContainerClass > KernelClass;
    typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass > KernelClass;
		typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

		// open particle file
		FFmaScanfLoader<ParticleClass> loader(filename);
		if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");
		
    // init cell class and oct-tree
    CellClass::Init(DevP);
		OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

    { // -----------------------------------------------------
			std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
			std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
			time.tic();
			loader.fillTree(tree);
			time.tac();
			std::cout << "Done  " << "(@Creating and Inserting Particles = " << time.elapsed() << "s)." << std::endl;
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

		{ // -----------------------------------------------------
			// read potential from particles and write to array p2
			p2 = new FReal [loader.getNumberOfParticles()];
			f2 = new FReal [loader.getNumberOfParticles() * 3];
			OctreeClass::Iterator iLeafs(&tree);
			iLeafs.gotoBottomLeft();
			unsigned int counter = 0;
			do {
				ContainerClass *const Targets = iLeafs.getCurrentListSrc();
				ContainerClass::BasicIterator iTarget(*Targets);
				while(iTarget.hasNotFinished()) {
					p2[counter] = iTarget.data().getPotential();
					f2[counter*3 + 0] = iTarget.data().getForces().getX();
					f2[counter*3 + 1] = iTarget.data().getForces().getY();
					f2[counter*3 + 2] = iTarget.data().getForces().getZ();
					counter++;
					iTarget.data().setPotential(FReal(0.));
					iTarget.data().setForces(FReal(0.), FReal(0.), FReal(0.));
					iTarget.gotoNext();
				}
			} while(iLeafs.moveRight());
		} // -----------------------------------------------------

		{ // -----------------------------------------------------
			// begin direct computation of first leaf cell
			OctreeClass::Iterator iLeafs(&tree);
			iLeafs.gotoBottomLeft();
			
			// retrieve targets
			const ContainerClass *const Targets = iLeafs.getCurrentListTargets();
			nt2 = Targets->getSize();
			p20 = new FReal [nt2];
			FBlas::setzero(nt2, p20);
			f20 = new FReal [nt2 * 3];
			FBlas::setzero(nt2 * 3, f20);
			
			std::cout << "\nDirect computation of " << nt2 << " target particles ..." << std::endl;
			do {
				const ContainerClass *const Sources = iLeafs.getCurrentListSrc();
				unsigned int counter = 0;
				ContainerClass::ConstBasicIterator iTarget(*Targets);
				while(iTarget.hasNotFinished()) {
					ContainerClass::ConstBasicIterator iSource(*Sources);
					while(iSource.hasNotFinished()) {
						if (&iTarget.data() != &iSource.data())
							kernels.directInteraction(const_cast<ParticleClass *>(&iTarget.data()), iSource.data());
						iSource.gotoNext();
					}
					p20[counter] = iTarget.data().getPotential();
					f20[counter*3 + 0] = iTarget.data().getForces().getX();
					f20[counter*3 + 1] = iTarget.data().getForces().getY();
					f20[counter*3 + 2] = iTarget.data().getForces().getZ();
					counter++;
					iTarget.gotoNext();
				}
			} while(iLeafs.moveRight());
		} // -----------------------------------------------------

//		for (unsigned int n=0; n<nt2; ++n)
//			std::cout << p20[n] << " - " << p2[n] << " = " << p20[n]-p2[n] << std::endl;

		std::cout << "\nPotential error:" << std::endl;
		std::cout << "Relative L2 error   = " << computeL2norm( nt2, p20, p2) << std::endl;
		std::cout << "Relative Lmax error = " << computeINFnorm(nt2, p20, p2) << std::endl;

		std::cout << "\nForces error:" << std::endl;
		std::cout << "Relative L2 error   = " << computeL2norm( nt2*3, f20, f2) << std::endl;
		std::cout << "Relative Lmax error = " << computeINFnorm(nt2*3, f20, f2) << std::endl;

		if (p20!=NULL) delete [] p20;
		if (p2 !=NULL) delete [] p2;
		if (f20!=NULL) delete [] f20;
		if (f2 !=NULL) delete [] f2;

	} // end FFmaBlas kernel




	// analyse error


	//std::cout << "\nRelative L2 error  = " << computeL2norm( NumParticles, p1, p2) << std::endl;
	//std::cout << "Relative Lmax error = "  << computeINFnorm(NumParticles, p1, p2) << "\n" << std::endl;

	
//	// free memory
//	if (p10!=NULL) delete [] p10;
//	if (p20!=NULL) delete [] p20;
//	if (p1!=NULL) delete [] p1;
//	if (p2!=NULL) delete [] p2;
	
	return 0;
}







/*
// Check if particles are strictly within its containing cells
const FReal BoxWidthLeaf = BoxWidth / FReal(FMath::pow(2, TreeHeight-1));
OctreeClass::Iterator octreeIterator(&tree);
octreeIterator.gotoBottomLeft();
do{
const CellClass *const LeafCell = octreeIterator.getCurrentCell();
const F3DPosition& LeafCellCenter = LeafCell -> getPosition();
const ContainerClass *const Particles = octreeIterator.getCurrentListSrc();
ContainerClass::ConstBasicIterator particleIterator(*Particles);
while(particleIterator.hasNotFinished()) {
const F3DPosition distance(LeafCellCenter-particleIterator.data().getPosition());
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
