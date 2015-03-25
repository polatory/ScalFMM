// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include  "ScalFmmConfig.h"
#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FIOVtk.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithmPeriodic.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"


#include "../../Src/Files/FDlpolyLoader.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#ifdef ScalFMM_USE_BLAS
// chebyshev kernel
#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymKernel.hpp"
#endif

/** Ewal particle is used in the gadget program
 * here we try to make the same simulation
 */

struct EwalParticle {
	FPoint<FReal> position;
	FReal forces[3];
	FReal physicalValue;
	FReal potential;
	int index;
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){

	typedef FP2PParticleContainerIndexed<FReal>                    ContainerClass;
	typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;

#ifdef  ScalFMM_USE_BLAS
	// begin Chebyshev kernel
	// accuracy
	const unsigned int ORDER = 13;
	// typedefs
	typedef FInterpMatrixKernelR<FReal>                                              MatrixKernelClass;
	typedef FChebCell<FReal,ORDER>                                                  CellClass;
	typedef FOctree<FReal, CellClass,ContainerClass,LeafClass>                       OctreeClass;
	typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER>  KernelClass;

#else
	typedef FSphericalCell<FReal>                                    CellClass;
	typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
	typedef FSphericalKernel< FReal, CellClass, ContainerClass >     KernelClass;
	const int DevP          = FParameters::getValue(argc,argv,"-P", 9);
#endif

	typedef FFmmAlgorithmPeriodic<FReal,OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
	typedef FFmmAlgorithm<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass >         FmmClassNoPer;
	///////////////////////What we do/////////////////////////////
	if( FParameters::existParameter(argc, argv, "-help")){
		std::cout << ">> This executable has to be used to compute direct interaction either for periodic or non periodic system.\n";
		std::cout << ">> options are -h H -sh SH  [-per perdeep,  -noper] -fin filenameIN (-bin)  -[no]scale \n";
		std::cout << ">> Recommended files : ../Data/EwalTest_Periodic.run ../Data/EwalTest_NoPeriodic.run\n";
		std::cout << " Options " << std::endl;
		std::cout << "     -per perDeep    " << std::endl;
		std::cout << "     -noper no periodic boundary conditions   " << std::endl;
		std::cout << "     -verbose : print index x y z fx fy fy Q and V" << std::endl;
		std::cout << "     -noscale no scaling and we don't remove the dipole term " << std::endl;
		exit(-1);

	}	//////////////////////////////////////////////////////////////

	const int NbLevels         = FParameters::getValue(argc,argv,"-h",   4);
	const int SizeSubLevels    = FParameters::getValue(argc,argv,"-sh",  2);
	const int PeriodicDeep     = FParameters::getValue(argc,argv,"-per", 3);
	const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/EwalTest_Periodic.run");
	//  file for -saveError option
	std::ofstream errorfile("outputEwaldError.txt",std::ios::out);

	FTic counter;
	//    c     conversion factor for coulombic terms in internal units
	//    c     i.e. (unit(charge)**2/(4 pi eps0 unit(length))/unit(energy)
	const FReal r4pie0 = FReal(138935.4835);
	FReal scaleEnergy, scaleForce ;
	// -----------------------------------------------------
	//  LOADER
	//  -----------------------------------------------------
	std::cout << "Opening : " << filename << "\n";
	FDlpolyLoader  *loader = nullptr ;
	if(FParameters::existParameter(argc, argv, "-bin")){
		loader  = new FDlpolyBinLoader(filename);
	}
	else {
		loader  = new FDlpolyAsciiLoader(filename);
	}

	if(! loader->isOpen()){
		std::cout << "Loader Error, " << filename << " is missing\n";
		return 1;
	}
	// End loading file
	// ---------------------------------------------------
	//        DL_POLY CONSTANT
	//  ---------------------------------------------------
//	bool scale = true ;
	if(FParameters::existParameter(argc, argv, "-noscale")){
//		scale = false ;
		scaleEnergy =  1.0;   // kcal mol^{-1}
		scaleForce   = 1.0 ;           // 10 J mol^{-1} A^{-1}
//		scalePot      = 1.0 ;           //
	}
	else {
		scaleEnergy =  r4pie0 / 418.4 ;   // kcal mol^{-1}
		scaleForce  = -r4pie0 ;           // 10 J mol^{-1} A^{-1}
//		scalePot      = 1.0 ;           //  TO BE DEFINED

	}


	//
#ifndef  ScalFMM_USE_BLAS
	CellClass::Init(DevP);
	std::cout << " $$$$$$$$$$$$$$$  SPHERICAL VERSION $$$$$$$$$$$$"<<std::endl;
	std::cout << " $$$$$$$$$$$$$$$  Order "<< DevP <<" $$$$$$$$$$$$"<<std::endl;
#else
	 const MatrixKernelClass MatrixKernel;

	std::cout << " $$$$$$$$$$$$$$$  CHEBYCHEV VERSION $$$$$$$$$$$$" <<std::endl;
	std::cout << " $$$$$$$$$$$$$$$  Order "<<ORDER <<" $$$$$$$$$$$$"<<std::endl;
#endif
	OctreeClass tree(NbLevels, SizeSubLevels, loader->getBoxWidth(), loader->getCenterOfBox());
	// ---------------------------------------------------------------------------------
	// Insert particles in the Octree
	// ---------------------------------------------------------------------------------   std::cout << "Creating & Inserting " << loader->getNumberOfParticles() << " particles ..." << std::endl;
	std::cout << "\tWidth : " << loader->getBoxWidth() << " \t center x : " << loader->getCenterOfBox().getX()
	    						<< " y : " << loader->getCenterOfBox().getY() << " z : " << loader->getCenterOfBox().getZ() << std::endl;
	std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

	counter.tic();
	FPoint<FReal> electricMoment(0.0,0.0,0.0) ;
	EwalParticle * const particles = new EwalParticle[loader->getNumberOfParticles()];
	memset(particles, 0, sizeof(EwalParticle) * loader->getNumberOfParticles());
	double totalCharge = 0.0;
	for(int idxPart = 0 ; idxPart < loader->getNumberOfParticles() ; ++idxPart){
		//
		loader->fillParticle(&particles[idxPart].position, particles[idxPart].forces,
				&particles[idxPart].physicalValue,&particles[idxPart].index);
		//
		totalCharge += particles[idxPart].physicalValue ;
		electricMoment.incX(particles[idxPart].physicalValue*particles[idxPart].position.getX() );
		electricMoment.incY(particles[idxPart].physicalValue*particles[idxPart].position.getY() );
		electricMoment.incZ(particles[idxPart].physicalValue*particles[idxPart].position.getZ() );
		// reset forces and insert in the tree
		tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue);
	}

	counter.tac();
	double dipoleNorm = electricMoment.norm2() ;
	double volume     =  loader->getBoxWidth()*loader->getBoxWidth()*loader->getBoxWidth() ;
	double coeffCorrection = 2.0*FMath::FPi<FReal>()/volume/3.0 ;

	std::cout << std::endl;
	std::cout << "Total Charge                 = "<< totalCharge <<std::endl;
	std::cout << "Electric Moment            = "<< electricMoment <<std::endl;
	std::cout << "Electric Moment norm2 = "<< dipoleNorm  <<std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << " s)." << std::endl;
	//
	// ---------------------------------------------------------------------------------
	//
	std::cout << "Create kernel & run simu ..." << std::endl;
	counter.tic();

	FTreeCoordinate min, max;

	if( FParameters::existParameter(argc, argv, "-noper") ){
#ifndef  ScalFMM_USE_BLAS

		KernelClass kernels( DevP, NbLevels, loader->getBoxWidth(), loader->getCenterOfBox());
#else

		KernelClass kernels( NbLevels, loader->getBoxWidth(), loader->getCenterOfBox(),&MatrixKernel);
#endif
		FmmClassNoPer algo(&tree,&kernels);
		algo.execute();
	}
	else{
		FmmClass algo(&tree,PeriodicDeep);
		std::cout << "The simulation box is repeated " << algo.theoricalRepetition() << " in X/Y/Z" << std::endl;
		std::cout << "Simulated box: " << algo.extendedBoxWidth()<<std::endl;

#ifndef  ScalFMM_USE_BLAS
		KernelClass kernels( DevP, algo.extendedTreeHeight(), algo.extendedBoxWidth(),algo.extendedBoxCenter());
#else
		KernelClass kernels(algo.extendedTreeHeight(), algo.extendedBoxWidth(),algo.extendedBoxCenter(),&MatrixKernel);
#endif
		algo.setKernel(&kernels);
		algo.execute();
		algo.repetitionsIntervals(&min, &max);
	}

	counter.tac();

	std::cout << "Done  " << "(@FMM Algorithm = " << counter.elapsed() << " s)." << std::endl;
	std::cout << "\n"<< "END FMM "
			<< "-------------------------------------------------------------------------" << std::endl << std::endl ;

	// ----------------------------------------------------------------------------------------------------------
	//                                  Ewald   -- FMM  Comparisons
	// ----------------------------------------------------------------------------------------------------------
	{
		FReal energy = 0.0, tmp;
		FMath::FAccurater<FReal> fx, fy, fz, fmmdfx, fmmdfy, fmmdfz,fmmpot;

		tree.forEachLeaf([&](LeafClass* leaf){
			const FReal*const positionsX = leaf->getTargets()->getPositions()[0];
			const FReal*const positionsY = leaf->getTargets()->getPositions()[1];
			const FReal*const positionsZ = leaf->getTargets()->getPositions()[2];
			//
			 FReal*const potentials  = leaf->getTargets()->getPotentials();
			//
			FReal*const          forcesX       = leaf->getTargets()->getForcesX();
			FReal*const          forcesY       = leaf->getTargets()->getForcesY();
			FReal*const          forcesZ       = leaf->getTargets()->getForcesZ();

			const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();

			const int nbParticlesInLeaf  = leaf->getTargets()->getNbParticles();
			const FVector<int>& indexes  = leaf->getTargets()->getIndexes();
			//
			for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
				const int indexPartOrig = indexes[idxPart];
				//				if(scale){
				// -----------------------------------------------------------------------------------------------------
				// remove polarization correction for DL_POLY
				//
				tmp = positionsX[idxPart]*electricMoment.getX()  + positionsY[idxPart]*electricMoment.getY() + positionsZ[idxPart]*electricMoment.getZ()  ;
				//
				potentials[idxPart] -= 2.0*tmp*coeffCorrection;
				forcesX[idxPart]     -= 2.0*physicalValues[idxPart]*coeffCorrection*electricMoment.getX() ;
				forcesY[idxPart]     -= 2.0*physicalValues[idxPart]*coeffCorrection*electricMoment.getY() ;
				forcesZ[idxPart]     -= 2.0*physicalValues[idxPart]*coeffCorrection*electricMoment.getZ() ;
				//
				// -----------------------------------------------------------------------------------------------------
				// Rescale FMM forces to compare with DL_POLY forces
				//
				forcesX[idxPart] *= scaleForce;
				forcesY[idxPart] *= scaleForce;
				forcesZ[idxPart] *= scaleForce;
				//				}

				fx.add(particles[indexPartOrig].forces[0],forcesX[idxPart]); // Ewald - FMM
				fy.add(particles[indexPartOrig].forces[1],forcesY[idxPart]);
				fz.add(particles[indexPartOrig].forces[2],forcesZ[idxPart]);

				energy +=  potentials[idxPart]*physicalValues[idxPart];
				//
				if(FParameters::existParameter(argc, argv, "-verbose")){
					std::cout << ">> index " << particles[indexPartOrig].index << std::endl;
					std::cout << "Good x " << particles[indexPartOrig].position.getX() << " y " << particles[indexPartOrig].position.getY() << " z " << particles[indexPartOrig].position.getZ() << std::endl;
					std::cout << std::fixed  << std::setprecision(5) ;

					std::cout << "FMM  x " << positionsX[idxPart] << " y " << positionsY[idxPart] << " z " << positionsZ[idxPart] << std::endl;
					std::cout << "Good fx " <<particles[indexPartOrig].forces[0] << " fy " << particles[indexPartOrig].forces[1] << " fz " << particles[indexPartOrig].forces[2] << std::endl;
					std::cout << "FMM  fx " << forcesX[idxPart] << " fy " <<forcesY[idxPart] << " fz " << forcesZ[idxPart] << std::endl;
					std::cout << std::scientific  << std::setprecision(5) ;
					std::cout << "Diff  fx " << particles[indexPartOrig].forces[0]-forcesX[idxPart] << " fy " <<particles[indexPartOrig].forces[1]-forcesY[idxPart] << " fz " << particles[indexPartOrig].forces[2]-forcesZ[idxPart] << std::endl;
					//					std::cout << "GOOD physical value " << particles[indexPartOrig].physicalValue << " potential " << particles[indexPartOrig].potential << std::endl;
					//					std::cout << "FMM  physical value " << physicalValues[idxPart] << " potential " << potentials[idxPart] <<  " energy cumul " << energy<<std::endl;
					std::cout << std::fixed  << std::setprecision(5) ;
					std::cout << "\n";
				}
			}
		});   // ForEachLeaf
		energy *= 0.5*scaleEnergy ;
		//
		// Check results

		printf("FmmEwald Difference between FMM and Ewald DL_POLY:\n");
		printf("FmmEwald Fx diff is = \n");
		printf("FmmEwald   L2Norm  %e\n",fx.getRelativeL2Norm());
		printf("FmmEwald   InfNorm %e\n",fx.getRelativeInfNorm());
		printf("FmmEwald Fy diff is = \n");
		printf("FmmEwald   L2Norm  %e\n",fy.getRelativeL2Norm());
		printf("FmmEwald   InfNorm %e\n",fy.getInfNorm());
		printf("FmmEwald Fz diff is = \n");
		printf("FmmEwald   L2Norm  %e\n",fz.getRelativeL2Norm());
		printf("FmmEwald   InfNorm %e\n",fz.getRelativeInfNorm());
		//
		double L2error = (fx.getL2Norm()*fx.getL2Norm() + fy.getL2Norm()*fy.getL2Norm()  + fz.getL2Norm() *fz.getL2Norm()  );
		printf("FmmEwald RMS Force Error= %e\n",FMath::Sqrt(L2error/static_cast<double>(loader->getNumberOfParticles()))) ;
		//
		printf("FmmEwald Energy FMM=   %.12e\n",energy);
		printf("FmmEwald Energy EWALD= %.12e\n",loader->getEnergy());
		printf("FmmEwald Energy EWALD - Energy FMM = %e\n",loader->getEnergy()-energy);
		printf("FmmEwald |Energy EWALD -Energy FMM|/Energy EWALD= %e\n",FMath::Abs((loader->getEnergy()-energy)/loader->getEnergy()));

	}
	//
	// ----------------------------------------------------------------------------
	//
	delete[] particles;

	return 0;
}



