// ===================================================================================
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
	FPoint position;
	FReal forces[3];
	FReal physicalValue;
	FReal potential;
	int index;
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){

	typedef FP2PParticleContainerIndexed                      ContainerClass;
	typedef FSimpleLeaf< ContainerClass >                     LeafClass;

#ifdef  ScalFMM_USE_BLAS
	// begin Chebyshev kernel
	// accuracy
	const unsigned int ORDER = 10;
	// typedefs
	typedef FInterpMatrixKernelR                                              MatrixKernelClass;
	typedef FChebCell<ORDER>                                                  CellClass;
	typedef FOctree<CellClass,ContainerClass,LeafClass>                       OctreeClass;
	typedef FChebSymKernel<CellClass,ContainerClass,MatrixKernelClass,ORDER>  KernelClass;
#else
	typedef FSphericalCell                                    CellClass;
	typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
	typedef FSphericalKernel< CellClass, ContainerClass >     KernelClass;
	const int DevP          = FParameters::getValue(argc,argv,"-P", 9);
#endif

	typedef FFmmAlgorithmPeriodic<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
	typedef FFmmAlgorithm<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass >         FmmClassNoPer;
	///////////////////////What we do/////////////////////////////
	std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
	std::cout << ">> options are -h H -sh SH -P p -per PER -f FILE -noper -verbose -gen -direct\n";
	std::cout << ">> Recommended files : ../Data/EwalTest_Periodic.run ../Data/EwalTest_NoPeriodic.run\n";
	//////////////////////////////////////////////////////////////

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
	FReal scaleEnergy, scaleForce  ;
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
	// ---------------------------------------------------
	//        DL_POLY CONSTANT
	//  ---------------------------------------------------
	scaleEnergy =  r4pie0 / 418.4 ;   // kcal mol^{-1}
	scaleForce  = -r4pie0 ;           // 10 J mol^{-1} A^{-1}

	//
#ifndef  ScalFMM_USE_BLAS
	CellClass::Init(DevP);
	std::cout << " $$$$$$$$$$$$$$$  SPHERICAL VERSION $$$$$$$$$$$$"<<std::endl;
#else
	std::cout << " $$$$$$$$$$$$$$$  CHEBYCHEV VERSION $$$$$$$$$$$$" <<std::endl;
#endif
	OctreeClass tree(NbLevels, SizeSubLevels, loader->getBoxWidth(), loader->getCenterOfBox());
	// ---------------------------------------------------------------------------------
	// Insert particles in the Octree
	// ---------------------------------------------------------------------------------   std::cout << "Creating & Inserting " << loader->getNumberOfParticles() << " particles ..." << std::endl;
	std::cout << "\tWidth : " << loader->getBoxWidth() << " \t center x : " << loader->getCenterOfBox().getX()
	    		<< " y : " << loader->getCenterOfBox().getY() << " z : " << loader->getCenterOfBox().getZ() << std::endl;
	std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

	counter.tic();
	FPoint electricMoment(0.0,0.0,0.0) ;
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
	double coeffCorrectionDLPOLY = 2.0*FMath::FPi/volume/3.0 ;

	std::cout << std::endl;
	std::cout << "Total Charge         = "<< totalCharge <<std::endl;
	std::cout << "Electric Moment      = "<< electricMoment <<std::endl;
	std::cout << "Electric Moment norm = "<< dipoleNorm  <<std::endl;
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;
	//
	// ---------------------------------------------------------------------------------
	// write Octree information in octreeData.txt file
	// ---------------------------------------------------------------------------------
	std::ofstream octreeData("octreeData.txt",std::ios::out);
	typename  OctreeClass::Iterator  octreeIterator(&tree);
	octreeIterator.gotoBottomLeft();
	int inTreeHeight = NbLevels ;
	double  inBoxWidth = loader->getBoxWidth() ;
	FPoint inBoxCenter(loader->getCenterOfBox()) ;
	double  widthAtLeafLevel(inBoxWidth/FReal(1 << (inTreeHeight-1))) , widthAtLeafLevelDiv2 = widthAtLeafLevel/2;
	FPoint  boxCorner(inBoxCenter.getX()-(inBoxWidth/2),inBoxCenter.getY()-(inBoxWidth/2),
			inBoxCenter.getZ()-(inBoxWidth/2));
	octreeData << "Box Width  " << inBoxWidth << std::endl ;
	octreeData << "Leaf width " << widthAtLeafLevel << std::endl ;
	octreeData << "Box Corner "<< boxCorner << std::endl<< std::endl ;

	do{
		auto * const FRestrict cell = octreeIterator.getCurrentCell();
		FTreeCoordinate coordinate  = cell->getCoordinate() ;
		FPoint leafCenter(FReal(coordinate.getX()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getX(),
				FReal(coordinate.getY()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getX(),
				FReal(coordinate.getZ()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getX());
		octreeData << "Leaf " << cell->getMortonIndex() << std::endl
				<< "    Center   "<<  coordinate <<    std::endl
				<< "    Center   "<<  leafCenter
				<< std::endl ;
	} while(octreeIterator.moveRight());
	//
	FIOVtk vtkfile ;
	vtkfile.writeOctree("octreeFile.vtk","Octree ", tree) ;
	//
	// ---------------------------------------------------------------------------------

	std::cout << "Create kernel & run simu ..." << std::endl;
	counter.tic();

	FTreeCoordinate min, max;

	if( FParameters::existParameter(argc, argv, "-noper") ){
#ifndef  ScalFMM_USE_BLAS
		KernelClass kernels( DevP, NbLevels, loader->getBoxWidth(), loader->getCenterOfBox());
#else
		KernelClass kernels( NbLevels, loader->getBoxWidth(), loader->getCenterOfBox());
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
		KernelClass kernels(algo.extendedTreeHeight(), algo.extendedBoxWidth(),algo.extendedBoxCenter());
#endif
		algo.setKernel(&kernels);
		algo.execute();
		algo.repetitionsIntervals(&min, &max);
	}

	counter.tac();

	std::cout << "Done  " << "(@FMM Algorithm = " << counter.elapsed() << "s)." << std::endl;
	std::cout << "\n"<< "END FMM "
			<< "-------------------------------------------------------------------------" << std::endl << std::endl ;

	// ----------------------------------------------------------------------------------------------------------
	//                                  DIRECT COMPUTATION
	// ----------------------------------------------------------------------------------------------------------
	EwalParticle* particlesDirect = nullptr;
	FReal directEnergy = 0.0;

	//
	// direct computation
	//
	bool direct = false ;
	if(FParameters::existParameter(argc, argv, "-direct")){
		direct = true ;
		printf("Compute direct:\n");
		printf("Box [%d;%d][%d;%d][%d;%d]\n", min.getX(), max.getX(), min.getY(),
				max.getY(), min.getZ(), max.getZ());

		particlesDirect = new EwalParticle[loader->getNumberOfParticles()];

		FReal denergy = 0.0;
		FMath::FAccurater dfx, dfy, dfz ;
		//
		//
		// particles       : Ewald results
		// particlesDirect : Direct results
		//
		counter.tic();
		for(int idxTarget = 0 ; idxTarget < loader->getNumberOfParticles() ; ++idxTarget){
			particlesDirect[idxTarget] = particles[idxTarget];
			EwalParticle & part        = particlesDirect[idxTarget];
			part.forces[0] = part.forces[1] = part.forces[2] = 0.0;
			part.potential = 0.0;
            //
			// compute with all other except itself
			//
			// Compute force and potential between  particles[idxTarget] and particles inside the box
			//
			for(int idxOther = 0; idxOther < loader->getNumberOfParticles() ; ++idxOther){
				if( idxOther != idxTarget ){
					FP2P::NonMutualParticles(
							particles[idxOther].position.getX(), particles[idxOther].position.getY(),
							particles[idxOther].position.getZ(),particles[idxOther].physicalValue,
							part.position.getX(), part.position.getY(),
							part.position.getZ(),part.physicalValue,
							&part.forces[0],&part.forces[1],
							&part.forces[2],&part.potential);
				}
			}
			//
			// compute with image boxes
			//
			for(int idxX = min.getX() ; idxX <= max.getX() ; ++idxX){
				for(int idxY = min.getY() ; idxY <= max.getY() ; ++idxY){
					for(int idxZ = min.getZ() ; idxZ <= max.getZ() ; ++idxZ){

						{
							//int nL = PeriodicDeep ;
							//	for(int idxX = -nL ; idxX <= nL -1; ++idxX){
							//	  for(int idxY = -nL ; idxY <=nL-1 ; ++idxY){
							//	    for(int idxZ = -nL ; idxZ <= nL-1; ++idxZ){
							if(idxX == 0 && idxY == 0 && idxZ == 0) continue;

							const FPoint offset(loader->getBoxWidth() * FReal(idxX),
									loader->getBoxWidth() * FReal(idxY),
									loader->getBoxWidth() * FReal(idxZ));
							//							std::cout <<" ( "<< idxX<<" , "<<idxY << " , "<< idxZ << " ) "<< offset <<std::endl;
							for(int idxSource = 0 ; idxSource < loader->getNumberOfParticles() ; ++idxSource){
								EwalParticle source = particles[idxSource];
								source.position += offset;
								//								std::cout << "Part "<<idxSource<< " " <<source.position.getX()<< " " << source.position.getY()<< " " <<source.position.getZ()<< " " <<source.physicalValue <<std::endl ;
								FP2P::NonMutualParticles(
										source.position.getX(), source.position.getY(),source.position.getZ(),source.physicalValue,
										part.position.getX(), part.position.getY(),part.position.getZ(),part.physicalValue,
										&part.forces[0],&part.forces[1],&part.forces[2],&part.potential
								);
							}
							//						std::cout <<std::endl<<std::endl<<std::endl;
						}
					}
				} // END
			}
			//
			// remove dipole correction for DL_POLY
			//
			part.forces[0] -= 2.0*part.physicalValue*coeffCorrectionDLPOLY*electricMoment.getX() ;
			part.forces[1] -= 2.0*part.physicalValue*coeffCorrectionDLPOLY*electricMoment.getY() ;
			part.forces[2] -= 2.0*part.physicalValue*coeffCorrectionDLPOLY*electricMoment.getZ() ;
			//
			//  Scaling
			//
			part.forces[0] *= scaleForce ;
			part.forces[1] *= scaleForce ;
			part.forces[2] *= scaleForce ;
			if(FParameters::existParameter(argc, argv, "-verbose")){
				std::cout << ">> index " << particles[idxTarget].index << std::endl;
				std::cout << "Good    x " << particles[idxTarget].position.getX() << " y " << particles[idxTarget].position.getY() << " z " << particles[idxTarget].position.getZ() << std::endl;
				std::cout << "Good    fx " <<particles[idxTarget].forces[0] << " fy " << particles[idxTarget].forces[1] << " fz " << particles[idxTarget].forces[2] << std::endl;
				std::cout << "DIRECT  fx " << part.forces[0]<< " fy " <<part.forces[1] << " fz " << part.forces[2] << std::endl;
				std::cout << "ratio  fx " << particles[idxTarget].forces[0]/part.forces[0] << " fy " <<particles[idxTarget].forces[1]/part.forces[1] << " fz " << particles[idxTarget].forces[2]/part.forces[2] << std::endl;
				std::cout << "DIRECT  physical value " << part.physicalValue << " potential " << part.potential << std::endl;

				std::cout << "\n";
			}
			//
			dfx.add(particles[idxTarget].forces[0],part.forces[0]);
			dfy.add(particles[idxTarget].forces[1],part.forces[1]);
			dfz.add(particles[idxTarget].forces[2],part.forces[2]);

			denergy += part.potential * part.physicalValue;

			//	particlesDirect[idxTarget] = part;
		}
		counter.tac();
		denergy *= 0.5*scaleEnergy ;
		denergy -= coeffCorrectionDLPOLY*dipoleNorm*scaleEnergy ;
		directEnergy	= denergy ;
		printf("Difference between direct and Ewald (DL_POLY)\n");
		printf("DirectEwald Fx diff is = \n");
		printf("DirectEwald   L2Norm  %e\n",dfx.getRelativeL2Norm() );
		printf("DirectEwald   InfNorm %e\n",dfx.getRelativeInfNorm());
		printf("DirectEwaldFy diff is = \n");
		printf("DirectEwald   L2Norm  %e\n",dfx.getRelativeL2Norm());
		printf("DirectEwald   InfNorm %e\n",dfy.getRelativeInfNorm());
		printf("DirectEwaldFz diff is = \n");
		printf("DirectEwald   L2Norm  %e\n",dfx.getRelativeL2Norm());
		printf("DirectEwald   InfNorm %e\n",dfz.getRelativeInfNorm());
		double L2error = (dfx.getRelativeL2Norm()*dfx.getRelativeL2Norm() + dfy.getRelativeL2Norm()*dfy.getRelativeL2Norm()  + dfz.getRelativeL2Norm() *dfz.getRelativeL2Norm()  );
		printf("DirectEwald L2 Force Error= %e\n",FMath::Sqrt(L2error)) ;
		printf("DirectEwald Energy DIRECT=   %.12e\n",denergy);
		printf("DirectEwald Energy EWALD =   %.12e\n",loader->getEnergy());
		printf("DirectEwald Energy EWALD - Energy DIRECT              = %e\n", loader->getEnergy()-denergy );
		printf("DirectEwald |Energy EWALD - Energy DIRECT|/directEnergy= %e\n",FMath::Abs((loader->getEnergy()-denergy)/loader->getEnergy()));

		//
		if(FParameters::existParameter(argc, argv, "-saveError")){
			errorfile << std::endl << "      END DIRECT " << std::endl ;
		}
		std::cout << "Done  " << "(@DIRECT Algorithm = " << counter.elapsed() << "s)." << std::endl;

		std::cout << "\n"<< "END DIRECT "
				<< "-------------------------------------------------------------------------" << std::endl << std::endl ;
	}
	// ----------------------------------------------------------------------------------------------------------
	//                                  DIRECT -- FMM  Comparisons
	//                                  Ewald  -- FMM  Comparisons
	// ----------------------------------------------------------------------------------------------------------
	{
		// particles       : Ewald results
		// particlesDirect : Direct results
		//

		FReal energy = 0.0;
		FMath::FAccurater fx, fy, fz, fmmdfx, fmmdfy, fmmdfz,fmmpot;

		tree.forEachLeaf([&](LeafClass* leaf){
			const FReal*const potentials = leaf->getTargets()->getPotentials();
			FReal*const forcesX = leaf->getTargets()->getForcesX();
			FReal*const forcesY = leaf->getTargets()->getForcesY();
			FReal*const forcesZ = leaf->getTargets()->getForcesZ();
			const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
			const FReal*const positionsX = leaf->getTargets()->getPositions()[0];
			const FReal*const positionsY = leaf->getTargets()->getPositions()[1];
			const FReal*const positionsZ = leaf->getTargets()->getPositions()[2];
			const int nbParticlesInLeaf  = leaf->getTargets()->getNbParticles();
			const FVector<int>& indexes  = leaf->getTargets()->getIndexes();
			//
			for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
				const int indexPartOrig = indexes[idxPart];
				//
				// remove dipole correction for DL_POLY
				//
				forcesX[idxPart] -= 2.0*physicalValues[idxPart]*coeffCorrectionDLPOLY*electricMoment.getX() ;
				forcesY[idxPart] -= 2.0*physicalValues[idxPart]*coeffCorrectionDLPOLY*electricMoment.getY() ;
				forcesZ[idxPart] -= 2.0*physicalValues[idxPart]*coeffCorrectionDLPOLY*electricMoment.getZ() ;
				//
				forcesX[idxPart] *= scaleForce;
				forcesY[idxPart] *= scaleForce;
				forcesZ[idxPart] *= scaleForce;
				if(direct) {
//					std::cout << "Good    x " << particlesDirect[idxPart].position.getX() << " y " << particlesDirect[idxPart].position.getY() << " z " << particlesDirect[idxPart].position.getZ() << std::endl;
//					std::cout << "Direct fx " <<particlesDirect[idxPart].forces[0]<< " fy " << particlesDirect[idxPart].forces[1] << " fz " << particlesDirect[idxPart].forces[2] << std::endl;
//					std::cout << "FMM  fx " << forcesX[idxPart] << " fy " <<forcesY[idxPart] << " fz " << forcesZ[idxPart] << std::endl;

					fmmdfx.add(particlesDirect[indexPartOrig].forces[0],forcesX[idxPart]);
					fmmdfy.add(particlesDirect[indexPartOrig].forces[1],forcesY[idxPart]);
					fmmdfz.add(particlesDirect[indexPartOrig].forces[2],forcesZ[idxPart]);
					fmmpot.add(particlesDirect[indexPartOrig].potential,potentials[idxPart]);
				} else {
					fx.add(particles[indexPartOrig].forces[0],forcesX[idxPart]); // Ewald - FMM
					fy.add(particles[indexPartOrig].forces[1],forcesY[idxPart]);
					fz.add(particles[indexPartOrig].forces[2],forcesZ[idxPart]);
				}
				energy +=  potentials[idxPart]*physicalValues[idxPart];
				//
				if(FParameters::existParameter(argc, argv, "-verbose")){
					std::cout << ">> index " << particles[indexPartOrig].index << std::endl;
					std::cout << "Good x " << particles[indexPartOrig].position.getX() << " y " << particles[indexPartOrig].position.getY() << " z " << particles[indexPartOrig].position.getZ() << std::endl;
					std::cout << "FMM  x " << positionsX[idxPart] << " y " << positionsY[idxPart] << " z " << positionsZ[idxPart] << std::endl;
					std::cout << "Good fx " <<particles[indexPartOrig].forces[0] << " fy " << particles[indexPartOrig].forces[1] << " fz " << particles[indexPartOrig].forces[2] << std::endl;
					std::cout << "FMM  fx " << forcesX[idxPart] << " fy " <<forcesY[idxPart] << " fz " << forcesZ[idxPart] << std::endl;
					std::cout << "ratio  fx " << particles[indexPartOrig].forces[0]/forcesX[idxPart] << " fy " <<particles[indexPartOrig].forces[1]/forcesY[idxPart] << " fz " << particles[indexPartOrig].forces[2]/forcesZ[idxPart] << std::endl;
//					std::cout << "GOOD physical value " << particles[indexPartOrig].physicalValue << " potential " << particles[indexPartOrig].potential << std::endl;
//					std::cout << "FMM  physical value " << physicalValues[idxPart] << " potential " << potentials[idxPart] <<  " energy cumul " << energy<<std::endl;

					std::cout << "\n";
				}
				if(FParameters::existParameter(argc, argv, "-saveError")){
					double ratio,tmp ;
					ratio = std::fabs(1-particles[indexPartOrig].forces[0]/forcesX[idxPart]) ;
					tmp    = std::fabs(1-particles[indexPartOrig].forces[1]/forcesY[idxPart]) ;
					ratio = std::max(tmp, ratio) ;
					ratio = std::max(std::fabs(1-particles[indexPartOrig].forces[2]/forcesZ[idxPart]) , ratio) ;
					if(ratio >=0.25) {
						errorfile << ">> index " << particles[indexPartOrig].index << std::endl;
						errorfile << "Good x " << particles[indexPartOrig].position.getX() << " y " << particles[indexPartOrig].position.getY() << " z " << particles[indexPartOrig].position.getZ() << std::endl;
						errorfile << "FMM  x " << positionsX[idxPart] << " y " << positionsY[idxPart] << " z " << positionsZ[idxPart] << std::endl;
						errorfile << "Good fx " <<particles[indexPartOrig].forces[0] << " fy " << particles[indexPartOrig].forces[1] << " fz " << particles[indexPartOrig].forces[2] << std::endl;
						errorfile << "FMM  fx " << forcesX[idxPart] << " fy " <<forcesY[idxPart] << " fz " << forcesZ[idxPart] << std::endl;
						errorfile << "ratio  fx " << particles[indexPartOrig].forces[0]/forcesX[idxPart] << " fy " <<particles[indexPartOrig].forces[1]/forcesY[idxPart] << " fz " << particles[indexPartOrig].forces[2]/forcesZ[idxPart] << std::endl;
						errorfile << "GOOD physical value " << particles[indexPartOrig].physicalValue << " potential " << particles[indexPartOrig].potential << std::endl;
						errorfile << "FMM  physical value " << physicalValues[idxPart] << " potential " << potentials[idxPart] << std::endl;
						errorfile << "\n";
					}
				}
			}
		});
		energy *= 0.5*scaleEnergy ;
		energy -= coeffCorrectionDLPOLY*dipoleNorm*scaleEnergy ;
		//    printf("(Energy EWALD - Energy FMM)/dipoleNorm = %e\n",(loader->getEnergy()-energy)*volume/(dipoleNorm));
		//    printf("(dipoleNorm /Volume = %e\n",correctEnergy);
		//
		if(direct) {
			printf("Difference between FMM and DiRECT:\n");
			printf("DirectFmm Fx diff is = \n");
			printf("DirectFmm   L2Norm  %e\n",fmmdfx.getRelativeL2Norm());
			printf("DirectFmm   InfNorm %e\n",fmmdfx.getRelativeInfNorm());
			printf("DirectFmm Fy diff is = \n");
			printf("DirectFmm   L2Norm  %e\n",fmmdfy.getRelativeL2Norm());
			printf("DirectFmm   InfNorm %e\n",fmmdfy.getRelativeInfNorm());
			printf("DirectFmm Fz diff is = \n");
			printf("DirectFmm   L2Norm  %e\n",fmmdfz.getRelativeL2Norm());
			printf("DirectFmm   InfNorm %e\n",fmmdfz.getRelativeInfNorm());
			double L2error = (fmmdfx.getRelativeL2Norm()*fmmdfx.getRelativeL2Norm() + fmmdfy.getRelativeL2Norm()*fmmdfy.getRelativeL2Norm()  + fmmdfz.getRelativeL2Norm() *fmmdfz.getRelativeL2Norm()  );
			printf("DirectFmm L2 Force Error= %e\n",FMath::Sqrt(L2error)) ;
			printf("DirectFmm L2 potential Error= %e\n",fmmpot.getRelativeL2Norm()) ;

			//
			printf("DirectFmm Energy FMM=    %.12e\n",energy);
			printf("DirectFmm Energy DIRECT= %.12e\n",directEnergy);
			printf("DirectFmm |Energy DIRECT - Energy DIRECT|/directEnergy= %e\n",FMath::Abs((directEnergy-energy)/directEnergy));
		}
		else {
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
			printf("FmmEwald Energy FMM=   %.12e\n",energy);
			printf("FmmEwald Energy EWALD= %.12e\n",loader->getEnergy());
			printf("FmmEwald Energy EWALD - Energy FMM = %e\n",loader->getEnergy()-energy);
			printf("FmmEwald |Energy EWALD -Energy FMM|/Energy EWALD= %e\n",FMath::Abs((loader->getEnergy()-energy)/loader->getEnergy()));

		}
		printf("vtkParticles\n");

		vtkfile.writeParticules(tree, true,1,1,"vtkParticles","FMM results");

	}
	//
	// -----------------------------------------------------
	// ReGenerate file
	//
	if( FParameters::existParameter(argc, argv, "-gen") ){
		std::cout << "Generate ewal.out from input file" << std::endl;
		std::ofstream fileout("ewal.out",std::ifstream::out);
		std::ifstream file(filename,std::ifstream::in);
		if(file.is_open()){
			const int bufferSize = 512;
			char buffer[bufferSize];
			file.getline(buffer, bufferSize);
			fileout << buffer << '\n';
			file.getline(buffer, bufferSize);
			fileout << buffer << '\n';
			if( !FParameters::existParameter(argc, argv, "-noper") ){
				file.getline(buffer, bufferSize);
				fileout << buffer << '\n';
				file.getline(buffer, bufferSize);
				fileout << buffer << '\n';
				file.getline(buffer, bufferSize);
				fileout << buffer << '\n';
			}
			for(int idxPart = 0 ; idxPart < loader->getNumberOfParticles() ; ++idxPart){
				file.getline(buffer, bufferSize);
				fileout << buffer << '\n';
				file.getline(buffer, bufferSize);
				fileout << buffer << '\n';
				file.getline(buffer, bufferSize);
				file.getline(buffer, bufferSize);
			}
		}
	}
	// end generate

	// -----------------------------------------------------

	delete[] particles;
	delete[] particlesDirect;

	return 0;
}



