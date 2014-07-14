// ===================================================================================
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
#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Files/FFmaGenericLoader.hpp"
#include "Kernels/P2P/FP2P.hpp"
#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Components/FBasicCell.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Containers/FOctree.hpp"


//
/// \file  DirectComputation.cpp
//!
//! \brief DirectComputation: Driver to compute direct interaction between N particles for 1/r kernel.
//!
//! DirectComputation: Driver to compute direct interaction between N particles for 1/r kernel.
//! the particles are read from file given by -fin argument and potential, forces are stored in FMA format.
//!  <b> General arguments:</b>
//!     \param   -help (-h)      to see the parameters available in this driver
//!     \param   -fin name:  file name  to convert (with extension .fma (ascii) or bfma (binary).
//!                             Only our FMA (.bma, .bfma) is allowed "
//!     \param    -fout filenameOUT   output file  with extension (default output.bfma)
//!      \param   -verbose : print index x y z Q V fx fy fz
//!

// Simply create particles and try the kernels
int main(int argc, char ** argv){
	//
	///////////////////////What we do/////////////////////////////
	if( FParameters::existParameter(argc, argv, "-help" ) || argc < 4){
		std::cout << ">> This executable has to be used to compute  interaction either for periodic or non periodic system.\n";
		std::cout << ">> Example -fin filenameIN.{fma or bfma)     -fout filenameOUT{fma or bfma) \n";
		std::cout << ">> Default input file : ../Data/unitCubeXYZQ20k.fma\n";
		std::cout << " Options " << std::endl;
		std::cout << "     -verbose : print index x y z Q V fx fy fz " << std::endl;
		std::cout << "     -fin filename. Extension specifies if the file is binary or not. " << std::endl;
		std::cout << "                            Only our FMA (.bma, .bfma) is allowed " << std::endl;
		std::cout << "     -fout filenameOUT   output file  with extension (default output.bfma)" << std::endl;
		exit(-1);

	}

	//////////////////////////////////////////////////////////////

	const std::string defaultFile(/*SCALFMMDataPath+*/"../Data/unitCubeXYZQ20k.fma");
	const std::string filenameIn(FParameters::getStr(argc,argv,"-fin",  defaultFile.c_str()));
	const std::string filenameOut(FParameters::getStr(argc,argv,"-fout", "output.bfma"));
	//
	FTic counter;

	// -----------------------------------------------------
	//  LOADER
	//  -----------------------------------------------------
	// ---------------------------------------------------------------------------------
	// Read  particles in the Octree
	// ---------------------------------------------------------------------------------
	std::cout << "Opening : " << filenameIn << "\n";
	//
	FFmaGenericLoader loader(filenameIn);
	//
	int nbParticles = static_cast<int>(loader.getNumberOfParticles());
	std::cout << "Read " << nbParticles << " particles ..." << std::endl;
	double BoxWith=loader.getBoxWidth();
	FPoint Centre(loader.getCenterOfBox().getX(), loader.getCenterOfBox().getY() , loader.getCenterOfBox().getZ());
	std::cout << "\tWidth : " <<BoxWith << " \t center x : " << loader.getCenterOfBox().getX()
																		<< " y : " << loader.getCenterOfBox().getY() << " z : " << loader.getCenterOfBox().getZ() << std::endl;

	counter.tic();
	FmaRParticle *  particles = new FmaRParticle[nbParticles];
	memset(particles, 0, sizeof(FmaRParticle) * nbParticles) ;
	//
	double totalCharge = 0.0;
	//
//	int nbDataToRead = particles[0].getReadDataNumber();
	for(int idx = 0 ; idx<nbParticles ; ++idx){
		//
	    loader.fillParticle(particles[idx].getPtrFirstData(), particles[idx].getReadDataNumber());
		//	loader.fillParticle(particles[idx].getPtrFirstData(), nbDataToRead);    // OK
		//  loader.fillParticle(particles[idx]); // OK
	//    std::cout << idx <<"  "<<  particles[idx].getPosition() << " "<<particles[idx].getPhysicalValue() << " "<<particles[idx].getPotential()
	//			<<"  " << particles[idx].getForces()[0]<<"  " <<particles[idx].getForces()[1]<<"  " <<particles[idx].getForces()[2]<<"  " <<std::endl;
		//
	    totalCharge += particles[idx].getPhysicalValue() ;
	}

	counter.tac();

	std::cout << std::endl;
	std::cout << "Total Charge         = "<< totalCharge <<std::endl;
	std::cout << std::endl;

	std::cout << "Done  " << "(@ reading Particles  " << counter.elapsed() << " s)." << std::endl;
	//Need to copy particles to ContainerClass -> FP2PParticleContainer
	typedef FBasicCell                 CellClass;
	typedef FP2PParticleContainer<>         ContainerClass;

	typedef FSimpleLeaf< ContainerClass >                     LeafClass;
	typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;

	OctreeClass tree(2, 1, loader.getBoxWidth(), loader.getCenterOfBox());
	for(int idxP=0 ; idxP<nbParticles ; ++idxP){
	    tree.insert(particles[idxP].getPosition());
	}
	int i=0;
	tree.forEachLeaf([&](LeafClass * leaf){
		printf("leaf : %d\n",i++ );
	    });
	//
	// ----------------------------------------------------------------------------------------------------------
	//                                   COMPUTATION
	// ----------------------------------------------------------------------------------------------------------
	FReal denergy = 0.0;
	//
	//  computation
	//
	{
	    typename OctreeClass::Iterator iterator(&tree);
	    iterator.gotoBottomLeft();

	    do{
		FTreeCoordinate coord = iterator.getCurrentGlobalCoordinate();
		ContainerClass** neighbors = new ContainerClass*[27];
		tree.getLeafsNeighbors(neighbors,coord,1);
		FP2P::FullMutual(iterator.getCurrentLeaf()->getTargets(),neighbors,27);

	    }while(iterator.moveRight());
	}

	delete[] particles;
	return 0;
}