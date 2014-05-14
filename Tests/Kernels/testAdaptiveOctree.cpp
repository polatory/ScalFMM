// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas
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
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>

#include "Utils/FParameters.hpp"
#include "Utils/FMath.hpp"
#include "Files/FFmaGenericLoader.hpp"
//
#include "Components/FSimpleIndexedLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
//
#include "Kernels/Rotation/FRotationKernel.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"
//
#include "Containers/FOctree.hpp"
//
#include "AdaptiveTree/FAdaptCell.hpp"
#include "AdaptiveTree/FAdaptTools.hpp"
#include "AdaptiveTree/FAdaptSeqAlgorithm.hpp"


/// \file  testadaptiveOnOctree.cpp
//!
//! \brief Driver to generate statistics on the octree for the particle distribution given by parameter -infile
//!  \authors B. Bramas, O. Coulaud
//!
//!  This code gives you some statistics (Particles, Cells, P2P and M2L number of operators) on the octree
//!    Those statistics are shown level by level.
//!
//! The statistics are
//!
//!  <b>for particles:</b> min/max particles per leaf, the average number, the variance and the average number of P2P neighbors.
//!
//! <b>For each level in the octree</b>
//!
//!  \arg       The number of cells, of adaptive cells (cell with more one child)
//!  \arg          The average, min and max numbers of M2L operators and also its  variance.

//!
//!  <b> General arguments:</b>
//!     \param   -help(-h)      to see the parameters available in this driver
//!     \param   -depth    The depth of the octree
//!     \param   -subdepth          Specifies the size of the sub octree
//!
//!     \param   -infile name   Name of the particles file. The file have to be in our FMA format
//!     \param    -bin              if the input file in binary mode
//!     \param   -outfile name Generic name  for output file  (without extension)
//!
//!  <b> Statistics options:</b>
//!   \param -stat to build the statistivs on the octree
//!   \param -histP build a file to generate histogram of particles per leaf. The data are store in file given by -outfile arguments and  .txt extension. (only if -stat is set)
//   \param   -sM    s_min^M threshold for Multipole expansion (l+1)^2 for Spherical harmonics"
//   \param   -sL    s_min^M threshold for local expansion  (l+1)^2 for Spherical harmonics"
// Simply create particles and try the kernels
//
void usage() {
	std::cout << "Driver to obtain statistics on the octree" << std::endl;
	std::cout <<	 "Options  "<< std::endl
			<<     "      -help       to see the parameters    " << std::endl
			<<	     "      -depth        the depth of the octree   "<< std::endl
			<<	     "      -subdepth   specifies the size of the sub octree   " << std::endl
			<<     "      -infile name specifies the name of the particle distribution" << std::endl
			<<    "        -bin              if the input file in binary mode"<< std::endl
			<<     "      -outfile name  specifies the file for the diagnostics" << std::endl
			<<     "      -histP   build the histogram of the particle number per leaf"<<std::endl
			<<     "      -sM    s_min^M threshold for Multipole (l+1)^2 for Spherical harmonics"<<std::endl;
}

int main(int argc, char ** argv){
	//

	// accuracy
	const unsigned int P =8;
	//
	typedef FP2PParticleContainerIndexed<>                             ContainerClass;
	typedef FSimpleIndexedLeaf<ContainerClass>                     LeafClass;
	typedef FRotationCell<P>                                                     KernelCellClass;
	//
	typedef FAdaptCell<KernelCellClass,LeafClass>                         CellClass;
	typedef FOctree<CellClass, ContainerClass, LeafClass >             OctreeClass;
	typedef FRotationKernel< KernelCellClass, ContainerClass , P>   KernelClass;
	//
	typedef FAdaptSeqAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
	//
	if(FParameters::existParameter(argc, argv, "-h")||FParameters::existParameter(argc, argv, "-help")|| (argc < 3 )){
		usage() ;
		exit(-1);
	}
	//
	//   Octree parameters
	//
	const int NbLevels        = FParameters::getValue(argc,argv,"-depth", 5);
	const int SizeSubLevels = FParameters::getValue(argc,argv,"subdepth", 3);
	const int sminM            = FParameters::getValue(argc,argv,"-sM", 0);
	const int sminL             = FParameters::getValue(argc,argv,"-sL", 0);
	//
	//  input and output  Files parameters
	//
	const char* const filename = FParameters::getStr(argc,argv,"-infile", "../Data/test20k.fma");
	const std::string genericFileName(FParameters::getStr(argc,argv,"-outfile",   "output"));
	//
	std::cout << "Opening : " << filename << "\n";
	bool binaryMode = false;
	if(FParameters::existParameter(argc, argv, "-bin")){
		binaryMode = true;
	}
	FFmaGenericLoader loader(filename,binaryMode);
	if(!loader.isOpen()){
		std::cout << "Loader Error, " << filename << " is missing\n";
		return 1;
	}
	//
	// -----------------------------------------------------
	OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
	//
	// -----------------------------------------------------
	//     Creating and Inserting particles in the tree
	// -----------------------------------------------------
	//
	std::cout << "Tree box is cubic "<<std::endl
			<< "         Centre:   "<< loader.getCenterOfBox() <<std::endl
			<< "         Length:  "<< loader.getBoxWidth()       <<std::endl <<std::endl;
	std::cout 		<< "         criteria SM:  "<< sminM     <<std::endl
	 		<< "         criteria SL:  "<< sminL     <<std::endl <<std::endl;
	//
	std::cout << "Creating and Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
	std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
	FPoint particlePosition, minPos, maxPos;
	FReal physicalValue;
	for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
		loader.fillParticle(&particlePosition,&physicalValue);
		//
		minPos.setX(FMath::Min(minPos.getX(),particlePosition.getX())) ;
		minPos.setY(FMath::Min(minPos.getY(),particlePosition.getY())) ;
		minPos.setZ(FMath::Min(minPos.getZ(),particlePosition.getZ())) ;
		maxPos.setX(FMath::Max(maxPos.getX(),particlePosition.getX())) ;
		maxPos.setX(FMath::Max(maxPos.getY(),particlePosition.getY())) ;
		maxPos.setX(FMath::Max(maxPos.getZ(),particlePosition.getZ())) ;
		//
		tree.insert(particlePosition, idxPart, physicalValue );
	}
	std::cout << "Data are inside the box delimited by "<<std::endl
			<< "         Min corner:  "<< minPos<<std::endl
			<< "         Max corner:  "<< maxPos<<std::endl <<std::endl;
	//
	OctreeClass::Iterator octreeIterator(&tree);

	//
	// -----------------------------------------------------
	//     Build information for adaptive tree
	// -----------------------------------------------------
	//
	{
		// Set global If for debug purpose
		long int idCell  = setGlobalID(tree);

		std::cout << " start build smin criteria " <<std::endl;
		//
		adaptiveTreeBuilSminC(tree,sminM,sminL) ;
		//
		//  Set Global id
		//
		//		long int idCell  = setGlobalID(tree);
		//
		//  Build CA and FA  lists
		std::cout << " start building CA and FA lists " <<std::endl;
		//
		adaptiveTreeBuildLists(tree) ;
		//
	}
	//
	////////////////////////////////////////////////////////////////////
	//             Execute adaptive algorithm
	////////////////////////////////////////////////////////////////////
	//
	// -----------------------------------------------------
	//
	KernelClass *kernels = new KernelClass(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
	//
	FmmClass algorithm(&tree, kernels);
	algorithm.execute();
	//
	////////////////////////////////////////////////////////////////////
	//              Export adaptive tree in tulip format
	////////////////////////////////////////////////////////////////////
	//
	// -----------------------------------------------------
	//

	//
	//  Set Global id for tulip export
	//
	long int idCell  = setGlobalID(tree);
	//
	//
	std::cout << " start export tulip " <<std::endl;

	//
	// Set Global indexes to save the octree in tulip format
	//
	// -----------------------------------------------------
	std::ofstream tlp("aa.tlp", std::ofstream::out );

	TulipExport( tlp, idCell, tree);

	std::cout << " NVCells " << idCell <<  std::endl ;


	octreeIterator.gotoTop() ;
	for(int idxLevel = 1 ; idxLevel < NbLevels ;  ++idxLevel){
		std::cout << "idxLevel: "<<idxLevel << "  Iterator Level    " << octreeIterator.level()<<  "  is leaves level: " << octreeIterator.isAtLeafLevel()	<<std::endl;
		octreeIterator.moveDown() ;octreeIterator.gotoLeft();
	}
	std::cout << "Level max " <<  NbLevels <<std::endl;



	std::ofstream file("aa.tree", std::ofstream::out );
	//
	octreeIterator.gotoTop() ;  // here we are at level 1 (first child)
	//
	////////////////////////////////////////////////////////////////////
	//              Export adaptive tree in our format
	////////////////////////////////////////////////////////////////////
	//
	// -----------------------------------------------------
	//
	std::cout << "Top of the octree " << octreeIterator.level() << std::endl ;
	for(int idxLevel = 1 ; idxLevel < NbLevels ;  ++idxLevel){
		file << std::endl << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<< std::endl;
		file << "  Level " << idxLevel <<" OLevel  "<<  octreeIterator.level()<<  "  -- leave level " <<   std::boolalpha <<  octreeIterator.isAtLeafLevel() << std::endl;
		do{
			file << *(octreeIterator.getCurrentCell())<< std::endl ;
		} while(octreeIterator.moveRight());
		octreeIterator.moveDown() ;
		octreeIterator.gotoLeft();
	}
	std::cout << "   END    " << std::endl;

	// Check
	octreeIterator.gotoBottomLeft();
	do {
		std::cout << " Level " <<octreeIterator.level() <<std::endl;
	}while(octreeIterator.moveUp() );
	//
	return 0;
}



