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
#include "Containers/FOctree.hpp"

#include "AdaptativeTree/FAdaptCell.hpp"
#include "AdaptativeTree/FAdaptTools.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Components/FBasicParticleContainer.hpp"


#include "Utils/FMath.hpp"
#include "Files/FFmaGenericLoader.hpp"


/// \file  statisticsOnOctree.cpp
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
	typedef FBasicParticleContainer<0>                                     ContainerClass;
	typedef FSimpleLeaf<ContainerClass>                                  LeafClass;
	typedef FAdaptCell<FBasicCell,LeafClass>                            CellClass;
	typedef FOctree<CellClass, ContainerClass, LeafClass >       OctreeClass;
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
		tree.insert(particlePosition );
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
		std::cout << " start build smin criteria " <<std::endl;
		//
		adaptativeTreeBuilSminC(tree,sminM,sminL) ;
		//
		//  Set Global id
		//
//		long int idCell  = setGlobalID(tree);
		//
		//  Build CA and FA  lists
		std::cout << " start building CA and FA lists " <<std::endl;
		//
		adaptativeTreeBuildLists(tree) ;
		//
	}
	//
	// -----------------------------------------------------
	//     Start statistics
	// -----------------------------------------------------
	//
	int removeM = 0 ;
	if(FParameters::existParameter(argc, argv, "-stat")){ // get stats
		{    // get stats on the leaf level (Particles)
			long int allLeaves =  (1 << (3* (NbLevels-1) )) ;
			std::cout << std::endl<< "[STAT] Leaf level "  << " is  " << NbLevels << std::endl;
			std::cout << "[STAT] potentials leafs number is " << allLeaves<< std::endl;

			FReal averageParticles = 0.0, varianceParticles = 0.0 ;
			int nbLeafs = 0,minParticles = 1000000.0, maxParticles = 0.0 ;
			//
			// Start to compute statistics on particles
			//
			int  nbPart ,nbTPart=0;
			octreeIterator.gotoBottomLeft();
			do{
				nbPart                   = octreeIterator.getCurrentListTargets()->getNbParticles() ;
				minParticles          =  FMath::Min(minParticles,nbPart) ;
				maxParticles         =  FMath::Max(maxParticles,nbPart) ;
				nbTPart              += nbPart;
				varianceParticles += FReal(nbPart*nbPart) ;
				++nbLeafs;
				if(nbPart < sminM){
					++removeM;
				}
			} while(octreeIterator.moveRight());
			averageParticles   = nbTPart/FReal(nbLeafs);
			varianceParticles  = varianceParticles/FReal(nbLeafs) - averageParticles*averageParticles;
			//
			std::cout.precision(4);
			std::cout << "[STAT] Non empty leafs: " << nbLeafs << " % of non empty leaves: "<<100*static_cast<FReal>(nbLeafs)/static_cast<FReal>(allLeaves)<<" %" << std::endl;
			std::cout << "[STAT]  Particles on leafs:"  << std::endl
					<<   "[STAT]           Min:         "<< minParticles << std::endl
					<<   "[STAT]           Max:        "<< maxParticles << std::endl
					<<   "[STAT]           Average:  "<< averageParticles << std::endl
					<<   "[STAT]           Variance: " << varianceParticles << std::endl;
			std::cout << "[STAT]  number of P2M to remove: " << 		removeM <<	 std::endl;
			//
			//  Histogram of particles per leaf
			//

			if(FParameters::existParameter(argc, argv, "-histP")){
				int size = maxParticles+1;
				int * hist = new int [size] ;
				memset(hist,0,(size)*sizeof(int));
				octreeIterator.gotoBottomLeft();
				do{
					nbPart  = octreeIterator.getCurrentListTargets()->getNbParticles() ;
					++hist[nbPart] ;
				} while(octreeIterator.moveRight());
				//
				// write data
				//
				std::ofstream outfile( genericFileName + ".txt", std::ofstream::out);
				if(!outfile) {
					std::cout << "Cannot open file "<< std::endl;
					exit(-1)	;
				}	//
				outfile << "# Particle histogram.   "<< size << " chunk" <<std::endl;
				for(int i=0 ; i < size ; ++i){
					outfile << i << "  " << hist[i] <<std::endl;
				}
				delete [] hist ;
			}
			FReal averageNeighbors = 0.0, varianceNeighbors =0.0 ;
			int nbBox,minBox=30,maxBox=0;
			octreeIterator.gotoBottomLeft();
			ContainerClass*  neighborsP2P[27];
			do{
				//
				//  P2P Neighbors
				//
				nbBox = tree.getLeafsNeighbors(neighborsP2P, octreeIterator.getCurrentGlobalCoordinate(),NbLevels-1) ;
				// need the current particles and neighbors particles
				minBox                      =  FMath::Min(minBox,nbBox) ;
				maxBox                     =  FMath::Max(maxBox,nbBox) ;
				averageNeighbors   += FReal(nbBox);
				varianceNeighbors  += FReal(nbBox*nbBox) ;
			} while(octreeIterator.moveRight());
			//
			averageNeighbors/=FReal(nbLeafs) ;
			varianceNeighbors = varianceNeighbors/nbLeafs-averageNeighbors*averageNeighbors;
			//
			std::cout << "[STAT]  P2P Neighbors for each leaf " << std::endl
					<< "[STAT]           Min:       " <<  minBox << std::endl
					<< "[STAT]           Max:      " <<  maxBox << std::endl
					<< "[STAT]           Average: " <<  averageNeighbors<< std::endl
					<< "[STAT]           Variance: " <<  varianceNeighbors << std::endl<< std::endl;
		}
		//
		//    ---------------  END LEAVES ---------------
		//
		{

			long long int totalCells = 0;
			long long int totalM2L = 0;
			long long int totalM2ML2L = 0;

			int nbCellsAtTop = 0;
			int nbCellsAtBottom = 0;

			octreeIterator.gotoBottomLeft();

			for(int idxLevel = NbLevels - 1 ; idxLevel >= 1 ; --idxLevel){
				removeM =0 ;
				int nbCellsAtLevel = 0;
				int nbChildAtLevel = 0, adaptiveCell=0 ,nbChildForMyCell;
				int nbNeighborsAtLevel = 0;
				//
				int nbM2LNeighbors, minM2L=500,maxM2L=-1;
				FReal averageM2LNeighbors=0.0, varianceM2LNeighbors=0.0	;
				//
				const CellClass* neighborsM2L[343];
				do{
					++nbCellsAtLevel;
					// Check number of
					if( idxLevel != NbLevels - 1 ){
						nbChildForMyCell=0 ;
						auto** child = octreeIterator.getCurrentChild();
						auto  &  myCell = *(octreeIterator.getCurrentCell());
						int nbPart = 0 ;
						//						std::cout << "NB: ";
						for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
							if(child[idxChild]) {
								++nbChildForMyCell;
								nbPart += child[idxChild]->getnbPart();
								//								std::cout << "  "<< child[idxChild]->getnbPart();
							}
						}
						//						std::cout << std::endl;
						octreeIterator.getCurrentCell()->addPart(nbPart);
						if(octreeIterator.getCurrentCell()->getnbPart() < sminM){
							++removeM;
						}
						nbChildAtLevel += nbChildForMyCell ;
						if(nbChildForMyCell>1) {
							++adaptiveCell ;
						}
						else
						{myCell.setCellAdaptative();}
					}
					const CellClass* neighbors[343];
					nbNeighborsAtLevel += tree.getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),idxLevel);
					//
					//  M2L Neighbors
					//
					nbM2LNeighbors = tree.getInteractionNeighbors(neighborsM2L, octreeIterator.getCurrentGlobalCoordinate(),idxLevel);
					minM2L                          =  FMath::Min(minM2L,nbM2LNeighbors) ;
					maxM2L                         =  FMath::Max(maxM2L,nbM2LNeighbors) ;
					averageM2LNeighbors  += FReal(nbM2LNeighbors) ;
					varianceM2LNeighbors += FReal(nbM2LNeighbors*nbM2LNeighbors) ;
				} while(octreeIterator.moveRight());
				//
				averageM2LNeighbors/=FReal(nbCellsAtLevel) ;
				varianceM2LNeighbors = varianceM2LNeighbors/nbCellsAtLevel-averageM2LNeighbors*averageM2LNeighbors;

				std::cout << "[STAT] Level = " << idxLevel << std::endl
						<< "[STAT]     >> Nb Cells =                                 \t " << nbCellsAtLevel << std::endl
						<< "[STAT]     >> Nb Adaptive Cells =                   \t" << adaptiveCell    << "  Non Adaptive (1 son): " <<100*FReal(nbCellsAtLevel-adaptiveCell)/nbCellsAtLevel<<std::endl
						<< "[STAT]     >> Number of M2M to remove: " << 		removeM <<	 std::endl
						<< "[STAT]     >> Nb M2M/L2L interactions =        \t" << nbChildAtLevel << std::endl
						<< "[STAT]     >> Average M2M/L2L interactions = \t" << FReal(nbChildAtLevel)/FReal(nbCellsAtLevel) << std::endl
						<< "[STAT]     >> Nb M2L interactions =                 \t" << nbNeighborsAtLevel << std::endl;
				std::cout << "[STAT]     >> M2L Neighbors for each leaf " << std::endl
						<< "[STAT]             >>  Min:        " <<  minM2L << std::endl
						<<  "[STAT]             >> Max:       " <<  maxM2L << std::endl
						<<  "[STAT]             >> Average: " <<  averageM2LNeighbors<< std::endl
						<<  "[STAT]             >> Variance: " <<  varianceM2LNeighbors << std::endl<< std::endl;

				totalCells += (long long int)(nbCellsAtLevel);
				totalM2L += (long long int)(nbNeighborsAtLevel);
				totalM2ML2L += (long long int)(nbChildAtLevel);
				nbCellsAtTop = nbCellsAtLevel;
				if( idxLevel == NbLevels - 1 ) nbCellsAtBottom = nbCellsAtLevel;
				std::cout << std::endl;
				//
				//  Go to next level
				octreeIterator.moveUp();
				octreeIterator.gotoLeft();
				//
			}
			//
			// Global statistics on the octree
			//
			std::cout << "[STAT] For all the tree\n";
			std::cout << "[STAT] >> Total Nb Cells = " << totalCells-nbCellsAtTop << "\n";
			std::cout << "[STAT] >> Total Nb M2M/L2L interactions = " << totalM2ML2L << "\n";
			std::cout << "[STAT] >> Total Average M2M/L2L interactions = " << FReal(totalM2ML2L-nbCellsAtTop)/FReal(totalCells-nbCellsAtBottom) << "\n";
			std::cout << "[STAT] >> Total Nb M2L interactions per cell = " << totalM2L << "\n";
			std::cout << "[STAT] >> Total Average M2L interactions per cell = " << FReal(totalM2L)/FReal(totalCells) << "\n";

			std::cout << "nbCellsAtTop " << nbCellsAtTop <<std::endl;
			//		idCell = totalCells ;
		}

	}
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



