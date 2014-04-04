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



#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"


#include "../../Src/Utils/FMath.hpp"

#include "../../Src/Files/FFmaLoader.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"

/// \file  statOnOctree.cpp
//!
//! \brief statOnOctree: Driver to generate statistics on the octree for the particle distribution given bt parameter -infile
//!  \authors
//!
//! The goa

//!
//!  <b> General arguments:</b>
//!     \param    -h     to see the parameters available in this driver
//!     \param  -N     The number of points in the distribution
//!     \param   -infile: Name of the particles file
//!     \param   -outfile Name of the file for output data  (without extension)

//!  <b> Statistics options:</b>
//!   \param -histP build a file to generate histogram of particles per leaf. The data are store in file with .txt extension (-coutfile)

// Simply create particles and try the kernels
int main(int argc, char ** argv){
	typedef FBasicParticleContainer<0>      ContainerClass;
	typedef FSimpleLeaf< ContainerClass >                     LeafClass;
	typedef FOctree< FBasicCell, ContainerClass , LeafClass >  OctreeClass;
	///////////////////////What we do/////////////////////////////
	std::cout << ">> This executable has to be used to show some stat about the tree.\n";
	//////////////////////////////////////////////////////////////

	const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
	const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
	FTic counter;

	const char* const filename = FParameters::getStr(argc,argv,"-infile", "../Data/test20k.fma");
	const std::string genericFileName(FParameters::getStr(argc,argv,"-outfile",   "output"));

	std::cout << "Opening : " << filename << "\n";

	FFmaLoader loader(filename);
	if(!loader.isOpen()){
		std::cout << "Loader Error, " << filename << " is missing\n";
		return 1;
	}

	// -----------------------------------------------------

	OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
	//
	// -----------------------------------------------------
	//     Creating and Inserting particles in the tree
	// -----------------------------------------------------
	//
	std::cout << "Creating and Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
	std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
	counter.tic();

	for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
		FPoint particlePosition;
		FReal physicalValue;
		loader.fillParticle(&particlePosition,&physicalValue);
		tree.insert(particlePosition );
	}
	counter.tac();
	std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;
	//
	// -----------------------------------------------------
	//     Start statistics
	// -----------------------------------------------------
	//
	{ // get stats
		{    // get stats on the leaf level (Particles)
			std::cout << std::endl<< "[STAT] Leaf level "  << " level is  " << NbLevels << std::endl;
			std::cout << "[STAT] potentials leafs number is " << (1 << (3* (NbLevels-1) )) << std::endl;

			FReal averageParticles = 0.0, varianceParticles = 0.0 ;
			int nbLeafs = 0,minParticles = 1000000.0, maxParticles = 0.0 ;
			OctreeClass::Iterator octreeIterator(&tree);
//
			// Start to cmpute statistics on particles
			//
			int  nbPart ,nbTPart=0;
			octreeIterator.gotoBottomLeft();
			do{
				nbPart                   = octreeIterator.getCurrentListTargets()->getNbParticles() ;
				minParticles          =  FMath::Min(minParticles,nbPart) ;
				maxParticles         =  FMath::Max(maxParticles,nbPart) ;
				nbTPart              += nbPart;
				varianceParticles += nbPart*nbPart ;
				//				std::cout <<   nbPart << "  "  << minParticles << "  "<< maxParticles <<std::endl;
				++nbLeafs;
			} while(octreeIterator.moveRight());
			averageParticles   = nbTPart/FReal(nbLeafs);
			varianceParticles /= FReal(nbLeafs);

			//
			std::cout << "[STAT] Nb (non empty) leafs: " << nbLeafs << std::endl;
			std::cout << "[STAT]  particles on leafs:"  << std::endl
					<<   "[STAT]           Average: "<< averageParticles << std::endl
					<<   "[STAT]           Min:        "<< minParticles << std::endl
					<<   "[STAT]           Max:       "<< maxParticles << std::endl
			       <<    "[STAT]           Variance: " << (varianceParticles - (averageParticles*averageParticles)) << std::endl;
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

			FReal averageNeighbors = 0;
			octreeIterator.gotoBottomLeft();
			do{
				ContainerClass* neighbors[27];
				// need the current particles and neighbors particles
				averageNeighbors += FReal(tree.getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),NbLevels-1));
				++nbLeafs;
			} while(octreeIterator.moveRight());
			std::cout << "[STAT] Average neighbors for each leafs = " << (averageNeighbors/FReal(nbLeafs)) << std::endl<< std::endl;
		}
		//    ---------------  END LEAVES ---------------

		{
			long long int totalCells = 0;
			long long int totalM2L = 0;
			long long int totalM2ML2L = 0;

			int nbCellsAtTop = 0;
			int nbCellsAtBottom = 0;

			OctreeClass::Iterator octreeIterator(&tree);
			octreeIterator.gotoBottomLeft();

			for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel){

				int nbCellsAtLevel = 0;
				int nbChildAtLevel = 0;
				int nbNeighborsAtLevel = 0;

				do{
					++nbCellsAtLevel;

					if( idxLevel != NbLevels - 1 ){
						FBasicCell** child = octreeIterator.getCurrentChild();
						for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
							if(child[idxChild]) ++nbChildAtLevel;
						}
					}

					const FBasicCell* neighbors[343];
					nbNeighborsAtLevel += tree.getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),idxLevel);

				} while(octreeIterator.moveRight());

				octreeIterator.moveUp();
				octreeIterator.gotoLeft();

				std::cout << "[STAT] Level = " << idxLevel << "\n";
				std::cout << "[STAT] >> Nb Cells = " << nbCellsAtLevel << "\n";
				std::cout << "[STAT] >> Nb M2M/L2L interactions = " << nbChildAtLevel << "\n";
				std::cout << "[STAT] >> Average M2M/L2L interactions = " << FReal(nbChildAtLevel)/FReal(nbCellsAtLevel) << "\n";
				std::cout << "[STAT] >> Nb M2L interactions = " << nbNeighborsAtLevel << "\n";
				std::cout << "[STAT] >> Average M2L interactions = " << FReal(nbNeighborsAtLevel)/FReal(nbCellsAtLevel) << "\n";

				totalCells += (long long int)(nbCellsAtLevel);
				totalM2L += (long long int)(nbNeighborsAtLevel);
				totalM2ML2L += (long long int)(nbChildAtLevel);
				nbCellsAtTop = nbCellsAtLevel;
				if( idxLevel == NbLevels - 1 ) nbCellsAtBottom = nbCellsAtLevel;
			}

			std::cout << "[STAT] For all the tree\n";
			std::cout << "[STAT] >> Total Nb Cells = " << totalCells-nbCellsAtTop << "\n";
			std::cout << "[STAT] >> Total Nb M2M/L2L interactions = " << totalM2ML2L << "\n";
			std::cout << "[STAT] >> Total Average M2M/L2L interactions = " << FReal(totalM2ML2L-nbCellsAtTop)/FReal(totalCells-nbCellsAtBottom) << "\n";
			std::cout << "[STAT] >> Total Nb M2L interactions per cell = " << totalM2L << "\n";
			std::cout << "[STAT] >> Total Average M2L interactions per cell = " << FReal(totalM2L)/FReal(totalCells) << "\n";

		}
	}

	// -----------------------------------------------------


	return 0;
}



