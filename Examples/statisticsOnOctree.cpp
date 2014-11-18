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

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Components/FBasicCell.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FBasicParticleContainer.hpp"


#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../Src/Utils/FParameterNames.hpp"

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
//!  <b> Statistics options:</b>
//!   \param -histP build a file to generate histogram of particles per leaf. The data are store in file given by -outfile arguments and  .txt extension

// Simply create particles and try the kernels
//



int main(int argc, char ** argv){
    const FParameterNames LocalOptionHist = {
        {"-histP", "--histogram-stat"} ,
         "Build the histogram of the particle number per leaf."
    };
    FHelpDescribeAndExit(argc, argv,
                         "Driver to obtain statistics on the octree.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::OutputFile, LocalOptionHist);


    typedef FBasicParticleContainer<0>   ContainerClass;
    typedef FSimpleLeaf< ContainerClass >  LeafClass;
	typedef FOctree< FBasicCell, ContainerClass , LeafClass >  OctreeClass;

	//
	//   Octree parameters
	//
    const int TreeHeight        = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SubTreeHeight  = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
	//
	//  input and output  Files parameters
	//
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    const std::string genericFileName(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,   "output"));
	//
	std::cout <<	 "Parameters  "<< std::endl
			<<     "      Octree Depth      "<< TreeHeight <<std::endl
			<<	  "      SubOctree depth " << SubTreeHeight <<std::endl
			<<std::endl;

	//
	FFmaGenericLoader loader(filename);
	// -----------------------------------------------------
	OctreeClass tree(TreeHeight, SubTreeHeight,loader.getBoxWidth(),loader.getCenterOfBox());
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
	std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;
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
	// -----------------------------------------------------
	//     Start statistics
	// -----------------------------------------------------
	//
	{ // get stats
		{    // get stats on the leaf level (Particles)
			long int allLeaves =  (1 << (3* (TreeHeight-1) )) ;
			std::cout << std::endl<< "[STAT] Leaf level "  << " is  " << TreeHeight -1<< std::endl;
			std::cout << "[STAT] potentials leafs number is " << allLeaves<< std::endl;

			FReal averageParticles = 0.0, varianceParticles = 0.0 ;
			int nbLeafs = 0,minParticles = 1000000.0, maxParticles = 0.0 ;
			OctreeClass::Iterator octreeIterator(&tree);
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
			//
			//  Histogram of particles per leaf
			//

            if(FParameters::existParameter(argc, argv, LocalOptionHist.options)){
				int size = maxParticles+1;
				int * hist = new int [size] ;
				memset(hist,0,(size)*sizeof(int));
				octreeIterator.gotoBottomLeft();
				do{
                    nbPart  = octreeIterator.getCurrentListSrc()->getNbParticles() ;
					++hist[nbPart] ;
				} while(octreeIterator.moveRight());
				//
				// write data
				//
				std::string name(genericFileName + ".txt");
				std::ofstream outfile( name.c_str(), std::ofstream::out);
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
				nbBox = tree.getLeafsNeighbors(neighborsP2P, octreeIterator.getCurrentGlobalCoordinate(),TreeHeight-1) ;
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

			OctreeClass::Iterator octreeIterator(&tree);
			octreeIterator.gotoBottomLeft();

			for(int idxLevel = TreeHeight - 1 ; idxLevel > 1 ; --idxLevel){

				int nbCellsAtLevel = 0;
				int nbChildAtLevel = 0, adaptiveCell=0 ,nbChildForMyCell;
				int nbNeighborsAtLevel = 0;
				//
				int nbM2LNeighbors, minM2L=500,maxM2L=-1;
				FReal averageM2LNeighbors=0.0, varianceM2LNeighbors=0.0	;
				//
				const FBasicCell* neighborsM2L[343];
				do{
					++nbCellsAtLevel;
					// Check number of
					if( idxLevel != TreeHeight - 1 ){
						nbChildForMyCell=0 ;
						FBasicCell** child = octreeIterator.getCurrentChild();
						for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
							if(child[idxChild]) ++nbChildForMyCell;
						}
						nbChildAtLevel += nbChildForMyCell ;
						if(nbChildForMyCell>1) ++adaptiveCell ;
					}
					const FBasicCell* neighbors[343];
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
				if( idxLevel == TreeHeight - 1 ) nbCellsAtBottom = nbCellsAtLevel;
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

		}
	}

	// -----------------------------------------------------


	return 0;
}



