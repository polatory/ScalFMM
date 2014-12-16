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
//
//
//  Utils/Release/FmmAlgorithmTsm -N 2 -h 4
//  plot "output.txt"  using ($1):(log10($3)) w l,"output.txt"  using ($1):(log10($4)) w l,"output.txt"  using ($1):(log10($5)) w l,"output.txt"  using ($1):(log10($6)) w l

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>


#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Utils/FMath.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FPoint.hpp"
#include "Utils/FParameterNames.hpp"
#include "Utils/FTic.hpp"
#include "Utils/FTemplate.hpp"

#include "Components/FTypedLeaf.hpp"


#include "Extensions/FExtendCellType.hpp"

#include "Core/FFmmAlgorithmTsm.hpp"
#include "Core/FFmmAlgorithmThreadTsm.hpp"


#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Files/FRandomLoader.hpp"


struct sourcePart {
	FPoint   position;
	double physicalValue ;
};
/** This program show an example of use of
 * the fmm basic algo
 * it also check that each particles is impacted each other particles
 */
template <int ORDER>
class FTestCellTsm: public FChebCell<ORDER> , public FExtendCellType{
};

// Simply create particles and try the kernels
struct TempMainStruct{
	template <const unsigned int ORDER>
	static  int Run(int argc, char ** argv){
		//	 int main(int argc, char ** argv){
		const FParameterNames LocalOptionPoints = {
				{"-MT"} ,
				"Number of target points to build."
		};
		FHelpDescribeAndExit(argc, argv,
				"Test FMM TSM (target source model) algorithm by counting the nb of interactions each particle receive.",
				FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
				FParameterDefinitions::NbParticles,LocalOptionPoints);

		//	const int ORDER=4 ;
		typedef FP2PParticleContainerIndexed<>                    ContainerClassTyped;
		typedef FTestCellTsm<ORDER>                                   CellClassTyped;

		typedef FTypedLeaf< ContainerClassTyped >                      LeafClassTyped;
		typedef FOctree< CellClassTyped, ContainerClassTyped , LeafClassTyped >  OctreeClassTyped;

		//	typedef FTestKernels< CellClassTyped, ContainerClassTyped >          KernelClassTyped;
		typedef FInterpMatrixKernelR                                                                              MatrixKernelClass;
		const MatrixKernelClass    MatrixKernel;
		typedef FChebSymKernel<CellClassTyped,ContainerClassTyped,MatrixKernelClass,ORDER>  KernelClassTyped;

		typedef FFmmAlgorithmThreadTsm<OctreeClassTyped, CellClassTyped, ContainerClassTyped, KernelClassTyped, LeafClassTyped > FmmClassTyped;
		///////////////////////What we do/////////////////////////////
		std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
		//////////////////////////////////////////////////////////////

		const int TreeHeight          = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
		const int SubTreeHeight    = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
		const int NbPart                = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 20);
		const std::string filename  = FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options,  "output");
		const int nbTargets           = FParameters::getValue(argc,argv,LocalOptionPoints.options,  100);

		FTic counter;

		//////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////
		//
		std::cout <<	 "Parameters  "<< std::endl
				<<     "      Octree Depth      " << TreeHeight <<std::endl
				<<	  "      SubOctree depth " << SubTreeHeight <<std::endl
				<<std::endl;
		//
		//
		//
		double boxSize = 1.0 ;
		FPoint CentreBox(0.5*boxSize,0.5*boxSize,0.5*boxSize) ;

		double cellSize  = 1.0/FMath::pow2(TreeHeight);
		std::cout << " Cell size: "<< cellSize << "  " << 1.0/FMath::pow2(TreeHeight-1)<<  std::endl;

		const int dimGrid = (1 << (TreeHeight-1));
		const FReal dimLeaf = (boxSize/FReal(dimGrid));
		const FReal quarterDimLeaf = (dimLeaf/4.0);
		std::cout << "dimLeaf: "<< dimLeaf << " quarterDimLeaf: "<< quarterDimLeaf << " nbGrid "<< dimGrid<< std::endl;
		//////////////////////////////////////////////////////////////////////////////////
		//
		// Sources are inside the cell number 1
		//
		FPoint CentreBoxSource(0.5*dimLeaf,0.5*dimLeaf,0.5*dimLeaf) ;
		FRandomLoaderTsm loader(NbPart, dimLeaf, CentreBoxSource, 1);
		OctreeClassTyped tree(TreeHeight, SubTreeHeight,boxSize,CentreBox);
		//
		//////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////

		std::cout << "Creating " << NbPart << " particles ..." << std::endl;
		counter.tic();

		FPoint particlePosition;
		double physicalValue = 1.0;
		int nbSRC = loader.getNumberOfParticles();
		sourcePart *tabSrc =  new sourcePart[nbSRC ];

		{
			// Insert sources
			FParticleType particleType, source = FParticleTypeSource;
			for(int idxPart = 0 ; idxPart < nbSRC ; ++idxPart){
				loader.fillParticle(&particlePosition, &particleType);
//				std::cout << idxPart << "  " << particlePosition << "  type " << particleType
//						<<  "  physicalValue: " << physicalValue<< std::endl;
				tree.insert(particlePosition, source, idxPart, physicalValue);
				tabSrc[idxPart].position          = particlePosition ;
				tabSrc[idxPart].physicalValue = physicalValue ;
				physicalValue                        =-physicalValue;
			}
		}
		{
			//
			// Insert target
			//
			physicalValue = 1.0;
			double dx = boxSize/(nbTargets-1) ;
//
			std::cout << "   TARGETS "  <<std::endl;
			FPoint particlePosition(-dx,dimLeaf+quarterDimLeaf,quarterDimLeaf);
			//		int nbTargets = 256;
			for(int idxPart = 0 ; idxPart < nbTargets; ++idxPart){
				particlePosition.incX(dx);
				std::cout << idxPart << "  " <<particlePosition.getX()/dimLeaf<< "   " <<  particlePosition << "  type " << FParticleTypeTarget
						<< "  " <<physicalValue<<std::endl;
				tree.insert(particlePosition, FParticleTypeTarget,idxPart,physicalValue );

			}
		}

		counter.tac();
		std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;


		//////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////

		std::cout << "Working on particles ..." << std::endl;
		counter.tic();

		////	KernelClassTyped kernels;
		KernelClassTyped kernels(TreeHeight, boxSize,CentreBox,&MatrixKernel);
		//
		FmmClassTyped algo(&tree,&kernels);
		algo.execute();

		counter.tac();
		std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;
		std::stringstream ss ; ss <<filename<< "-"<<ORDER <<".txt" ;
		std::ofstream out( ss.str(), std::ofstream::out);
		std::ofstream outMax("outputMax.txt",std::fstream::out | std::fstream::app);


		//////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////
		double errorMax =0.0,  errorNorm2 =0.0 , d2=0.0;
		tree.forEachLeaf([&](LeafClassTyped* leaf){
			const FReal*const posX = leaf->getTargets()->getPositions()[0];
			const FReal*const posY = leaf->getTargets()->getPositions()[1];
			const FReal*const posZ = leaf->getTargets()->getPositions()[2];
			//
			const FReal*const potentials = leaf->getTargets()->getPotentials();
			const FReal*const forcesX = leaf->getTargets()->getForcesX();
			const FReal*const forcesY = leaf->getTargets()->getForcesY();
			const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
			const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
			//
			const FVector<int>& indexes = leaf->getTargets()->getIndexes();
			for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
				double pot = 0.0 , xx,yy,zz;
				FPoint  FF ;
				for(int idxSrc = 0 ; idxSrc < nbSRC ; ++idxSrc){
					xx = posX[idxPart]-tabSrc[idxSrc].position.getX();
					yy  = posY[idxPart]-tabSrc[idxSrc].position.getY();
					zz = posZ[idxPart]-tabSrc[idxSrc].position.getZ();
					FPoint dx(xx,yy,zz) ;
					double dd = 1.0/dx.norm() ;
					double coeff = tabSrc[idxSrc].physicalValue*dd;
					pot  += coeff;
					dx *= coeff*dd*dd ;
					FF -= dx;
				}
				d2           += pot*pot;
				double d = FMath::Abs(potentials[idxPart] - pot);
				errorMax = FMath::Max(errorMax, d);
				errorNorm2 +=  d*d ;
				out << posX[idxPart] << "  " << potentials[idxPart] << "  " << d
				<< "  " << FMath::Abs(forcesX[idxPart] - FF.getX())
				<< "  " << FMath::Abs(forcesY[idxPart] - FF.getY())
				<< "  " << FMath::Abs(forcesZ[idxPart] - FF.getZ())<<std::endl;
			}
		});
		outMax << ORDER << "    "  << errorMax << "    "  << FMath::Sqrt(errorNorm2/d2) << std::endl;
		out.close();
		return 0;
	}
};
int main(int argc, char** argv){
	const unsigned int order = FParameters::getValue(argc, argv, "-order", 5);
	std::cout << "Order given by user is : " << order << "\n";
	FRunIf::Run<unsigned int, 2, 13, 1, TempMainStruct>(order, argc, argv);
	return 0;
}

