// ===================================================================================
// Copyright ScalFmm 2016 INRIA, Olivier Coulaud, Bérenger Bramas,
// Matthias Messner olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the
// FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
// An extension to the license is given to allow static linking of scalfmm
// inside a proprietary application (no matter its license).
// See the main license file for more details.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
//
// ==== CMAKE =====
//
// ================

// ==== CMAKE =====
// @FUSE_BLAS
// ================
// Keep in private GIT
// @SCALFMM_PRIVATE


#include <iostream>
#include <stdexcept>
#include <string>
#include <cstdlib>
#include <cstdio>
//
#include "ScalFmmConfig.h"
#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

//#include "Core/FFmmAlgorithmNoP2P.hpp"
#include "Core/FFmmAlgorithm.hpp"

#ifdef SCALFMM_USE_BLAS
// chebyshev kernel

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
//#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#endif
//
#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#ifdef SCALFMM_USE_FFT
// Uniform grid kernel
#include "Kernels/Uniform/FUnifCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "Kernels/Uniform/FUnifKernel.hpp"
#endif

#include "../../Src/Utils/FParameterNames.hpp"

/**
 * This program compute the potential and the energy in non cubic geometry
 */
template <class classTree>
void compare(std::string& val, const classTree& tree ) {

}
void usage() {
	std::cout << "Driver for testing different approximations  for the  1/r kernel" << std::endl;
	std::cout <<	 "Options  "<< std::endl
			<<     "      -help         to see the parameters    " << std::endl
			<<	  "      -depth       the depth of the octree   "<< std::endl
			<<	  "      -subdepth  specifies the size of the sub octree   " << std::endl
			<<     "      -f   name    name specifies the name of the particle distribution" << std::endl
			<<     "      -t  n  specifies the number of threads used in the computations" << std::endl;
}

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    typedef double FReal;

    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OctreeHeight.options)||FParameters::existParameter(argc, argv, "-help")){
		usage() ;
		std::cout << "Driver for testing different approximations  for the  1/r kernel" << std::endl;

		exit(EXIT_SUCCESS);
	}

	// get info from commande line
    const std::string  filename(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/UTest/unitCubeRef20kDouble.bfma"));
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);

			std::cout <<	 "Parameters  "<< std::endl
			<<     "      Octree Depth      \t"<< TreeHeight <<std::endl
			<<	  "      SubOctree depth \t"<< SubTreeHeight <<std::endl
			<<     "      Input file  name: \t" <<filename <<std::endl
			<<std::endl;

	// init timer
	FTic time;

	FFmaGenericLoader<FReal> loader(filename);
	//

	FSize nbParticles = loader.getNumberOfParticles() ;
	FmaRWParticle<FReal, 8,8>* const particles = new FmaRWParticle<FReal, 8,8>[nbParticles];
	//
	loader.fillParticle(particles,nbParticles);
	FReal LL = loader.getBoxWidth() ;
    FPoint<FReal> maxPos(-LL,-LL,-LL),minPos(LL,LL,LL), BoxWidth;
	for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
        const FPoint<FReal> PP(particles[idxPart].getPosition() ) ;
		//
		minPos.setX(FMath::Min(minPos.getX(),PP.getX())) ;
		minPos.setY(FMath::Min(minPos.getY(),PP.getY())) ;
		minPos.setZ(FMath::Min(minPos.getZ(),PP.getZ())) ;
		maxPos.setX(FMath::Max(maxPos.getX(),PP.getX())) ;
		maxPos.setY(FMath::Max(maxPos.getY(),PP.getY())) ;
		maxPos.setZ(FMath::Max(maxPos.getZ(),PP.getZ())) ;
		//
	}
	BoxWidth = maxPos-minPos;
	BoxWidth.setX(ceil(BoxWidth.getX()));
	BoxWidth.setY(ceil(BoxWidth.getY()));
	BoxWidth.setZ(ceil(BoxWidth.getZ()));

//	BoxWidth = FPoint<FReal>(loader.getBoxWidth(),loader.getBoxWidth(),loader.getBoxWidth() );

	FReal LX = ceil(BoxWidth.getX()),  LY = ceil(BoxWidth.getY()),  LZ = ceil(BoxWidth.getZ()) ;
	std::cout << "Data are inside the box delimited by "<<std::endl
			<< "         Min corner:  "<< minPos<<std::endl
			<< "         Max corner:  "<< maxPos<<std::endl <<std::endl
			<< "         Box size        "  << BoxWidth <<"  " <<  LX
			<<"  " <<  LY  <<"  " <<  LZ << std::endl;
	//
	////////////////////////////////////////////////////////////////////
	//  Compute direct energy
	FReal energyD =0.0, totPhysicalValue =0.0;

	//#pragma omp parallel for reduction(+:energyD,totPhysicalValue)
    for(FSize idx = 0 ; idx < loader.getNumberOfParticles()  ; ++idx){
		energyD             +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
		totPhysicalValue += particles[idx].getPhysicalValue() ;
	}
	std::cout << " Total Physical values: "<< totPhysicalValue <<std::endl;
	std::cout << " Energy of the system: "<< energyD <<std::endl;
	////////////////////////////////////////////////////////////////////
	// Rescale the data into the unit cube
	//
    FPoint<FReal> CenterOfBox(0.5,0.5,0.5) ; //(loader.getCenterOfBox()) ;
	FReal   boxWidth = 1.0 ;// BoxWidth.getX();
	//
	// Scale the Points
	for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
		// put in tree
		//	std::cout << idxPart<< " "<< particles[idxPart].getPosition()<< "  ";
        FPoint<FReal> PP(particles[idxPart].getPosition()/= BoxWidth) ;
		particles[idxPart].setPosition(PP);
		//	std::cout <<particles[idxPart].getPosition()<< std::endl;
	}
    typedef FInterpMatrixKernelRH<FReal> MatrixKernelClass;
	MatrixKernelClass MatrixKernel;
	MatrixKernel.setCoeff(LX,LY,LZ);
	//
    FPoint<FReal> A(1.0,0.0,0.0), B(0.0,1.0,0.0);
	//
	//
    FPoint<FReal> xy(A-B);
	FReal rr = FReal(1.) / FMath::Sqrt(xy.getX()*xy.getX() +
			xy.getY()*xy.getY() + xy.getZ()*xy.getZ());
	//
	A /= BoxWidth; 			B /= BoxWidth;
	std::cout << "Good eval: " << rr  <<  "     "  <<  MatrixKernel.evaluate(A,B) <<std::endl;
	std::cout << " Coeff " <<  LX 	<<"  " <<  LY  <<"  " <<  LZ << std::endl;
	//
#ifdef  SCALFMM_USE_BLAS
	{	// begin ChebSymKernel kernel

		// accuracy
		const unsigned int ORDER = 7;
		std::cout << "\nFChebSymKernel FMM ... ORDER: " << ORDER <<std::endl;

		// typedefs
        typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
        typedef FChebCell<FReal,ORDER> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;

    //	typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FChebKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
		typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

		//
		// init oct-tree on rescaled data
		//
		OctreeClass tree(TreeHeight, SubTreeHeight, boxWidth /*loader.getBoxWidth()*/, CenterOfBox);

		{ // -----------------------------------------------------
			time.tic();
			// Insert and scale the Points
			for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
				// put in tree
				tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
			}

			time.tac();
			std::cout <<  "(FChebSymKernel @Inserting Particles = "<< time.elapsed() << " s)." << std::endl;
		} // -----------------------------------------------------

		{ // -----------------------------------------------------
			time.tic();
			//		KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
			KernelClass kernels(TreeHeight, boxWidth, CenterOfBox,&MatrixKernel);
			//
			FmmClass algorithm(&tree, &kernels);
			algorithm.execute();
			time.tac();
			std::cout << "(FChebSymKernel @Algorithm = " << time.elapsed() << " s)." << std::endl;
		} // -----------------------------------------------------
		FReal energy = 0.0;
		FMath::FAccurater<FReal> potentialDiff;
		FMath::FAccurater<FReal> fx, fy, fz;
		{ // Check that each particle has been summed with all other

			tree.forEachLeaf([&](LeafClass* leaf){
				const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
				const FReal*const potentials        = leaf->getTargets()->getPotentials();
				FReal*const forcesX            = leaf->getTargets()->getForcesX();
				FReal*const forcesY            = leaf->getTargets()->getForcesY();
				FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
				const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

				for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
					const FSize indexPartOrig = indexes[idxPart];
					// Rescale the forces
				 	forcesX[idxPart] /= LX;
					forcesY[idxPart] /= LY;
					forcesZ[idxPart] /= LZ;
					//
					potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
					fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
					fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
					fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
					energy += potentials[idxPart]*physicalValues[idxPart] ;
				}
			});
		}

		// Print for information
		std::cout << "FChebSymKernel Energy "  << FMath::Abs(energy-energyD) /energyD << "   " << energy/energyD<< std::endl;
		std::cout << "FChebSymKernel Potential " << potentialDiff << std::endl;
		std::cout << "FChebSymKernel Fx " << fx << std::endl;
		std::cout << "FChebSymKernel Fy " << fy << std::endl;
		std::cout << "FChebSymKernel Fz " << fz << std::endl;

	} // end Chebyshev kernel
#endif

#ifdef  SCALFMM_USE_FFT
	//
	////////////////////////////////////////////////////////////////////
	//
	{	// begin Lagrange/Uniform Grid kernel

		// TODO

		// accuracy
		const unsigned int ORDER = 10;
		std::cout << "\nLagrange FMM ... ORDER " << ORDER <<std::endl;

		// typedefs

        typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass >  LeafClass;
        //		typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        typedef FUnifCell<FReal,ORDER> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FUnifKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
		typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
		//		typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;


		// init oct-tree
		OctreeClass tree(TreeHeight, SubTreeHeight, boxWidth, CenterOfBox);

		{ // -----------------------------------------------------
			time.tic();

			for(FSize idxPart = 0 ; idxPart <nbParticles ; ++idxPart){
				// put in tree
				tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
			}

			time.tac();
			std::cout << "(Lagrange @Inserting Particles = " << time.elapsed() << " s)." << std::endl;
		} // -----------------------------------------------------

		{ // -----------------------------------------------------
			time.tic();
			KernelClass kernels(TreeHeight, boxWidth, CenterOfBox,&MatrixKernel);
			FmmClass algorithm(&tree, &kernels);
			algorithm.execute();
			time.tac();
			std::cout <<  "(Lagrange @Algorithm = " << time.elapsed() << " s)." << std::endl;
		} // -----------------------------------------------------

		FReal energy = 0.0;
		FMath::FAccurater<FReal> potentialDiff;
		FMath::FAccurater<FReal> fx, fy, fz;
		{ // Check that each particle has been summed with all other

			tree.forEachLeaf([&](LeafClass* leaf){
				const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
				const FReal*const potentials        = leaf->getTargets()->getPotentials();
				FReal*const forcesX            = leaf->getTargets()->getForcesX();
				FReal*const forcesY            = leaf->getTargets()->getForcesY();
				FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
				const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

				for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
					const FSize indexPartOrig = indexes[idxPart];
					forcesX[idxPart] /= LX;
					forcesY[idxPart] /= LY;
					forcesZ[idxPart] /= LZ;
					potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
					fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
					fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
					fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
					energy += potentials[idxPart]*physicalValues[idxPart] ;
				}
			});
		}

		// Print for information
		std::cout << "Lagrange Energy "  << FMath::Abs(energy-energyD) /energyD << std::endl;
		std::cout << "Lagrange Potential " << potentialDiff << std::endl;
		std::cout << "Lagrange Fx " << fx << std::endl;
		std::cout << "Lagrange Fy " << fy << std::endl;
		std::cout << "Lagrange Fz " << fz << std::endl;

	} // end Lagrange/Uniform Grid kernel
#endif
	delete[] particles;

	return 0;

}
