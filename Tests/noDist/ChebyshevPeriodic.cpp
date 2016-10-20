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

// ==== CMAKE =====
// @FUSE_BLAS
// ================
// Keep in private GIT
// @SCALFMM_PRIVATE

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <string>

#include "ScalFmmConfig.h"
#include "Utils/FGlobal.hpp"

#include "Files/FFmaGenericLoader.hpp"

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"

#include "Core/FFmmAlgorithmPeriodic.hpp"


#include "../../Src/Utils/FParameterNames.hpp"

template <class Output>
void Print(const Output& value){
    std::cout<< "--- Output from program : " << value << "\n";
}

/**
 * This program runs the FMM Algorithm with the Chebyshev kernel and compares the results with a direct computation.
 */
/// \file  ChebyshevInterpolationFMM.cpp
//!
//! \brief This program runs the FMM Algorithm with the interpolation kernel based on Chebyshev interpolation (1/r kernel)
//!  \authors B. Bramas, O. Coulaud
//!
//!  This code is a short example to use the Chebyshev Interpolation approach for the 1/r kernel
//!


// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Driver for Chebyshev interpolation kernel  (1/r kernel) with periodicity.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OutputFile,
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::NbThreads, FParameterDefinitions::PeriodicityNbLevels);

    typedef double FReal;

    const std::string defaultFile(/*SCALFMMDataPath+*/"../Data/test20k.fma" );
    const std::string filename                = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str());
    const std::string filenameOut          = FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options, "resultPer.fma");
    const unsigned int TreeHeight        = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight  = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const unsigned int NbThreads        = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);
    const int PeriodicDeep                     = FParameters::getValue(argc,argv,FParameterDefinitions::PeriodicityNbLevels.options, 3);

    //
    std::cout <<	 "Parameters  "<< std::endl
                  <<     "\t      Octree Depth      \t"<< TreeHeight <<std::endl
                  <<	    "\t      SubOctree depth \t" << SubTreeHeight <<std::endl
                  <<	    "\t      Periodic depth    \t" << PeriodicDeep <<std::endl
                  <<     "\t      Input file  name: \t" <<filename <<std::endl
                   <<    "\t      Thread number: \t " << NbThreads <<std::endl
                                   <<std::endl;
    //
    // init timer
    FTic time;

    // open particle file
    ////////////////////////////////////////////////////////////////////
    //
    FFmaGenericLoader<FReal> loader(filename);
    FSize nbParticles = loader.getNumberOfParticles() ;
    FmaRWParticle<FReal,8,8>* const particles = new FmaRWParticle<FReal,8,8>[nbParticles];

    //
    ////////////////////////////////////////////////////////////////////
    // begin Chebyshev kernel

    // accuracy
    const unsigned int ORDER = 7;
    // typedefs
    typedef FP2PParticleContainerIndexed<FReal>                     ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                        LeafClass;
    typedef FChebCell<FReal,ORDER>                                         CellClass;
    typedef FOctree<FReal,CellClass,ContainerClass,LeafClass>  OctreeClass;
    //
    typedef FInterpMatrixKernelR<FReal>                                       MatrixKernelClass;
    typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER>  KernelClass;
    //
    typedef FFmmAlgorithmPeriodic<FReal, OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

    // init oct-tree
    OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


    { // -----------------------------------------------------
        std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                  << " particles ..." << std::endl;
        time.tic();
        //
        FPoint<FReal> position;
        //
        loader.fillParticle(particles,nbParticles);

        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue() );
        }

        time.tac();
        std::cout << "Done  " << "(@Creating and Inserting Particles = "
                  << time.elapsed() << " s) ." << std::endl;
    } // -----------------------------------------------------

    /////////////////////////////////////////////////////////////////////////////////////////////////
    { // -----------------------------------------------------
        std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ") ... " << std::endl;

        time.tic();

        FmmClass algo(&tree,PeriodicDeep );
        const MatrixKernelClass MatrixKernel;
        KernelClass kernels(algo.extendedTreeHeight(), algo.extendedBoxWidth(), algo.extendedBoxCenter(),&MatrixKernel);
        algo.setKernel(&kernels);
        algo.execute();
        //
        time.tac();
        std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s) ." << std::endl;

    }
    // -----------------------------------------------------
    //
    // Some output
    //
    //
    FmaRWParticle<FReal, 8,8>* const particlesOut = new FmaRWParticle<FReal, 8,8>[nbParticles];

    { // -----------------------------------------------------
        //
        FReal energy= 0.0 , energyD = 0.0 ;
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Compute direct energy
        /////////////////////////////////////////////////////////////////////////////////////////////////

        for(FSize idx = 0 ; idx < loader.getNumberOfParticles()  ; ++idx){
            energyD +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
        }
        //
        //   Loop over all leaves
        //
        std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl;
        std::cout << std::scientific;
        std::cout.precision(10) ;
        /////////////////////////////////////////////////////////////////////////////////////////////////
        // Compare
        /////////////////////////////////////////////////////////////////////////////////////////////////
        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;

        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const posX = leaf->getTargets()->getPositions()[0];
            const FReal*const posY = leaf->getTargets()->getPositions()[1];
            const FReal*const posZ = leaf->getTargets()->getPositions()[2];
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
            const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();

            const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

            for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                const FSize indexPartOrig = indexes[idxPart];
                //
                particlesOut[indexPartOrig].setPosition(posX[idxPart],posY[idxPart],posZ[idxPart]) ;
                particlesOut[indexPartOrig].setPhysicalValue(physicalValues[idxPart]) ;
                particlesOut[indexPartOrig].setPotential (potentials[idxPart]) ;
                particlesOut[indexPartOrig].setForces(forcesX[idxPart],forcesY[idxPart],forcesZ[idxPart]) ;

                potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                energy+=potentials[idxPart]*physicalValues[idxPart];
            }
        });
        energy *= 0.5;
        std::cout <<std::endl<<"Energy: "<< energy<<std::endl;
        std::cout <<std::endl<<" &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& "<<std::endl<<std::endl;


        // remove index
        Print("Test1 - Error Relative L2 norm Potential ");
        printf("         Pot L2Norm     %e\n",potentialDiff.getL2Norm());
        printf("         Pot RL2Norm   %e\n",potentialDiff.getRelativeL2Norm());
        printf("         Pot RMSError   %e\n",potentialDiff.getRMSError());
        Print("Fx diff is = ");
        printf("         Fx L2Norm     %e\n",fx.getL2Norm());
        printf("         Fx RL2Norm   %e\n",fx.getRelativeL2Norm());
        printf("         Fx RMSError   %e\n",fx.getRMSError());
        Print("Fy diff is = ");
        printf("        Fy L2Norm     %e\n",fy.getL2Norm());
        printf("        Fy RL2Norm   %e\n",fy.getRelativeL2Norm());
        printf("        Fy RMSError   %e\n",fy.getRMSError());
        Print("Fz diff is = ");
        printf("        Fz L2Norm     %e\n",fz.getL2Norm());
        printf("        Fz RL2Norm   %e\n",fz.getRelativeL2Norm());
        printf("        Fz RMSError   %e\n",fz.getRMSError());
        FReal L2error = (fx.getRelativeL2Norm()*fx.getRelativeL2Norm() + fy.getRelativeL2Norm()*fy.getRelativeL2Norm()  + fz.getRelativeL2Norm() *fz.getRelativeL2Norm()  );
        printf(" Total L2 Force Error= %e\n",FMath::Sqrt(L2error)) ;
        printf("  Energy Error  =   %.12e\n",FMath::Abs(energy-energyD));
        printf("  Energy FMM    =   %.12e\n",FMath::Abs(energy));
        printf("  Energy DIRECT =   %.12e\n",FMath::Abs(energyD));

    }
    // -----------------------------------------------------
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::OutputFile.options) ){

        std::cout << "Generate " << filenameOut <<"  for output file" << std::endl;
        //
        std::cout << " numberofParticles: " << nbParticles <<"  " << sizeof(nbParticles) <<std::endl;
        std::cout << " Box size: " << loader.getBoxWidth() << "  " << sizeof(loader.getBoxWidth())<<std::endl;
        //
        FFmaGenericWriter<FReal> writer(filenameOut) ;
        writer.writeHeader(loader.getCenterOfBox(), loader.getBoxWidth() , nbParticles,*particlesOut) ;
        writer.writeArrayOfParticles(particlesOut, nbParticles);

    }
    delete [] particlesOut;
    return 0;
}
