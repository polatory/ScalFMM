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

// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

// ==== CMAKE =====
// @FUSE_BLAS
// ================

#define BALANCED_PERF

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <string>

#include "ScalFmmConfig.h"

#include "Files/FFmaGenericLoader.hpp"

#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"

#ifdef BALANCED_PERF
#include "BalanceTree/FChebBalanceSymKernel.hpp"
#include "BalanceTree/FCostZones.hpp"
#endif

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"

#ifdef BALANCED_PERF
#include "Core/FFmmAlgorithm.hpp"
#include "BalanceTree/FFmmAlgorithmThreadBalanced.hpp"
#else
#include "Core/FFmmAlgorithmThread.hpp"
#endif

#include "Utils/FParameterNames.hpp"

/**
 * \file
 * \authors B. Bramas, O. Coulaud
 * \brief This program runs the balanced FMM Algorithm with the interpolation
 * kernel based on Chebyshev interpolation (1/r kernel)
 *
 * This program runs the FMM Algorithm with the Chebyshev kernel and compares
 * the results with a direct computation.
 *
 *
 *  This code is a short example to use the Chebyshev Interpolation approach for the 1/r kernel
 */


//
////////////////////////////////////////////////////////////////////
// Chebyshev accuracy
const unsigned int ORDER = 7;
// typedefs
using FReal = double;

#ifdef BALANCED_PERF
using ContainerClass = FP2PParticleContainerIndexed<FReal>;
using LeafClass      = FSimpleLeaf<FReal, ContainerClass >;
using CellClass      = FCostCell<FChebCell<FReal,ORDER>>;
using OctreeClass    = FOctree<FReal,CellClass,ContainerClass,LeafClass>;
//
using MatrixKernelClass  = FInterpMatrixKernelR<FReal>;
using KernelClass        = FChebSymKernel<FReal, CellClass,ContainerClass,MatrixKernelClass,ORDER>;
//
using CostKernelClass = FChebBalanceSymKernel<FReal, CellClass, ContainerClass, 
                                              MatrixKernelClass, ORDER, OctreeClass>;

template <template <typename...> class T, typename Kernel>
using FmmAlgoClass = T<OctreeClass, CellClass, ContainerClass, Kernel, LeafClass>;

using CostFmmClass = FmmAlgoClass<FFmmAlgorithm, CostKernelClass>;
using FmmClass     = FmmAlgoClass<FFmmAlgorithmThreadBalanced, KernelClass>;

const FReal epsilon = 1e-4;

#else

using ContainerClass = FP2PParticleContainerIndexed<FReal>;
using LeafClass      = FSimpleLeaf<FReal, ContainerClass >;
using CellClass      = FChebCell<FReal,ORDER>;
using OctreeClass    = FOctree<FReal,CellClass,ContainerClass,LeafClass>;
using MatrixKernelClass  = FInterpMatrixKernelR<FReal>;
using KernelClass        = FChebSymKernel<FReal, CellClass,ContainerClass,MatrixKernelClass,ORDER>;
using FmmClass = FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass>;

#endif


// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Driver for Chebyshev interpolation kernel  (1/r kernel).",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile,
                         FParameterDefinitions::NbThreads);


    const std::string defaultFile("../Data/unitCubeXYZQ100.bfma" );
    const std::string filename =
        FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str());
    const unsigned int TreeHeight =
        FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight =
        FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const unsigned int NbThreads =
        FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);


    omp_set_num_threads(NbThreads);
    std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;
    //
    std::cout <<     "Parameters  "<< std::endl
              <<     "\tOctree Depth     :" << TreeHeight <<std::endl
              <<     "\tSubOctree depth  :" << SubTreeHeight <<std::endl
              <<     "\tInput file  name :" << filename <<std::endl
              <<     "\tThread number    :" << NbThreads <<std::endl
              <<std::endl;
    ////////////////////////////////////////////////////////////////////

    // init timer
    FTic time;

    // open particle file
    FFmaGenericLoader<FReal> loader(filename);

    // init oct-tree
    OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
    const MatrixKernelClass MatrixKernel;


    { // -----------------------------------------------------
        std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                  << " particles ..."
                  << std::endl
                  << "\tHeight : " << TreeHeight
                  << " \t sub-height : " << SubTreeHeight
                  << std::endl;
        time.tic();
        //
        FPoint<FReal> position;
        FReal physicalValue = 0.0;
        //
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            //
            // Read particle per particle from file
            loader.fillParticle(&position,&physicalValue);
            //
            // put particle in octree
            tree.insert(position, idxPart, physicalValue);
        }

        time.tac();
        std::cout << "Done  " << "(@Creating and Inserting Particles = "
                  << time.elapsed() << " s) ." << std::endl;
    } // -----------------------------------------------------

    { // -----------------------------------------------------
        std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ") ... " << std::endl;

#ifdef BALANCED_PERF
        CostKernelClass balanceKernel(&tree, epsilon);
        CostFmmClass costAlgo(&tree, &balanceKernel);

        time.tic();
        costAlgo.execute();
        time.tac();
        std::cout << "Generating tree cost: " << time.elapsed() << "s.\n";

        FCostZones<OctreeClass, CellClass> costzones(&tree, omp_get_max_threads());

        time.tic();
        costzones.run();
        time.tac();
        std::cout << "Generating cost zones: " << time.elapsed() << "s.\n";
#endif

        time.tic();
        //
        KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
        //
#ifdef BALANCED_PERF
        FmmClass algo(&tree, &kernels, costzones.getZoneBounds(), costzones.getLeafZoneBounds());
#else
        FmmClass algo(&tree, &kernels);
#endif
        //
        algo.execute();   // Here the call of the FMM algorithm
        //
        time.tac();
        std::cout << "Timers Far Field \n"
                  << "P2M " << algo.getTime(FAlgorithmTimers::P2MTimer) << " seconds\n"
                  << "M2M " << algo.getTime(FAlgorithmTimers::M2MTimer) << " seconds\n"
                  << "M2L " << algo.getTime(FAlgorithmTimers::M2LTimer) << " seconds\n"
                  << "L2L " << algo.getTime(FAlgorithmTimers::L2LTimer) << " seconds\n"
                  << "P2P and L2P " << algo.getTime(FAlgorithmTimers::NearTimer) << " seconds\n"
                  << std::endl;


        std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s) ." << std::endl;
    }
    // -----------------------------------------------------
    //
    // Some output
    //
    //
    { // -----------------------------------------------------
        FSize N1 = 0, N2 = loader.getNumberOfParticles()/2, N3 = loader.getNumberOfParticles() - 1;
        FReal energy = 0.0;
        //
        //   Loop over all leaves
        //
        std::cout << std::endl;
        std::cout << std::scientific;
        std::cout.precision(10) ;

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
                    if ((indexPartOrig == N1) || (indexPartOrig == N2) || (indexPartOrig == N3)  ) {
                        std::cout << "Index "<< indexPartOrig <<"  potential  " << potentials[idxPart]
                                  << " Pos "<<posX[idxPart]<<" "<<posY[idxPart]<<" "<<posZ[idxPart]
                                  << "   Forces: " << forcesX[idxPart] << " " << forcesY[idxPart] << " "<< forcesZ[idxPart] <<std::endl;
                    }
                    energy += potentials[idxPart]*physicalValues[idxPart] ;
                }
            });
        std::cout <<std::endl<<"Energy: "<< energy<<std::endl;

    }
    // -----------------------------------------------------


    return 0;
}
