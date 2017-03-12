// ===================================================================================
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

// Keep in private GIT
// @SCALFMM_PRIVATE

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
#include "Kernels/P2P/FP2PR.hpp"
#include "Kernels/P2P/FP2PParticleContainer.hpp"

#include "Components/FBasicCell.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Containers/FOctree.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

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
//!

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    //
    ///////////////////////What we do/////////////////////////////
    FHelpDescribeAndExit(argc, argv,
                         "This executable has to be used to compute  interaction either for periodic or non periodic system.\n"
                         "Example -fin filenameIN.{fma or bfma)\n"
                         "Default input file : ../Data/unitCubeXYZQ20k.fma\n"
                         "For the input file, the extension specifies if the file is binary or not.\n",
                         FParameterDefinitions::InputFile);

    //////////////////////////////////////////////////////////////

    const std::string defaultFile(/*SCALFMMDataPath+*/"../Data/unitCubeXYZQ20k.fma");
    const std::string filenameIn(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options,  defaultFile.c_str()));

    typedef double FReal;

    FTic counter;

    // -----------------------------------------------------
    //  LOADER
    //  -----------------------------------------------------
    // ---------------------------------------------------------------------------------
    // Read  particles in the Octree
    // ---------------------------------------------------------------------------------
    std::cout << "Opening : " << filenameIn << "\n";
    //
    FFmaGenericLoader<FReal> loader(filenameIn);
    //
    FSize nbParticles = static_cast<int>(loader.getNumberOfParticles());
    std::cout << "Read " << nbParticles << " particles ..." << std::endl;
    double BoxWith=loader.getBoxWidth();
    FPoint<FReal> Centre(loader.getCenterOfBox().getX(), loader.getCenterOfBox().getY() , loader.getCenterOfBox().getZ());
    std::cout << "\tWidth : " <<BoxWith << " \t center x : " << loader.getCenterOfBox().getX()
              << " y : " << loader.getCenterOfBox().getY() << " z : " << loader.getCenterOfBox().getZ() << std::endl;

    counter.tic();
    FmaRWParticle<FReal, 4,8> *  particles = new FmaRWParticle<FReal, 4,8>[nbParticles];
    memset(particles, 0, sizeof(FmaRWParticle<FReal, 4,8>) * nbParticles) ;
    //
    double totalCharge = 0.0;
    //
    //	int nbDataToRead = particles[0].getReadDataNumber();
    for(int idx = 0 ; idx<nbParticles ; ++idx){
        //
        loader.fillParticle(particles[idx].getPtrFirstData(), particles[idx].getReadDataNumber());
        totalCharge += particles[idx].getPhysicalValue() ;
    }

    counter.tac();

    std::cout << std::endl;
    std::cout << "Total Charge         = "<< totalCharge <<std::endl;
    std::cout << std::endl;

    std::cout << "Done  " << "(@ reading Particles  " << counter.elapsed() << " s)." << std::endl;
    //Need to copy particles to ContainerClass -> FP2PParticleContainer
    typedef FBasicCell                 CellClass;
    typedef FP2PParticleContainer<FReal>         ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;

    OctreeClass tree(2, 1, loader.getBoxWidth(), loader.getCenterOfBox());
    for(int idxP=0 ; idxP<nbParticles ; ++idxP){
        tree.insert(particles[idxP].getPosition(),particles[idxP].getPhysicalValue());
    }
    int k=0;
    // tree.forEachLeaf([&](LeafClass * leaf){
    //	printf("leaf : %d\n",k++ );
    //     });
    //
    // ----------------------------------------------------------------------------------------------------------
    //                                   COMPUTATION
    // ----------------------------------------------------------------------------------------------------------
    //
    //  computation
    //
    printf("Precision, sizeof Real %lu\n", sizeof(FReal));

#ifdef SCALFMM_USE_AVX
    printf("AVX incomming .......\n\n");
#endif

#ifdef SCALFMM_USE_SSE
    printf("SSE incomming .......\n\n");
#endif

#ifndef SCALFMM_USE_SSE
#ifndef SCALFMM_USE_AVX
    printf("Classic incomming ...\n\n");
#endif
#endif
    counter.tic();
    {
        typename OctreeClass::Iterator iterator(&tree);
        iterator.gotoBottomLeft();

        do{
            FTreeCoordinate coord = iterator.getCurrentGlobalCoordinate();
            ContainerClass** neighbors = new ContainerClass*[27];
            tree.getLeafsNeighbors(neighbors,coord,1);
            FP2PRT<FReal>::FullMutual<ContainerClass>(iterator.getCurrentLeaf()->getTargets(),neighbors,27);

        }while(iterator.moveRight());
    }
    counter.tac();
    std::cout << "Done  " << "(@ Computation  " << counter.elapsed() << " s)." << std::endl;
    FReal cumulPot = 0.0;
    k=0;
    tree.forEachLeaf([&](LeafClass * leaf){
        FSize maxParts = leaf->getSrc()->getNbParticles();
        FReal* datas = leaf->getSrc()->getPotentials();
        for(int i=0 ; i<maxParts ; ++i){
            cumulPot += datas[i];
        }
        printf("leaf : %d --> cumul pot %e : \n",k++, cumulPot);
    });


    delete[] particles;
    return 0;
}
