// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Components/FTypedLeaf.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithmTsm.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../../Src/Files/FFmaTsmLoader.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


// Simply create particles and try the kernels
template <class CellClass, class ContainerClass, class LeafClass, class OctreeClass,
          class KernelClass, class FmmClass, typename... Args>
int testFunction(int argc, char ** argv, Args... kernelPreArgs){
    FTic counter;
    // Retrieve parameters
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    // Get working file
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.tsm.fma");
    std::cout << "Opening : " << filename << "\n";
    // Create particles loader
    FFmaTsmLoader loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------
    // Build the tree
    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint particlePosition;
        FReal physicalValue = 0.0;
        FParticleType particleType;
        loader.fillParticle(&particlePosition,&physicalValue,&particleType);
        tree.insert(particlePosition, particleType, physicalValue );
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Create kernel ..." << std::endl;
    counter.tic();

    KernelClass kernels( kernelPreArgs... , NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    counter.tac();
    std::cout << "Done  " << " in " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;

    FmmClass algo(&tree,&kernels);

    counter.tic();
    algo.execute();
    counter.tac();

    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    { // get sum forces&potential
        FReal potential = 0;
        FReal fx = 0.0, fy = 0.0, fz = 0.0;

        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();

            for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                potential += potentials[idxPart];
                fx += forcesX[idxPart];
                fy += forcesY[idxPart];
                fz += forcesZ[idxPart];
            }
        });

        std::cout << "Foces Sum  x = " << fx << " y = " << fy << " z = " << fz << std::endl;
        std::cout << "Potential = " << potential << std::endl;
    }


    // -----------------------------------------------------

    return 0;
}

// This is the real main!
int main(int argc, char ** argv){
    std::cout << "[PARAM] Use Parameters -spherical -rotation -chebyshev\n";

    if( FParameters::existParameter(argc,argv,"-spherical") ){
        std::cout << "[INFO] -spherical is used\n";
        // Create template
        typedef FTypedSphericalCell            CellClass;
        typedef FP2PParticleContainer         ContainerClass;

        typedef FTypedLeaf< ContainerClass >                      LeafClass;
        typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
        typedef FSphericalKernel< CellClass, ContainerClass >          KernelClass;

        typedef FFmmAlgorithmTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        const int DevP = FParameters::getValue(argc,argv,"-p", 8);
        CellClass::Init(DevP);

        // Call Main function
        testFunction< CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass>(argc, argv, DevP);
    }

    if( FParameters::existParameter(argc,argv,"-rotation") ){
        std::cout << "[INFO] -rotation is used\n";
        // Create template
        static const int P = 9;
        typedef FTypedRotationCell<P>            CellClass;
        typedef FP2PParticleContainer         ContainerClass;

        typedef FTypedLeaf< ContainerClass >                      LeafClass;
        typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
        typedef FRotationKernel< CellClass, ContainerClass, P >          KernelClass;

        typedef FFmmAlgorithmTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        // Call Main function
        testFunction< CellClass, ContainerClass, LeafClass, OctreeClass, KernelClass, FmmClass>(argc, argv);
    }

    return 0;
}