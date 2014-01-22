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
#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/Rotation/FRotationCell.hpp"
#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationOriginalKernel.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Core/FFmmAlgorithmPeriodic.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMath.hpp"


int main(int argc, char** argv){
    /////////////////////////////////////////////////
    // Types
    /////////////////////////////////////////////////
    static const int P = 6;
    typedef FRotationCell<P>            CellClass;
    typedef FP2PParticleContainerIndexed  ContainerClass;
    typedef FRotationKernel< CellClass, ContainerClass, P >   KernelClass;
    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree<CellClass,ContainerClass,LeafClass> OctreeClass;
    typedef FFmmAlgorithmPeriodic<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassPer;
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

    /////////////////////////////////////////////////
    // Parameters
    /////////////////////////////////////////////////

    // Parameters
    const int NbLevels       = 4;
    const int SizeSubLevels  = 2;
    const int MinLevelAbove  = FParameters::getValue(argc, argv, "-min",-1);
    const int MaxLevelAbove  = FParameters::getValue(argc, argv, "-max",3);
    const int IncLevelAbove  = 1;
    const int NbParticles    = FParameters::getValue(argc, argv, "-nb",6);
    const int NbParticlesPrint    = FParameters::getValue(argc, argv, "-nbprint", FMath::Min(6, NbParticles));
    FAssertLF(NbParticlesPrint <= NbParticles , "The number of printer particles cannot be higer than the number of particles.");

    std::cout << "The application will use " << NbParticles << " but studies only " << NbParticlesPrint << " of them." << std::endl;

    /////////////////////////////////////////////////
    // Insert particlePositions in tree and copy into array
    /////////////////////////////////////////////////

    FRandomLoader loader(NbParticles);
    FPoint* const particlePositions = new FPoint[NbParticles];

    for(int idxPart = 0 ; idxPart < NbParticles ; ++idxPart){
        loader.fillParticle(&particlePositions[idxPart]);
    }

    FReal* allPotential = new FReal[(MaxLevelAbove-MinLevelAbove)*NbParticlesPrint];

    /////////////////////////////////////////////////
    // Test for different periodicity for FMM!!
    /////////////////////////////////////////////////
    for(int idxLevelAbove = MinLevelAbove; idxLevelAbove < MaxLevelAbove ; idxLevelAbove += IncLevelAbove){
        OctreeClass treePer(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

        // insert in tree
        for(int idxPart = 0 ; idxPart < NbParticles ; ++idxPart){
            // put in tree
            treePer.insert(particlePositions[idxPart], idxPart, (idxPart&1?0.00010:-0.00010));
        }

        // Run FMM
        FmmClassPer algoPer(&treePer,idxLevelAbove);
        KernelClass kernelsPer( algoPer.extendedTreeHeight(), algoPer.extendedBoxWidth(), algoPer.extendedBoxCenter());
        algoPer.setKernel(&kernelsPer);
        algoPer.execute();

        FTreeCoordinate min, max;
        algoPer.repetitionsIntervals(&min,&max);
        std::cout << "Nb Levels above root = " << idxLevelAbove << "\n";
        std::cout << "\t repetitions       = " << algoPer.theoricalRepetition() << std::endl;
        std::cout << "\t properties was height " << treePer.getHeight() << ", width " << loader.getBoxWidth()
                  << ", center " << loader.getCenterOfBox() << "\n";
        std::cout << "\t new properties height " << algoPer.extendedTreeHeight() << ", width " << algoPer.extendedBoxWidth()
                  << ", center " << algoPer.extendedBoxCenter() << "\n";
        std::cout << "min = " << min << "\n";
        std::cout << "max = " << max << "\n";
        std::cout << std::endl;

        treePer.forEachLeaf([&](LeafClass* leaf){
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FVector<int>& indexes =  leaf->getTargets()->getIndexes();
            for(int idxPart = 0 ; idxPart < leaf->getTargets()->getNbParticles() ; ++idxPart){
                if( indexes[idxPart] < NbParticlesPrint){
                    allPotential[(idxLevelAbove - MinLevelAbove)*NbParticlesPrint + indexes[idxPart]] = potentials[idxPart];
                }
            }
        });
    }

    /////////////////////////////////////////////////
    // Print results
    /////////////////////////////////////////////////


    std::cout << "Part" << "\t";
    for(int idxLevelAbove = MinLevelAbove; idxLevelAbove < MaxLevelAbove ; idxLevelAbove += IncLevelAbove){
        std::cout << idxLevelAbove << "\t";
    }
    std::cout << "\n";
    for(int idxPart = 0 ; idxPart < NbParticlesPrint ; ++idxPart){
        std::cout << idxPart << "\t";
        for(int idxLevelAbove = MinLevelAbove; idxLevelAbove < MaxLevelAbove ; idxLevelAbove += IncLevelAbove){
            std::cout << allPotential[(idxLevelAbove - MinLevelAbove)*NbParticlesPrint + idxPart] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "Part" << "\t";
    for(int idxLevelAbove = MinLevelAbove+1; idxLevelAbove < MaxLevelAbove ; idxLevelAbove += IncLevelAbove){
        std::cout << idxLevelAbove << "/" << idxLevelAbove-1 << "\t";
    }
    std::cout << "\n";
    for(int idxPart = 0 ; idxPart < NbParticlesPrint ; ++idxPart){
        std::cout << idxPart << "\t";
        for(int idxLevelAbove = MinLevelAbove +1; idxLevelAbove < MaxLevelAbove ; idxLevelAbove += IncLevelAbove){
            std::cout << FMath::Abs((allPotential[(idxLevelAbove - MinLevelAbove)*NbParticlesPrint + idxPart]
                         / allPotential[(idxLevelAbove - MinLevelAbove - 1)*NbParticlesPrint + idxPart])-1.0) << "\t";
        }
        std::cout << "\n";
    }

    delete[] allPotential;
    delete[] particlePositions;

    return 0;
}



