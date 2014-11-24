#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"
#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/GroupTree/FGroupSeqAlgorithm.hpp"
#include "../../Src/GroupTree/FGroupTaskAlgorithm.hpp"
#include "../../Src/GroupTree/FGroupTaskDepAlgorithm.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../Src/GroupTree/FP2PGroupParticleContainer.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

#include <memory>




int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize { {"-bs"}, "The size of the block of the blocked tree"};
    const FParameterNames LocalOptionNoValidate { {"-no-validation"}, "To avoid comparing with direct computation"};
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,FParameterDefinitions::InputFile,
                          FParameterDefinitions::NbThreads,
                         FParameterDefinitions::NbParticles, LocalOptionBlocSize, LocalOptionNoValidate);

    // Initialize the types
    static const int P = 9;
    typedef FRotationCell<P>               GroupCellClass;
    typedef FP2PGroupParticleContainer<>          GroupContainerClass;
    typedef FGroupTree< GroupCellClass, GroupContainerClass, 5, FReal>  GroupOctreeClass;
    typedef FRotationKernel< GroupCellClass, GroupContainerClass , P>   GroupKernelClass;
#ifdef ScalFMM_USE_STARPU
    typedef FGroupStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#elif defined(ScalFMM_USE_OMP4)
    // Set the number of threads
    omp_set_num_threads(FParameters::getValue(argc,argv,FParameterDefinitions::NbThreads.options, omp_get_max_threads()));
    typedef FGroupTaskDepAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#else
    //typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
    typedef FGroupTaskAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#endif


    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int groupSize     = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");

    // Load the particles
    //FRandomLoader loader(FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 20), 1.0, FPoint(0,0,0), 0);
    FFmaGenericLoader loader(filename);
    FAssertLF(loader.isOpen());
    FTic timer;

    FP2PParticleContainer<> allParticles;
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FReal physicalValue;
        FPoint particlePosition;
        loader.fillParticle(&particlePosition, &physicalValue);
        allParticles.push(particlePosition, physicalValue);
    }
    std::cout << "Particles loaded in " << timer.tacAndElapsed() << "s\n";

    // Put the data into the tree
    timer.tic();
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    groupedTree.printInfoBlocks();
    std::cout << "Tree created in " << timer.tacAndElapsed() << "s\n";

    // Run the algorithm
    GroupKernelClass groupkernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
    GroupAlgorithm groupalgo(&groupedTree,&groupkernel);

    timer.tic();
    groupalgo.execute();
    std::cout << "Kernel executed in in " << timer.tacAndElapsed() << "s\n";

    // Validate the result
    if(FParameters::existParameter(argc, argv, LocalOptionNoValidate.options) == false){
        int offsetParticles = 0;
        FReal*const allPhysicalValues = allParticles.getPhysicalValues();
        FReal*const allPosX = const_cast<FReal*>( allParticles.getPositions()[0]);
        FReal*const allPosY = const_cast<FReal*>( allParticles.getPositions()[1]);
        FReal*const allPosZ = const_cast<FReal*>( allParticles.getPositions()[2]);

        groupedTree.forEachCellLeaf<FP2PGroupParticleContainer<> >([&](GroupCellClass* cellTarget, FP2PGroupParticleContainer<> * leafTarget){
            const FReal*const physicalValues = leafTarget->getPhysicalValues();
            const FReal*const posX = leafTarget->getPositions()[0];
            const FReal*const posY = leafTarget->getPositions()[1];
            const FReal*const posZ = leafTarget->getPositions()[2];
            const int nbPartsInLeafTarget = leafTarget->getNbParticles();

            for(int idxPart = 0 ; idxPart < nbPartsInLeafTarget ; ++idxPart){
                allPhysicalValues[offsetParticles + idxPart] = physicalValues[idxPart];
                allPosX[offsetParticles + idxPart] = posX[idxPart];
                allPosY[offsetParticles + idxPart] = posY[idxPart];
                allPosZ[offsetParticles + idxPart] = posZ[idxPart];
            }

            offsetParticles += nbPartsInLeafTarget;
        });

        FAssertLF(offsetParticles == loader.getNumberOfParticles());

        FReal*const allDirectPotentials = allParticles.getPotentials();
        FReal*const allDirectforcesX = allParticles.getForcesX();
        FReal*const allDirectforcesY = allParticles.getForcesY();
        FReal*const allDirectforcesZ = allParticles.getForcesZ();

        for(int idxTgt = 0 ; idxTgt < offsetParticles ; ++idxTgt){
            for(int idxMutual = idxTgt + 1 ; idxMutual < offsetParticles ; ++idxMutual){
                FP2PR::MutualParticles(
                    allPosX[idxTgt],allPosY[idxTgt],allPosZ[idxTgt], allPhysicalValues[idxTgt],
                    &allDirectforcesX[idxTgt], &allDirectforcesY[idxTgt], &allDirectforcesZ[idxTgt], &allDirectPotentials[idxTgt],
                    allPosX[idxMutual],allPosY[idxMutual],allPosZ[idxMutual], allPhysicalValues[idxMutual],
                    &allDirectforcesX[idxMutual], &allDirectforcesY[idxMutual], &allDirectforcesZ[idxMutual], &allDirectPotentials[idxMutual]
                );
            }
        }

        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;
        offsetParticles = 0;
        groupedTree.forEachCellLeaf<FP2PGroupParticleContainer<> >([&](GroupCellClass* cellTarget, FP2PGroupParticleContainer<> * leafTarget){
            const FReal*const potentials = leafTarget->getPotentials();
            const FReal*const forcesX = leafTarget->getForcesX();
            const FReal*const forcesY = leafTarget->getForcesY();
            const FReal*const forcesZ = leafTarget->getForcesZ();
            const int nbPartsInLeafTarget = leafTarget->getNbParticles();

            for(int idxTgt = 0 ; idxTgt < nbPartsInLeafTarget ; ++idxTgt){
                potentialDiff.add(allDirectPotentials[idxTgt + offsetParticles], potentials[idxTgt]);
                fx.add(allDirectforcesX[idxTgt + offsetParticles], forcesX[idxTgt]);
                fy.add(allDirectforcesY[idxTgt + offsetParticles], forcesY[idxTgt]);
                fz.add(allDirectforcesZ[idxTgt + offsetParticles], forcesZ[idxTgt]);
            }

            offsetParticles += nbPartsInLeafTarget;
        });

        std::cout << "Error : Potential " << potentialDiff << "\n";
        std::cout << "Error : fx " << fx << "\n";
        std::cout << "Error : fy " << fy << "\n";
        std::cout << "Error : fz " << fz << "\n";
    }

    return 0;
}

