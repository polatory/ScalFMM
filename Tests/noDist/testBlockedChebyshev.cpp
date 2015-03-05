// ==== CMAKE =====
// @FUSE_BLAS
// ================
// Keep in private GIT
// @SCALFMM_PRIVATE


#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Kernels/Chebyshev/FChebSymKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"
#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/GroupTree/FGroupSeqAlgorithm.hpp"
#include "../../Src/GroupTree/FGroupTaskAlgorithm.hpp"
#ifdef ScalFMM_USE_OMP4
#include "../../Src/GroupTree/FGroupTaskDepAlgorithm.hpp"
#endif
#ifdef ScalFMM_USE_STARPU
#include "../../Src/GroupTree/FGroupTaskStarpuAlgorithm.hpp"
#include "../../Src/GroupTree/FStarPUKernelCapacities.hpp"
#endif
#include "../../Src/GroupTree/FP2PGroupParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../Src/GroupTree/FP2PGroupParticleContainer.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

#include <memory>


//#define RANDOM_PARTICLES

int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize { {"-bs"}, "The size of the block of the blocked tree"};
    const FParameterNames LocalOptionNoValidate { {"-no-validation"}, "To avoid comparing with direct computation"};
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,
#ifdef RANDOM_PARTICLES
                         FParameterDefinitions::NbParticles,
#else
                         FParameterDefinitions::InputFile,
#endif
                         FParameterDefinitions::NbThreads,
                         LocalOptionBlocSize, LocalOptionNoValidate);

    // Initialize the types
    static const int ORDER = 6;
    typedef FInterpMatrixKernelR MatrixKernelClass;
    typedef FChebCell<ORDER>               GroupCellClass;
    typedef FP2PGroupParticleContainer<>          GroupContainerClass;
    typedef FGroupTree< GroupCellClass, GroupContainerClass, 5, FReal>  GroupOctreeClass;
#ifdef ScalFMM_USE_STARPU
    typedef FStarPUAllYesCapacities<FChebSymKernel<GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER>> GroupKernelClass;
    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#elif defined(ScalFMM_USE_OMP4)
    typedef FChebSymKernel<GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER> GroupKernelClass;
    // Set the number of threads
    omp_set_num_threads(FParameters::getValue(argc,argv,FParameterDefinitions::NbThreads.options, omp_get_max_threads()));
    typedef FGroupTaskDepAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#else
    typedef FChebSymKernel<GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER> GroupKernelClass;
    //typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
    typedef FGroupTaskAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#endif

    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int groupSize     = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

    // Load the particles
#ifdef RANDOM_PARTICLES
    FRandomLoader loader(FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 2000), 1.0, FPoint(0,0,0), 0);
#else
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    FFmaGenericLoader loader(filename);
#endif
    FAssertLF(loader.isOpen());
    FTic timer;

    FP2PParticleContainer<> allParticles;
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint particlePosition;
        FReal physicalValue;
#ifdef RANDOM_PARTICLES
        physicalValue = 0.10;
        loader.fillParticle(&particlePosition);
#else
        loader.fillParticle(&particlePosition, &physicalValue);
#endif
        allParticles.push(particlePosition, physicalValue);
    }
    std::cout << "Particles loaded in " << timer.tacAndElapsed() << "s\n";

    // Put the data into the tree
    timer.tic();
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    groupedTree.printInfoBlocks();
    std::cout << "Tree created in " << timer.tacAndElapsed() << "s\n";

    // Run the algorithm
    const MatrixKernelClass MatrixKernel;
    GroupKernelClass groupkernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);
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

