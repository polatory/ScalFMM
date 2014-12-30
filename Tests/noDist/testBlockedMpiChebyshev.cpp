// ==== CMAKE =====
// @FUSE_BLAS
// ================
// Keep in private GIT
// @SCALFMM_PRIVATE
// @SCALFMM_USE_MPI
// @SCALFMM_USE_STARPU


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
#include "../../Src/GroupTree/FGroupTaskStarpuAlgorithm.hpp"
#include "../../Src/GroupTree/FP2PGroupParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../Src/GroupTree/FP2PGroupParticleContainer.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"
#include "../../Src/GroupTree/FGroupTaskStarpuMpiAlgorithm.hpp"

#include "../../Src/Files/FMpiFmaGenericLoader.hpp"

#include <memory>


int getTreeCoordinate(const FReal inRelativePosition, const FReal boxWidth,
                      const FReal boxWidthAtLeafLevel, const int treeHeight) {
    FAssertLF( (inRelativePosition >= 0 && inRelativePosition <= boxWidth), "inRelativePosition : ",inRelativePosition );
    if(inRelativePosition == boxWidth){
        return FMath::pow2(treeHeight-1)-1;
    }
    const FReal indexFReal = inRelativePosition / boxWidthAtLeafLevel;
    return static_cast<int>(indexFReal);
}

FTreeCoordinate getCoordinateFromPosition(const FPoint& centerOfBox, const FReal boxWidth, const int treeHeight,
                                          const FPoint& pos) {
    const FPoint boxCorner(centerOfBox,-(boxWidth/2));
    const FReal boxWidthAtLeafLevel(boxWidth/FReal(1<<(treeHeight-1)));

    // box coordinate to host the particle
    FTreeCoordinate host;
    // position has to be relative to corner not center
    host.setX( getTreeCoordinate( pos.getX() - boxCorner.getX(), boxWidth, boxWidthAtLeafLevel, treeHeight));
    host.setY( getTreeCoordinate( pos.getY() - boxCorner.getY(), boxWidth, boxWidthAtLeafLevel, treeHeight));
    host.setZ( getTreeCoordinate( pos.getZ() - boxCorner.getZ(), boxWidth, boxWidthAtLeafLevel, treeHeight));
    return host;
}


int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize { {"-bs"}, "The size of the block of the blocked tree"};
    const FParameterNames LocalOptionNoValidate { {"-no-validation"}, "To avoid comparing with direct computation"};
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,FParameterDefinitions::InputFile,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::NbThreads,
                         LocalOptionBlocSize, LocalOptionNoValidate);

    // Initialize the types
    static const int ORDER = 6;
    typedef FInterpMatrixKernelR MatrixKernelClass;
    typedef FChebCell<ORDER>               GroupCellClass;
    typedef FP2PGroupParticleContainer<>          GroupContainerClass;
    typedef FGroupTree< GroupCellClass, GroupContainerClass, 5, FReal>  GroupOctreeClass;
    typedef FChebSymKernel<GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER> GroupKernelClass;
    typedef FGroupTaskStarPUMpiAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;

    // Get params
    FTic timer;
    const int groupSize     = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

    FMpi mpiComm(argc,argv);

    const char* const filename       = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const unsigned int NbThreads     = FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);

    // init particles position and physical value
    struct TestParticle{
        FPoint position;
        FReal physicalValue;
        const FPoint& getPosition(){
            return position;
        }
    };

    // open particle file
    std::cout << "Opening : " << filename << "\n" << std::endl;
    FMpiFmaGenericLoader loader(filename,mpiComm.global());
    FAssertLF(loader.isOpen());

    TestParticle* allParticles = new TestParticle[loader.getMyNumberOfParticles()];
    memset(allParticles,0,(unsigned int) (sizeof(TestParticle)* loader.getMyNumberOfParticles()));
    for(int idxPart = 0 ; idxPart < loader.getMyNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&allParticles[idxPart].position,&allParticles[idxPart].physicalValue);
    }

    FVector<TestParticle> myParticles;
    FLeafBalance balancer;
    FMpiTreeBuilder< TestParticle >::DistributeArrayToContainer(mpiComm.global(),allParticles,
                                                                loader.getMyNumberOfParticles(),
                                                                loader.getCenterOfBox(),
                                                                loader.getBoxWidth(),TreeHeight,
                                                                &myParticles, &balancer);

    std::cout << "Creating & Inserting " << loader.getMyNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "For a total of " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;

    // Each proc need to know the righest morton index
    const FTreeCoordinate host = getCoordinateFromPosition(
                loader.getCenterOfBox(),
                loader.getBoxWidth(),
                TreeHeight,
                myParticles[myParticles.getSize()-1].position );
    const MortonIndex myLeftLimite = host.getMortonIndex(TreeHeight-1);
    MortonIndex leftLimite = -1;
    if(mpiComm.global().processId() != 0){
        FMpi::Assert(MPI_Recv(&leftLimite, sizeof(leftLimite), MPI_BYTE,
                              mpiComm.global().processId()-1, 0,
                              mpiComm.global().getComm(), MPI_STATUS_IGNORE), __LINE__);
    }
    if(mpiComm.global().processId() != mpiComm.global().processCount()-1){
        FMpi::Assert(MPI_Send((void*)&myLeftLimite, sizeof(myLeftLimite), MPI_BYTE,
                              mpiComm.global().processId()+1, 0,
                              mpiComm.global().getComm()), __LINE__);
    }
    FLOG(std::cout << "My last index is " << leftLimite << "\n");
    FLOG(std::cout << "My left limite is " << myLeftLimite << "\n");


    // Put the data into the tree
    FP2PParticleContainer<> myParticlesInContainer;
    for(int idxPart = 0 ; idxPart < myParticles.getSize() ; ++idxPart){
        myParticlesInContainer.push(myParticles[idxPart].position);
    }
    GroupOctreeClass groupedTree(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
                                 &myParticlesInContainer, true, leftLimite);
    groupedTree.printInfoBlocks();

    timer.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = "
              << timer.elapsed() << "s)." << std::endl;

    { // -----------------------------------------------------
        std::cout << "\nChebyshev FMM (ORDER="<< ORDER << ") ... " << std::endl;
        timer.tic();

        // Create Matrix Kernel
        const MatrixKernelClass MatrixKernel;
        GroupKernelClass groupkernel(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);
        // Run the algorithm
        GroupAlgorithm groupalgo(mpiComm.global(), &groupedTree,&groupkernel,NbThreads);
        groupalgo.execute();

        timer.tac();
        std::cout << "Done  " << "(@Algorithm = " << timer.elapsed() << "s)." << std::endl;
    } // -----------------------------------------------------


    if(FParameters::existParameter(argc, argv, LocalOptionNoValidate.options) == false){
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;

        groupedTree.forEachCellLeaf<FP2PGroupParticleContainer<> >([&](GroupCellClass* /*cellTarget*/, FP2PGroupParticleContainer<> * leafTarget){
            const FReal*const potentials = leafTarget->getPotentials();
            const FReal*const forcesX = leafTarget->getForcesX();
            const FReal*const forcesY = leafTarget->getForcesY();
            const FReal*const forcesZ = leafTarget->getForcesZ();
            const int nbPartsInLeafTarget = leafTarget->getNbParticles();

            for(int idxTgt = 0 ; idxTgt < nbPartsInLeafTarget ; ++idxTgt){
//                potentialDiff.add(allDirectPotentials[idxTgt + offsetParticles], potentials[idxTgt]);
//                fx.add(allDirectforcesX[idxTgt + offsetParticles], forcesX[idxTgt]);
//                fy.add(allDirectforcesY[idxTgt + offsetParticles], forcesY[idxTgt]);
//                fz.add(allDirectforcesZ[idxTgt + offsetParticles], forcesZ[idxTgt]);
            }
        });

        std::cout << "Error : Potential " << potentialDiff << "\n";
        std::cout << "Error : fx " << fx << "\n";
        std::cout << "Error : fy " << fy << "\n";
        std::cout << "Error : fz " << fz << "\n";
    }

    return 0;
}


