
// Keep in private GIT
// @SCALFMM_PRIVATE
// @SCALFMM_USE_MPI
// @SCALFMM_USE_STARPU

#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FMpi.hpp"

#include "../../Src/GroupTree/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/GroupTree/FGroupTaskStarpuMpiAlgorithm.hpp"

#include "../../Src/GroupTree/FP2PGroupParticleContainer.hpp"
#include "../../Src/GroupTree/FGroupTaskAlgorithm.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../Src/GroupTree/FGroupTestParticleContainer.hpp"

#include "../../Src/BalanceTree/FLeafBalance.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"


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
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbThreads,
                         FParameterDefinitions::NbParticles, FParameterDefinitions::NbThreads,
                         LocalOptionBlocSize);
    // Initialize the types
    typedef FTestCell                                                       GroupCellClass;
    typedef FGroupTestParticleContainer                                     GroupContainerClass;
    typedef FGroupTree< GroupCellClass, GroupContainerClass, 2, long long int>  GroupOctreeClass;
    typedef FTestKernels< GroupCellClass, GroupContainerClass >                       GroupKernelClass;
    typedef FGroupTaskStarPUMpiAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;

    FMpi mpiComm(argc, argv);
    // Get params
    const int maxThreads    = FParameters::getValue(argc,argv,FParameterDefinitions::NbThreads.options, -1);
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 20);
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);
    const int totalNbParticles = (NbParticles*mpiComm.global().processCount());

    // Load the particles
    FRandomLoader loader(NbParticles, 1.0, FPoint(0,0,0), mpiComm.global().processId());
    FAssertLF(loader.isOpen());

    // Fill the particles
    struct TestParticle{
        FPoint position;
        const FPoint& getPosition(){
            return position;
        }
    };

    std::unique_ptr<TestParticle[]> particles(new TestParticle[loader.getNumberOfParticles()]);
    memset(particles.get(), 0, sizeof(TestParticle) * loader.getNumberOfParticles());
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&particles[idxPart].position);
    }
    // Sort in parallel
    FVector<TestParticle> myParticles;
    FLeafBalance balancer;
    FMpiTreeBuilder< TestParticle >::DistributeArrayToContainer(mpiComm.global(),
                                                                particles.get(),
                                                                loader.getNumberOfParticles(),
                                                                loader.getCenterOfBox(),
                                                                loader.getBoxWidth(),
                                                                NbLevels,
                                                                &myParticles,
                                                                &balancer);

    FP2PParticleContainer<> allParticles;
    for(int idxPart = 0 ; idxPart < myParticles.getSize() ; ++idxPart){
        allParticles.push(myParticles[idxPart].position);
    }

    // Each proc need to know the righest morton index
    const FTreeCoordinate host = getCoordinateFromPosition(
                loader.getCenterOfBox(),
                loader.getBoxWidth(),
                NbLevels,
                myParticles[myParticles.getSize()-1].position );
    const MortonIndex myLeftLimite = host.getMortonIndex(NbLevels-1);
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
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
                                 &allParticles, true, leftLimite);
    groupedTree.printInfoBlocks();

    // Run the algorithm
    GroupKernelClass groupkernel;
    GroupAlgorithm groupalgo(mpiComm.global(), &groupedTree,&groupkernel,maxThreads);
    groupalgo.execute();

    groupedTree.forEachCellLeaf<FGroupTestParticleContainer>([&](GroupCellClass* cell, FGroupTestParticleContainer* leaf){
        const int nbPartsInLeaf = leaf->getNbParticles();
        const long long int* dataDown = leaf->getDataDown();
        for(int idxPart = 0 ; idxPart < nbPartsInLeaf ; ++idxPart){
            if(dataDown[idxPart] != totalNbParticles-1){
                std::cout << "[Full] Error a particle has " << dataDown[idxPart] << " (it should be " << (totalNbParticles-1) << ") at index " << cell->getMortonIndex() << "\n";
            }
        }
    });

    return 0;
}

