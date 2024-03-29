// ==== CMAKE =====
// @FUSE_BLAS
// ================
// Keep in private GIT
// @FUSE_MPI
// @FUSE_STARPU
// @FUSE_CUDA


#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/Core/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"
#include "../../Src/GroupTree/Rotation/FRotationCellPOD.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"
#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/GroupTree/Core/FGroupSeqAlgorithm.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskAlgorithm.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"
#include "../../Src/GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"

#include "../../Src/GroupTree/Core/FP2PGroupParticleContainer.hpp"

#include "../../Src/GroupTree/Cuda/FCudaDeviceWrapper.hpp"
#include "../../Src/GroupTree/Cuda/FCudaEmptyCellSymb.hpp"
#include "../../Src/GroupTree/Cuda/FCudaGroupOfParticles.hpp"
#include "../../Src/GroupTree/Cuda/FCudaGroupOfCells.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskStarpuMpiAlgorithm.hpp"

#include "../../Src/Files/FMpiFmaGenericLoader.hpp"
#include "../../Src/Containers/FCoordinateComputer.hpp"

#include "../../Src/GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"

#include <memory>

template <class FReal>
class FCudaP2P;

int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize { {"-bs"}, "The size of the block of the blocked tree"};
    const FParameterNames LocalOptionNoValidate { {"-validation"}, "To comparing with direct computation"};
    const FParameterNames LocalOptionStarpu { {"-starpu"}, "To compute with openmp"};
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight,FParameterDefinitions::NbParticles,
                         LocalOptionBlocSize, LocalOptionNoValidate,
                         LocalOptionStarpu);

    const bool execStarpu = FParameters::existParameter(argc, argv, LocalOptionStarpu.options);

    typedef double FReal;
    // Initialize the types
    static const int ORDER = 5;

    // Get params
    FTic timer;

    FMpi mpiComm(argc,argv);

    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);

    // init particles position and physical value
    struct TestParticle{
        FPoint<FReal> position;
        FReal physicalValue;
        const FPoint<FReal>& getPosition(){
            return position;
        }
    };

    // open particle file
    const int totalParticles = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 2000);
    const int localParticles = totalParticles/mpiComm.global().processCount();

    std::cout << " totalParticles " << totalParticles << std::endl;
    std::cout << " localParticles " << localParticles << std::endl;

    FRandomLoader<FReal> loader(localParticles,
                                1.0, FPoint<FReal>(0,0,0), mpiComm.global().processId());
    FAssertLF(loader.isOpen());

    timer.tic();
    TestParticle* allParticles = new TestParticle[loader.getNumberOfParticles()];
    memset(allParticles,0,(unsigned int) (sizeof(TestParticle)* loader.getNumberOfParticles()));
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&allParticles[idxPart].position);
        allParticles[idxPart].physicalValue = 1.0;
    }
    timer.tac();
    std::cout << "(@Create = " << timer.elapsed() << "s)." << std::endl;

    FVector<TestParticle> myParticles;
    FLeafBalance balancer;
    timer.tic();
    FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(mpiComm.global(),allParticles,
                                                                      loader.getNumberOfParticles(),
                                                                      loader.getCenterOfBox(),
                                                                      loader.getBoxWidth(),TreeHeight,
                                                                      &myParticles, &balancer);
    timer.tac();
    std::cout << "(@Distribute = " << timer.elapsed() << "s)." << std::endl;

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "For a total of " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;

    timer.tic();
    // Each proc need to know the righest morton index
    const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPosition<FReal>(
                loader.getCenterOfBox(),
                loader.getBoxWidth(),
                TreeHeight,
                myParticles[myParticles.getSize()-1].position );
    const MortonIndex myLeftLimite = host.getMortonIndex();
    MortonIndex leftLimite = -1;
    if(mpiComm.global().processId() != 0){
        FMpi::Assert(MPI_Recv(&leftLimite, sizeof(leftLimite), MPI_BYTE,
                              mpiComm.global().processId()-1, 0,
                              mpiComm.global().getComm(), MPI_STATUS_IGNORE), __LINE__);
    }
    if(mpiComm.global().processId() != mpiComm.global().processCount()-1){
        FMpi::Assert(MPI_Send(const_cast<MortonIndex*>(&myLeftLimite), sizeof(myLeftLimite), MPI_BYTE,
                              mpiComm.global().processId()+1, 0,
                              mpiComm.global().getComm()), __LINE__);
    }
    timer.tac();
    std::cout << "(@exchang border = " << timer.elapsed() << "s)." << std::endl;

    FLOG(std::cout << "My last index is " << leftLimite << "\n");
    FLOG(std::cout << "My left limite is " << myLeftLimite << "\n");

    if(execStarpu){
        std::cout << " Using Starpu " << std::endl;

        typedef FRotationCellPODCore         GroupCellSymbClass;
        typedef FRotationCellPODPole<FReal,ORDER>  GroupCellUpClass;
        typedef FRotationCellPODLocal<FReal,ORDER> GroupCellDownClass;
        typedef FRotationCellPOD<FReal,ORDER>      GroupCellClass;


        typedef FP2PGroupParticleContainer<FReal>          GroupContainerClass;
        typedef FGroupTree< FReal, GroupCellClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass, 1, 4, FReal>  GroupOctreeClass;

        typedef FStarPUCudaP2PCapacities<FRotationKernel<FReal,GroupCellClass,GroupContainerClass,ORDER>> GroupKernelClass;
        typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;

        typedef FStarPUCudaWrapper<GroupKernelClass,
                FCudaEmptyCellSymb, int, int,
                FCudaGroupOfCells<FCudaEmptyCellSymb, int, int>,
                FCudaGroupOfParticles<FReal, 1, 4, FReal>, FCudaGroupAttachedLeaf<FReal, 1, 4, FReal>, FCudaP2P<FReal> > GroupCudaWrapper;

        typedef FStarPUCudaP2PCapacities<FRotationKernel<FReal,GroupCellClass,GroupContainerClass,ORDER>> GroupKernelClass;
        typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
        typedef FGroupTaskStarPUMpiAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper, GroupCudaWrapper> GroupAlgorithm;

        const int groupSize     = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

        // Put the data into the tree
        FP2PParticleContainer<FReal> myParticlesInContainer;
        for(FSize idxPart = 0 ; idxPart < myParticles.getSize() ; ++idxPart){
            myParticlesInContainer.push(myParticles[idxPart].position,
                                        myParticles[idxPart].physicalValue);
        }
        GroupOctreeClass groupedTree(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
                                     &myParticlesInContainer, true, leftLimite);
        groupedTree.printInfoBlocks();

        timer.tac();
        std::cout << "Done  " << "(@Creating and Inserting Particles = "
                  << timer.elapsed() << "s)." << std::endl;

        { // -----------------------------------------------------
            std::cout << "\nRotationyshev FMM (ORDER="<< ORDER << ") ... " << std::endl;

            // Create Matrix Kernel
            timer.tic();
            GroupKernelClass groupkernel(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
            timer.tac();
            std::cout << "Done  " << "(@Creating kernel = " << timer.elapsed() << "s)." << std::endl;

            // Run the algorithm
            GroupAlgorithm groupalgo(mpiComm.global(), &groupedTree,&groupkernel);
            timer.tic();
            groupalgo.execute();
            timer.tac();
            std::cout << "Done  " << "(@Algorithm = " << timer.elapsed() << "s)." << std::endl;
        } // -----------------------------------------------------
    }


    return 0;
}




