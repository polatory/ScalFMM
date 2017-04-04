// ==== CMAKE =====
// @FUSE_BLAS
// ================
// Keep in private GIT
// @FUSE_MPI


#include "../../Src/Utils/FGlobal.hpp"

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

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskStarpuMpiAlgorithm.hpp"

#include "../../Src/Files/FMpiFmaGenericLoader.hpp"
#include "../../Src/Files/FGenerateDistribution.hpp"
#include "../../Src/Containers/FCoordinateComputer.hpp"


#include <memory>

void timeAverage(int mpi_rank, int nproc, double elapsedTime);
FSize getNbParticlesPerNode(FSize mpi_count, FSize mpi_rank, FSize total);

int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize { {"-bs"}, "The size of the block of the blocked tree"};
    const FParameterNames LocalOptionNoValidate { {"-no-validation"}, "To avoid comparing with direct computation"};
    const FParameterNames LocalOptionEllipsoid = {{"-ellipsoid"} , " non uniform distribution on  an ellipsoid of aspect ratio given by a=0.5 b=0.25 c=0.125"};
    const FParameterNames LocalOptionEllipsoidv2 = {{"-ellipsoidv2"} , " non uniform distribution on  an ellipsoid of aspect ratio given by a=0.5 b=0.25 c=0.125"};
    const FParameterNames LocalOptionPlummer = {{"-plummer"} , " (Highly non uniform) plummer distribution (astrophysics)"};
    const FParameterNames LocalOptionCube = {{"-cube", "-uniform"} , " uniform distribution on cube (default)"};
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,FParameterDefinitions::InputFile,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::NbParticles,
                         LocalOptionBlocSize, LocalOptionNoValidate, LocalOptionEllipsoid, LocalOptionPlummer, LocalOptionCube,
                         LocalOptionEllipsoidv2);

    typedef double FReal;
    // Initialize the types
    static const int ORDER = 6;
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;

    // Get params
    FTic timer;

    FMpi mpiComm(argc,argv);

    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const FSize totalNbParticles = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(20));
    const FSize NbParticles   = getNbParticlesPerNode(mpiComm.global().processCount(), mpiComm.global().processId(), totalNbParticles);

    // init particles position and physical value
    struct TestParticle{
        FPoint<FReal> position;
        FReal physicalValue;
        const FPoint<FReal>& getPosition(){
            return position;
        }
		const unsigned int getWriteDataSize(void) const {
			return sizeof(FReal);
		}
		const unsigned int getWriteDataNumber(void) const {
			return 3;
		}
		const FReal* getPtrFirstData(void) const {
			return position.data();
		}
    };

//#define LOAD_FILE
#ifndef LOAD_FILE
    srand48(0);
    FReal boxWidth = 1.0;
    FRandomLoader<FReal> loader(NbParticles, boxWidth, FPoint<FReal>(0,0,0), mpiComm.global().processId());
    FAssertLF(loader.isOpen());

	setSeed(mpiComm.global().processId()+1);
    TestParticle* allParticles = new TestParticle[loader.getNumberOfParticles()];
	FReal * tmpParticles = new FReal[4*loader.getNumberOfParticles()];
    memset(allParticles,0,(unsigned int) (sizeof(TestParticle)* loader.getNumberOfParticles()));
	memset(tmpParticles,0,(unsigned int) (sizeof(FReal)* loader.getNumberOfParticles() * 4));
	if(FParameters::existParameter(argc, argv, "-ellipsoid")) {
                std::cout << "ellipsoid\n";
                nonunifRandomPointsOnElipsoid(loader.getNumberOfParticles(), 0.5, 0.1, tmpParticles);
	}
        else if(FParameters::existParameter(argc, argv, LocalOptionEllipsoidv2.options)) {
            std::cout << "ellipsoidv2\n";
            unifRandomPointsOnProlate(loader.getNumberOfParticles(), boxWidth/2, boxWidth/8, tmpParticles);
        }
	else if(FParameters::existParameter(argc, argv, "-plummer")) {
		//The M argument is not used in the algorithm of the plummer distribution
                std::cout << "plummer\n";
		unifRandomPlummer(loader.getNumberOfParticles(), boxWidth/2, tmpParticles) ;
	}
	else { //Uniform cube
                std::cout << "cube\n";
		unifRandomPointsInCube(loader.getNumberOfParticles(), boxWidth/2, boxWidth/2, boxWidth/2, tmpParticles);
	}

    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
		allParticles[idxPart].position.setPosition(tmpParticles[idxPart*4], tmpParticles[idxPart*4+1], tmpParticles[idxPart*4+2]);
		allParticles[idxPart].physicalValue = 0.1;
    }
	delete[] tmpParticles;
#else
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    FMpiFmaGenericLoader<FReal> loader(filename,mpiComm.global());
    FAssertLF(loader.isOpen());

    TestParticle* allParticles = new TestParticle[loader.getMyNumberOfParticles()];
    memset(allParticles,0,(unsigned int) (sizeof(TestParticle)* loader.getMyNumberOfParticles()));
    for(FSize idxPart = 0 ; idxPart < loader.getMyNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&allParticles[idxPart].position,&allParticles[idxPart].physicalValue);
    }
#endif

    FVector<TestParticle> myParticles;
    FLeafBalance balancer;
    FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(mpiComm.global(),allParticles,
                                                                loader.getNumberOfParticles(),
                                                                loader.getCenterOfBox(),
                                                                loader.getBoxWidth(),TreeHeight,
                                                                &myParticles, &balancer);

    //std::cout << "\tHeight : " << TreeHeight << " \t sub-height : " << SubTreeHeight << std::endl;

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

    // Put the data into the tree
    FP2PParticleContainer<FReal> myParticlesInContainer;
    for(FSize idxPart = 0 ; idxPart < myParticles.getSize() ; ++idxPart){
        myParticlesInContainer.push(myParticles[idxPart].position,
                                    myParticles[idxPart].physicalValue);
    }
    timer.tac();

	typedef FP2PParticleContainer<FReal> ContainerClass;
	typedef FSimpleLeaf<FReal, ContainerClass >  LeafClass;
	typedef FChebCell<FReal,ORDER> CellClass;
	typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
	typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
	typedef FFmmAlgorithmThreadProc<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

	const FReal epsi = 1E-10;

	OctreeClass treeCheck(TreeHeight, SubTreeHeight,loader.getBoxWidth(),loader.getCenterOfBox());

	for(FSize idxPart = 0 ; idxPart < myParticles.getSize() ; ++idxPart){
		// put in tree
		treeCheck.insert(myParticles[idxPart].position,
				myParticles[idxPart].physicalValue);
	}

	MatrixKernelClass MatrixKernel;
	KernelClass kernels(TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);
	FmmClass algorithm(mpiComm.global(),&treeCheck, &kernels);
	mpiComm.global().barrier();
	timer.tic();
	algorithm.execute();
	mpiComm.global().barrier();
	timer.tac();
	timeAverage(mpiComm.global().processId(), mpiComm.global().processCount(), timer.elapsed());

    if(FParameters::existParameter(argc, argv, LocalOptionNoValidate.options) == false){
    }

    return 0;
}


void timeAverage(int mpi_rank, int nproc, double elapsedTime)
{
	if(mpi_rank == 0)
	{
                double sumElapsedTimeMin = elapsedTime;
                double sumElapsedTimeMax = elapsedTime;
		for(int i = 1; i < nproc; ++i)
                {
			double tmp;
                        MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        if(tmp < sumElapsedTimeMin)
                                sumElapsedTimeMin = tmp;
                        if(tmp > sumElapsedTimeMax)
                                sumElapsedTimeMax = tmp;
		}
                std::cout << "Min time per node (MPI)  : " << sumElapsedTimeMin << "s" << std::endl;
                std::cout << "Max time per node (MPI)  : " << sumElapsedTimeMax << "s" << std::endl;
	}
	else
	{
		MPI_Send(&elapsedTime, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
FSize getNbParticlesPerNode(FSize mpi_count, FSize mpi_rank, FSize total){
	if(mpi_rank < (total%mpi_count))
		return ((total - (total%mpi_count))/mpi_count)+1;
	return ((total - (total%mpi_count))/mpi_count);
}

