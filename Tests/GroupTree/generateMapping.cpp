
// Keep in private GIT
// @FUSE_MPI
// @FUSE_STARPU

#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FMpi.hpp"

#include "../../Src/GroupTree/Core/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"


#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"
#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/GroupTree/Core/FGroupTaskStarpuMpiAlgorithm.hpp"

#include "../../Src/GroupTree/Core/FP2PGroupParticleContainer.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskAlgorithm.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/GroupTree/TestKernel/FGroupTestParticleContainer.hpp"

#include "../../Src/GroupTree/TestKernel/FTestCellPOD.hpp"

#include "../../Src/BalanceTree/FLeafBalance.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Containers/FCoordinateComputer.hpp"

#include "../../Src/GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"
#include "../../Src/GroupTree/StarPUUtils/FStarPUCpuWrapper.hpp"

#include <vector>
#include <iostream>
#include <fstream>


typedef double FReal;
// Initialize the types
typedef FTestCellPODCore  GroupCellSymbClass;
typedef FTestCellPODData  GroupCellUpClass;
typedef FTestCellPODData  GroupCellDownClass;
typedef FTestCellPOD      GroupCellClass;

typedef FGroupTestParticleContainer<FReal>                                     GroupContainerClass;
typedef FGroupTree< FReal, GroupCellClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass,
		GroupContainerClass, 0, 1, long long int>  GroupOctreeClass;
typedef FStarPUAllCpuCapacities<FTestKernels< GroupCellClass, GroupContainerClass >>  GroupKernelClass;
typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
typedef FGroupTaskStarPUMpiAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper> GroupAlgorithm;

void dumpTreeInfo(GroupOctreeClass& groupedTree, int rank);

int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::NbParticles,
                         LocalOptionBlocSize);
	

    FMpi mpiComm(argc, argv);
    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const FSize NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(20));
    const FSize totalNbParticles = (NbParticles*mpiComm.global().processCount());
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 8);

    // Load the particles
    FRandomLoader<FReal> loader(NbParticles, 1.0, FPoint<FReal>(0,0,0), mpiComm.global().processId());
    FAssertLF(loader.isOpen());

    // Fill the particles
    struct TestParticle{
        FPoint<FReal> position;
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

	//Create array of particle and fill it with 0
    std::unique_ptr<TestParticle[]> particles(new TestParticle[loader.getNumberOfParticles()]);
    memset(particles.get(), 0, sizeof(TestParticle) * loader.getNumberOfParticles());

	//Fill with particles
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
		MortonIndex id = mpiComm.global().processId()*loader.getNumberOfParticles() + idxPart;
        loader.fillParticleAtMortonIndex(&particles[idxPart].position, id, NbLevels);
		const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPosition<FReal>(loader.getCenterOfBox(), loader.getBoxWidth(), NbLevels, particles[idxPart].position );
		const MortonIndex particleIndex = host.getMortonIndex(NbLevels);
		if(particleIndex != id )
			std::cout << "Error (" << mpiComm.global().processId() << ") " << id << " - " << particleIndex << std::endl;
    }
    // Sort in parallel
    FVector<TestParticle> myParticles;
    FLeafBalance balancer;
    FMpiTreeBuilder<FReal, TestParticle >::DistributeArrayToContainer(mpiComm.global(),
                                                                particles.get(),
                                                                loader.getNumberOfParticles(),
                                                                loader.getCenterOfBox(),
                                                                loader.getBoxWidth(),
                                                                NbLevels,
                                                                &myParticles,
                                                                &balancer);

    FTestParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < myParticles.getSize() ; ++idxPart){
        allParticles.push(myParticles[idxPart].position);
    }

	for(int i = 0; i < mpiComm.global().processCount(); ++i)
	{
		if(i == mpiComm.global().processId())
			std::cout << mpiComm.global().processId() << " has " << myParticles.getSize() << " particles" << std::endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}


    // Each proc need to know the righest morton index
    const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPosition<FReal>(
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
        FMpi::Assert(MPI_Send(const_cast<MortonIndex*>(&myLeftLimite), sizeof(myLeftLimite), MPI_BYTE,
                              mpiComm.global().processId()+1, 0,
                              mpiComm.global().getComm()), __LINE__);
    }
	if(mpiComm.global().processId() == 0){
		std::vector<TestParticle*> particlesGathered;
		std::vector<int> sizeGathered;
		std::vector<MortonIndex> mortonIndex;
		
		//Ajout des mes particules
		int sizeofParticle = sizeof(TestParticle)*myParticles.getSize();
		sizeGathered.push_back(sizeofParticle);
		particlesGathered.push_back(new TestParticle[sizeofParticle]);
		memcpy(particlesGathered.back(), myParticles.data(), sizeofParticle);
		mortonIndex.push_back(leftLimite);
		mortonIndex.push_back(myLeftLimite);
		//Recupération des particules des autres
		for(int i = 1; i < mpiComm.global().processCount(); ++i)
		{
			int sizeReceive;
			MortonIndex mortonStart, mortonEnd;
			MPI_Recv(&sizeReceive, sizeof(sizeReceive), MPI_BYTE, i, 0, mpiComm.global().getComm(), MPI_STATUS_IGNORE);
			sizeGathered.push_back(sizeReceive);
			particlesGathered.push_back(new TestParticle[sizeReceive]);
			MPI_Recv(particlesGathered.back(), sizeReceive, MPI_BYTE, i, 0, mpiComm.global().getComm(), MPI_STATUS_IGNORE);
			MPI_Recv(&mortonStart, sizeof(mortonStart), MPI_BYTE, i, 0, mpiComm.global().getComm(), MPI_STATUS_IGNORE);
			MPI_Recv(&mortonEnd, sizeof(mortonEnd), MPI_BYTE, i, 0, mpiComm.global().getComm(), MPI_STATUS_IGNORE);
			mortonIndex.push_back(mortonStart);
			mortonIndex.push_back(mortonEnd);
		}
		int sum = 0;
		for(int a : sizeGathered)
			sum += a/sizeof(TestParticle);
		if(sum != totalNbParticles)
			std::cout << "Erreur sum : " << sum << " instead of " << totalNbParticles << std::endl;
		//Store in that bloody file
		FFmaGenericWriter<FReal> writer("canard.fma");
		writer.writeHeader(loader.getCenterOfBox(), loader.getBoxWidth(),totalNbParticles, particles[0]);
		for(unsigned int i = 0; i < particlesGathered.size(); ++i)
			writer.writeArrayOfParticles(particlesGathered[i], sizeGathered[i]/sizeof(TestParticle));
		std::ofstream fichier("mapping", std::ios::out | std::ios::trunc);  //déclaration du flux et ouverture du fichier
		fichier << mortonIndex.size()/2 << std::endl;
		for(unsigned int i = 0; i < mortonIndex.size(); i+=2)
			fichier << mortonIndex[i] << " " << mortonIndex[i+1] << std::endl;
		for(TestParticle* ptr : particlesGathered)
			delete ptr;
	}
	else{
		int sizeofParticle = sizeof(TestParticle)*myParticles.getSize();
        MPI_Send(&sizeofParticle, sizeof(sizeofParticle), MPI_BYTE, 0, 0, mpiComm.global().getComm());//Send size
        MPI_Send(myParticles.data(), sizeofParticle, MPI_BYTE, 0, 0, mpiComm.global().getComm());
        MPI_Send(const_cast<MortonIndex*>(&leftLimite), sizeof(leftLimite), MPI_BYTE, 0, 0, mpiComm.global().getComm());
        MPI_Send(const_cast<MortonIndex*>(&myLeftLimite), sizeof(myLeftLimite), MPI_BYTE, 0, 0, mpiComm.global().getComm());
	}
	// Put the data into the tree
	GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize,
			&allParticles, true, leftLimite);

	//Print morton index for each node
	//for(int i = 0; i < mpiComm.global().processCount(); ++i)
	//{
		//if(i == mpiComm.global().processId())
		//{
			//FLOG(std::cout << "My last index is " << leftLimite << "\n");
			//FLOG(std::cout << "My left limite is " << myLeftLimite << "\n");
			//std::cout << "My last index is (" << mpiComm.global().processId() << ") " << leftLimite << "\n";
			//std::cout << "My left limite is (" << mpiComm.global().processId() << ") " << myLeftLimite << "\n";
			//std::cout << "Size (" << mpiComm.global().processId() << ") " << allParticles.getNbParticles()  << "\n";
		//}
		//MPI_Barrier(MPI_COMM_WORLD);
	//}
	//if(mpiComm.global().processId() == 0)
	//{
		//std::cout << "NbLevel " << NbLevels << std::endl;
		//remove("mapping_morton");
	//}
	//for(int i = 0; i < mpiComm.global().processCount(); ++i)
	//{
		//if(i == mpiComm.global().processId())
			//dumpTreeInfo(groupedTree, i);
		//MPI_Barrier(MPI_COMM_WORLD);
	//}
    // Run the algorithm
    GroupKernelClass groupkernel;
    GroupAlgorithm groupalgo(mpiComm.global(), &groupedTree,&groupkernel);
    groupalgo.execute();

    std::cout << "Wait Others... " << std::endl;
    mpiComm.global().barrier();

    groupedTree.forEachCellLeaf<GroupContainerClass>([&](GroupCellClass cell, GroupContainerClass* leaf){
        const FSize nbPartsInLeaf = leaf->getNbParticles();
        const long long int* dataDown = leaf->getDataDown();
        for(FSize idxPart = 0 ; idxPart < nbPartsInLeaf ; ++idxPart){
            if(dataDown[idxPart] != totalNbParticles-1){
                std::cout << "[Full] Error a particle has " << dataDown[idxPart] << " (it should be " << (totalNbParticles-1) << ") at index " << cell.getMortonIndex() << "\n";
            }
        }
    });


	//Check if the algorithm was good
    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer<FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    {
        // Usual octree
        OctreeClass tree(NbLevels, 2, loader.getBoxWidth(), loader.getCenterOfBox());
        for(int idxProc = 0 ; idxProc < mpiComm.global().processCount() ; ++idxProc){
            FRandomLoader<FReal> loaderAll(NbParticles, 1.0, FPoint<FReal>(0,0,0), idxProc);
            for(FSize idxPart = 0 ; idxPart < loaderAll.getNumberOfParticles() ; ++idxPart){
                FPoint<FReal> pos;
				MortonIndex id = idxProc*loader.getNumberOfParticles() + idxPart;
                loaderAll.fillParticleAtMortonIndex(&pos, id, NbLevels);
                tree.insert(pos);
            }
        }
        // Usual algorithm
        KernelClass kernels;            // FTestKernels FBasicKernels
        FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
        algo.execute();

        // Compare the results
        groupedTree.forEachCellWithLevel([&](GroupCellClass gcell, const int level){
            const CellClass* cell = tree.getCell(gcell.getMortonIndex(), level);
            if(cell == nullptr){
                std::cout << "[Empty] Error cell should not exist " << gcell.getMortonIndex() << "\n";
            }
            else {
                if(gcell.getDataUp() != cell->getDataUp()){
                    std::cout << "[Up] Up is different at index " << gcell.getMortonIndex() << " level " << level << " is " << gcell.getDataUp() << " should be " << cell->getDataUp() << "\n";
                }
                if(gcell.getDataDown() != cell->getDataDown()){
                    std::cout << "[Down] Down is different at index " << gcell.getMortonIndex() << " level " << level << " is " << gcell.getDataDown() << " should be " << cell->getDataDown() << "\n";
                }
            }
        });
    }

    return 0;
}
void dumpTreeInfo(GroupOctreeClass& groupedTree, int rank)
{
	std::cout << "Dump tree " << rank << std::endl;
	std::cout << "NbGrPar " << groupedTree.getNbParticleGroup() << std::endl;
	std::ofstream fichier("mapping_morton", std::ios::out | std::ios::ate);  //déclaration du flux et ouverture du fichier
	groupedTree.forEachCellWithLevel([&](GroupCellClass gcell, int idxLevel){
					fichier << "c " <<  idxLevel << " " << gcell.getMortonIndex() << " " << rank << std::endl;
			});
	for(int idxParticleGroup = 0 ; idxParticleGroup < groupedTree.getNbParticleGroup() ; ++idxParticleGroup){
		GroupOctreeClass::ParticleGroupClass* currentGroupParticle = groupedTree.getParticleGroup(idxParticleGroup);
		for(MortonIndex idxParticle = currentGroupParticle->getStartingIndex(); idxParticle < currentGroupParticle->getEndingIndex(); ++idxParticle) {
			if(currentGroupParticle->exists(idxParticle))
				fichier << "p " << idxParticle << " " << rank << std::endl;
		}
	}
	fichier << std::flush;
	fichier.close();
}

