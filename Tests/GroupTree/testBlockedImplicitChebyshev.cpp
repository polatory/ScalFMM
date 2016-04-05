
// Keep in private GIT
#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
using namespace std;

#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/Core/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"


#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/GroupTree/Core/FGroupSeqAlgorithm.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskStarpuImplicitAlgorithm.hpp"
#include "../../Src/GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"

#include "../../Src/GroupTree/StarPUUtils/FStarPUCpuWrapper.hpp"
#include "../../Src/GroupTree/Core/FP2PGroupParticleContainer.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskAlgorithm.hpp"
#include "../../Src/GroupTree/Chebyshev/FChebCellPOD.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymKernel.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "../../Src/BalanceTree/FLeafBalance.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"

    typedef double FReal;

    // Initialize the types
    static const int ORDER = 6;
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;

    typedef FChebCellPODCore         GroupCellSymbClass;
    typedef FChebCellPODPole<FReal,ORDER>  GroupCellUpClass;
    typedef FChebCellPODLocal<FReal,ORDER> GroupCellDownClass;
    typedef FChebCellPOD<FReal,ORDER>      GroupCellClass;


    typedef FP2PGroupParticleContainer<FReal>          GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass, 1, 4, FReal>  GroupOctreeClass;
    typedef FStarPUAllCpuCapacities<FChebSymKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER>> GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUImplicitAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper > GroupAlgorithm;

#define LOAD_FILE
#ifndef LOAD_FILE
	typedef FRandomLoader<FReal> LoaderClass;
#else
	typedef FFmaGenericLoader<FReal> LoaderClass;
#endif
void timeAverage(int mpi_rank, int nproc, double elapsedTime);
void sortParticle(FPoint<FReal> * allParticlesToSort, int treeHeight, int groupSize, vector<vector<int>> & sizeForEachGroup, vector<MortonIndex> & distributedMortonIndex, LoaderClass& loader, int nproc);
void createNodeRepartition(std::vector<MortonIndex> distributedMortonIndex, std::vector<std::vector<std::vector<MortonIndex>>>& nodeRepartition, int nproc, int treeHeight);

int main(int argc, char* argv[]){
    setenv("STARPU_NCPU","1",1);
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbParticles,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile, LocalOptionBlocSize);

    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 8);

#ifndef STARPU_USE_MPI
		cout << "Pas de mpi -_-\" " << endl;
#endif
#ifndef LOAD_FILE
    const FSize NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(10000));
	LoaderClass loader(NbParticles, 1.0, FPoint<FReal>(0,0,0), 0);
#else
    // Load the particles
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    LoaderClass loader(filename);
#endif
	int mpi_rank, nproc = 8, provided;
	MPI_Init_thread(nullptr,nullptr, MPI_THREAD_SERIALIZED, &provided);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    FAssertLF(loader.isOpen());

	FPoint<FReal> * allParticlesToSort = new FPoint<FReal>[loader.getNumberOfParticles()];
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&allParticlesToSort[idxPart]);//Same with file or not
    }

	std::vector<MortonIndex> distributedMortonIndex;
	vector<vector<int>> sizeForEachGroup;
	sortParticle(allParticlesToSort, NbLevels, groupSize, sizeForEachGroup, distributedMortonIndex, loader, nproc);

    FP2PParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
		FReal physicalValue = 0.1;
        allParticles.push(allParticlesToSort[idxPart], physicalValue);
	}
	delete allParticlesToSort;
	allParticlesToSort = nullptr;
    // Put the data into the tree
	
	//GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles, false);
	GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles, sizeForEachGroup, false);

    // Run the algorithm
	FTic timerExecute;
    const MatrixKernelClass MatrixKernel;
    GroupKernelClass groupkernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);
    GroupAlgorithm groupalgo(&groupedTree,&groupkernel, distributedMortonIndex);
	groupalgo.execute();
	double elapsedTime = timerExecute.tacAndElapsed();
	timeAverage(mpi_rank, nproc, elapsedTime);
	
    // Validate the result
    return 0;
}
void timeAverage(int mpi_rank, int nproc, double elapsedTime)
{
	if(mpi_rank == 0)
	{
		double sumElapsedTime = elapsedTime;
		std::cout << "Executing time node 0 (implicit Cheby) : " << sumElapsedTime << "s" << std::endl;
		for(int i = 1; i < nproc; ++i)
		{
			double tmp;
			MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, 0);
			sumElapsedTime += tmp;
			std::cout << "Executing time node " << i << " (implicit Cheby) : " << tmp << "s" << std::endl;
		}
		sumElapsedTime = sumElapsedTime / (double)nproc;
		std::cout << "Average time per node (implicit Cheby) : " << sumElapsedTime << "s" << std::endl;
	}
	else
	{
		MPI_Send(&elapsedTime, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
void sortParticle(FPoint<FReal> * allParticles, int treeHeight, int groupSize, vector<vector<int>> & sizeForEachGroup, vector<MortonIndex> & distributedMortonIndex, LoaderClass& loader, int nproc)
{
	//Structure pour trier
	struct ParticleSortingStruct{
		FPoint<FReal> position;
		MortonIndex mindex;
	};
	// Création d'un tableau de la structure pour trier puis remplissage du tableau
	const FSize nbParticles = loader.getNumberOfParticles();
	ParticleSortingStruct* particlesToSort = new ParticleSortingStruct[nbParticles];
	for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
		const FTreeCoordinate host = FCoordinateComputer::GetCoordinateFromPosition<FReal>(loader.getCenterOfBox(), loader.getBoxWidth(),
				treeHeight,
				allParticles[idxPart]);
		const MortonIndex particleIndex = host.getMortonIndex(treeHeight-1);
		particlesToSort[idxPart].mindex = particleIndex;
		particlesToSort[idxPart].position = allParticles[idxPart];
	}

	//Trie du nouveau tableau
	FQuickSort<ParticleSortingStruct, FSize>::QsOmp(particlesToSort, nbParticles, [](const ParticleSortingStruct& v1, const ParticleSortingStruct& v2){
			return v1.mindex <= v2.mindex;
			});
	//Replace tout dans l'ordre dans le tableau d'origine
	for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
		allParticles[idxPart] = particlesToSort[idxPart].position;
	}
	
	//Compte le nombre de feuilles
	sizeForEachGroup.resize(treeHeight);
	MortonIndex previousLeaf = -1;
	int numberOfLeaf = 0;
	for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart)
	{
		if(particlesToSort[idxPart].mindex != previousLeaf)
		{
			previousLeaf = particlesToSort[idxPart].mindex;
			++numberOfLeaf;
		}
	}

	//Calcul de la taille des groupes au niveau des feuilles
    FLeafBalance balancer;
	for(int processId = 0; processId < nproc; ++processId)
	{
		FSize size_last;
		FSize countGroup;
		FSize leafOnProcess = balancer.getRight(numberOfLeaf, nproc, processId) - balancer.getLeft(numberOfLeaf, nproc, processId);
		size_last = leafOnProcess%groupSize;
		countGroup = (leafOnProcess - size_last)/groupSize;
		for(int i = 0; i < countGroup; ++i)
			sizeForEachGroup[treeHeight-1].push_back(groupSize);
		if(size_last > 0)
			sizeForEachGroup[treeHeight-1].push_back(size_last);
	}
	
	//Calcul du working interval au niveau des feuilles
	previousLeaf = -1;
	int countLeaf = 0;
	int processId = 0;
	FSize leafOnProcess = balancer.getRight(numberOfLeaf, nproc, 0) - balancer.getLeft(numberOfLeaf, nproc, 0);
	distributedMortonIndex.push_back(previousLeaf);
	for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart)
	{
		if(particlesToSort[idxPart].mindex != previousLeaf)
		{
			previousLeaf = particlesToSort[idxPart].mindex;
			++countLeaf;
			if(countLeaf == leafOnProcess)
			{
				distributedMortonIndex.push_back(previousLeaf);
				distributedMortonIndex.push_back(previousLeaf);
				countLeaf = 0;
				++processId;
				leafOnProcess = balancer.getRight(numberOfLeaf, nproc, processId) - balancer.getLeft(numberOfLeaf, nproc, processId);
			}
		}
	}
	distributedMortonIndex.push_back(particlesToSort[nbParticles - 1].mindex);

	//Calcul des working interval à chaque niveau
	std::vector<std::vector<std::vector<MortonIndex>>> nodeRepartition;
	createNodeRepartition(distributedMortonIndex, nodeRepartition, nproc, treeHeight);

	//Pour chaque niveau calcul de la taille des groupe
	for(int idxLevel = treeHeight - 2; idxLevel >= 0; --idxLevel)
	{
		processId = 0;
		int countParticleInTheGroup = 0;
		MortonIndex previousMortonCell = -1;

		//cout << "Compute Level " << idxLevel << endl;
		for(int idxPart = 0; idxPart < nbParticles; ++idxPart)
		{
			MortonIndex mortonCell = (particlesToSort[idxPart].mindex) >> (3*(treeHeight - 1 - idxLevel));
			if(mortonCell <= nodeRepartition[idxLevel][processId][1]) //Si l'indice est dans le working interval
			{
				if(mortonCell != previousMortonCell) //Si c'est un nouvelle indice
				{
					++countParticleInTheGroup; //On le compte dans le groupe
					previousMortonCell = mortonCell;
					if(countParticleInTheGroup == groupSize) //Si le groupe est plein on ajoute le compte
					{
						sizeForEachGroup[idxLevel].push_back(groupSize);
						countParticleInTheGroup = 0;
					}
				}
			}
			else //Si l'on change d'interval de process on ajoute ce que l'on a compté
			{
				if(countParticleInTheGroup > 0)
					sizeForEachGroup[idxLevel].push_back(countParticleInTheGroup);
				countParticleInTheGroup = 1;
				previousMortonCell = mortonCell;
				++processId;
			}
		}
		if(countParticleInTheGroup > 0)
			sizeForEachGroup[idxLevel].push_back(countParticleInTheGroup);
	}
}
void createNodeRepartition(std::vector<MortonIndex> distributedMortonIndex, std::vector<std::vector<std::vector<MortonIndex>>>& nodeRepartition, int nproc, int treeHeight) {
	nodeRepartition.resize(treeHeight, std::vector<std::vector<MortonIndex>>(nproc, std::vector<MortonIndex>(2)));
	for(int node_id = 0; node_id < nproc; ++node_id){
		nodeRepartition[treeHeight-1][node_id][0] = distributedMortonIndex[node_id*2];
		nodeRepartition[treeHeight-1][node_id][1] = distributedMortonIndex[node_id*2+1];
	}
	for(int idxLevel = treeHeight - 2; idxLevel >= 0  ; --idxLevel){
		nodeRepartition[idxLevel][0][0] = nodeRepartition[idxLevel+1][0][0] >> 3;
		nodeRepartition[idxLevel][0][1] = nodeRepartition[idxLevel+1][0][1] >> 3;
		for(int node_id = 1; node_id < nproc; ++node_id){
			nodeRepartition[idxLevel][node_id][0] = FMath::Max(nodeRepartition[idxLevel+1][node_id][0] >> 3, nodeRepartition[idxLevel][node_id-1][0]+1); //Berenger phd :)
			nodeRepartition[idxLevel][node_id][1] = nodeRepartition[idxLevel+1][node_id][1] >> 3;
		}
	}
}
