// @FUSE_BLAS
// @FUSE_MPI
// @FUSE_STARPU
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
#include "../../Src/Kernels/Chebyshev/FChebCell.hpp" //For validation
#include "../../Src/Core/FFmmAlgorithm.hpp" //For validation

#include "../../Src/Utils/FLeafBalance.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Files/FGenerateDistribution.hpp"
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

//#define LOAD_FILE
#ifndef LOAD_FILE
	typedef FRandomLoader<FReal> LoaderClass;
#else
	typedef FFmaGenericLoader<FReal> LoaderClass;
#endif

void timeAverage(int mpi_rank, int nproc, double elapsedTime);
void sortParticle(FPoint<FReal> * allParticlesToSort, int treeHeight, int groupSize, vector<vector<int>> & sizeForEachGroup, vector<MortonIndex> & distributedMortonIndex, LoaderClass& loader, int nproc);
void createNodeRepartition(std::vector<MortonIndex> distributedMortonIndex, std::vector<std::vector<std::vector<MortonIndex>>>& nodeRepartition, int nproc, int treeHeight);
FSize getNbParticlesPerNode(FSize mpi_count, FSize mpi_rank, FSize total);

int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    const FParameterNames LocalOptionNoValidate { {"-no-validation"}, "To avoid comparing with direct computation"};
    const FParameterNames LocalOptionEllipsoid = {{"-ellipsoid"} , " non uniform distribution on  an ellipsoid of aspect ratio given by a=0.5 b=0.25 c=0.125"};
    const FParameterNames LocalOptionEllipsoidv2 = {{"-ellipsoidv2"} , " non uniform distribution on  an ellipsoid of aspect ratio given by a=0.5 b=0.25 c=0.125"};
    const FParameterNames LocalOptionPlummer = {{"-plummer"} , " (Highly non uniform) plummer distribution (astrophysics)"};
    const FParameterNames LocalOptionCube = {{"-cube", "-uniform"} , " uniform distribution on cube (default)"};
    FHelpDescribeAndExit(argc, argv, "Loutre",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbParticles,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile, LocalOptionBlocSize,
                         LocalOptionNoValidate, LocalOptionEllipsoid, LocalOptionPlummer, LocalOptionCube,
                         LocalOptionEllipsoidv2);

    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 8);

#ifndef STARPU_USE_MPI
		cout << "Pas de mpi -_-\" " << endl;
#endif
	int mpi_rank, nproc;
    FMpi mpiComm(argc,argv);
	mpi_rank = mpiComm.global().processId();
	nproc = mpiComm.global().processCount();

#ifndef LOAD_FILE
    const FSize NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(10000));
#else
    // Load the particles
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    LoaderClass loader(filename);
    FAssertLF(loader.isOpen());
	const FSize NbParticles   = loader.getNumberOfParticles();
#endif

	FPoint<FReal> * allParticlesToSort = new FPoint<FReal>[NbParticles];

	//Fill particles
#ifndef LOAD_FILE
    srand48(0);
	FReal boxWidth = 1.0;
	{
		FSize idxPart = 0;
		for(int i = 0; i < mpiComm.global().processCount(); ++i){
			FSize NbParticlesPerNode = getNbParticlesPerNode(mpiComm.global().processCount(), i, NbParticles);
			setSeed(i+1);//Add +1 so the seed given is never 0 which correspond to random
			FReal * tmpParticles = new FReal[4*NbParticlesPerNode];
            if(FParameters::existParameter(argc, argv, "-ellipsoid")) {
                std::cout << "ellipsoid\n";
                nonunifRandomPointsOnElipsoid(NbParticlesPerNode, 0.5, 0.1, tmpParticles);
            }
            else if(FParameters::existParameter(argc, argv, LocalOptionEllipsoidv2.options)) {
                std::cout << "ellipsoidv2\n";
                unifRandomPointsOnProlate(NbParticlesPerNode, boxWidth/2, boxWidth/8, tmpParticles);
            }
            else if(FParameters::existParameter(argc, argv, "-plummer")) {
                std::cout << "plummer\n";
                //The M argument is not used in the algorithm of the plummer distribution
                unifRandomPlummer(NbParticlesPerNode, boxWidth/2, tmpParticles) ;
            }
            else { //Uniform cube
                std::cout << "cube\n";
                unifRandomPointsInCube(NbParticlesPerNode, boxWidth/2, boxWidth/2, boxWidth/2, tmpParticles);
            }
			for(FSize j= 0 ; j < NbParticlesPerNode ; ++j){
				allParticlesToSort[idxPart].setPosition(tmpParticles[j*4], tmpParticles[j*4+1], tmpParticles[j*4+2]);//Same with file or not
				++idxPart;
			}
			delete[] tmpParticles;
		}
	}
	LoaderClass loader(NbParticles, boxWidth, FPoint<FReal>(0,0,0));
#else
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
		FReal physicalValue = 0.1;
        loader.fillParticle(&allParticlesToSort[idxPart], &physicalValue);//Same with file or not
    }
#endif

	std::vector<MortonIndex> distributedMortonIndex;
	vector<vector<int>> sizeForEachGroup;
	sortParticle(allParticlesToSort, NbLevels, groupSize, sizeForEachGroup, distributedMortonIndex, loader, nproc);

    FP2PParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
		FReal physicalValue = 0.1;
        allParticles.push(allParticlesToSort[idxPart], physicalValue);
	}
    // Put the data into the tree
	
	//GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles, true);
	GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles, sizeForEachGroup, true);
	//groupedTree.printInfoBlocks();

    // Run the algorithm
	int operationsToProceed = FFmmP2M | FFmmM2M | FFmmM2L | FFmmL2L | FFmmL2P | FFmmP2P;
    const MatrixKernelClass MatrixKernel;
    GroupKernelClass groupkernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);
    GroupAlgorithm groupalgo(&groupedTree,&groupkernel, distributedMortonIndex);
	mpiComm.global().barrier();
	FTic timerExecute;
	starpu_fxt_start_profiling();
	groupalgo.execute(operationsToProceed);
	mpiComm.global().barrier();
	starpu_fxt_stop_profiling();
	timerExecute.tac();
	timeAverage(mpi_rank, nproc, timerExecute.elapsed());
	
    // Validate the result
    if(FParameters::existParameter(argc, argv, LocalOptionNoValidate.options) == false){
        typedef FP2PParticleContainer<FReal> ContainerClass;
        typedef FSimpleLeaf<FReal, ContainerClass >  LeafClass;
        typedef FChebCell<FReal,ORDER> CellClass;
        typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
        typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FFmmAlgorithm<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;

        const FReal epsi = 1E-10;

        OctreeClass treeCheck(NbLevels, 2,loader.getBoxWidth(),loader.getCenterOfBox());

        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            // put in tree
            treeCheck.insert(allParticlesToSort[idxPart], 0.1);
        }

        MatrixKernelClass MatrixKernelValidation;
        KernelClass kernels(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernelValidation);
        FmmClass algorithm(&treeCheck, &kernels);
        algorithm.execute(operationsToProceed);

        groupedTree.forEachCellWithLevel([&](GroupCellClass gcell, const int level){
			if(groupalgo.isDataOwnedBerenger(gcell.getMortonIndex(), level))
			{
				const CellClass* cell = treeCheck.getCell(gcell.getMortonIndex(), level);
				if(cell == nullptr){
					std::cout << "[Empty] Error cell should exist " << gcell.getMortonIndex() << "\n";
				}
				else {
					FMath::FAccurater<FReal> diffUp;
					diffUp.add(cell->getMultipole(0), gcell.getMultipole(0), gcell.getVectorSize());
					if(diffUp.getRelativeInfNorm() > epsi || diffUp.getRelativeL2Norm() > epsi){
						std::cout << "[Up] Up is different at index " << gcell.getMortonIndex() << " level " << level << " is " << diffUp << "\n";
					}
					FMath::FAccurater<FReal> diffDown;
					diffDown.add(cell->getLocal(0), gcell.getLocal(0), gcell.getVectorSize());
					if(diffDown.getRelativeInfNorm() > epsi || diffDown.getRelativeL2Norm() > epsi){
						std::cout << "[Down] Down is different at index " << gcell.getMortonIndex() << " level " << level << " is " << diffDown << "\n";
					}
				}
			}
        });
		groupedTree.forEachCellLeaf<FP2PGroupParticleContainer<FReal> >([&](GroupCellClass gcell, FP2PGroupParticleContainer<FReal> * leafTarget){
			if(groupalgo.isDataOwnedBerenger(gcell.getMortonIndex(), NbLevels-1))
			{
				const ContainerClass* targets = treeCheck.getLeafSrc(gcell.getMortonIndex());
				if(targets == nullptr){
					std::cout << "[Empty] Error leaf should exist " << gcell.getMortonIndex() << "\n";
				}
				else{
					const FReal*const gposX = leafTarget->getPositions()[0];
					const FReal*const gposY = leafTarget->getPositions()[1];
					const FReal*const gposZ = leafTarget->getPositions()[2];
					const FSize gnbPartsInLeafTarget = leafTarget->getNbParticles();
					const FReal*const gforceX = leafTarget->getForcesX();
					const FReal*const gforceY = leafTarget->getForcesY();
					const FReal*const gforceZ = leafTarget->getForcesZ();
					const FReal*const gpotential = leafTarget->getPotentials();

					const FReal*const posX = targets->getPositions()[0];
					const FReal*const posY = targets->getPositions()[1];
					const FReal*const posZ = targets->getPositions()[2];
					const FSize nbPartsInLeafTarget = targets->getNbParticles();
					const FReal*const forceX = targets->getForcesX();
					const FReal*const forceY = targets->getForcesY();
					const FReal*const forceZ = targets->getForcesZ();
					const FReal*const potential = targets->getPotentials();

					if(gnbPartsInLeafTarget != nbPartsInLeafTarget){
						std::cout << "[Empty] Not the same number of particles at " << gcell.getMortonIndex()
								  << " gnbPartsInLeafTarget " << gnbPartsInLeafTarget << " nbPartsInLeafTarget " << nbPartsInLeafTarget << "\n";
					}else{
						FMath::FAccurater<FReal> potentialDiff;
						FMath::FAccurater<FReal> fx, fy, fz;
						for(FSize idxPart = 0 ; idxPart < nbPartsInLeafTarget ; ++idxPart){
							if(gposX[idxPart] != posX[idxPart] || gposY[idxPart] != posY[idxPart] || gposZ[idxPart] != posZ[idxPart]){
								std::cout << "[Empty] Not the same particlea at " << gcell.getMortonIndex() << " idx " << idxPart << " "
										  << gposX[idxPart] << " " << posX[idxPart] << " " << gposY[idxPart] << " " << posY[idxPart]
										  << " " << gposZ[idxPart] << " " << posZ[idxPart] << "\n";
							}
							else{
								potentialDiff.add(potential[idxPart], gpotential[idxPart]);
								fx.add(forceX[idxPart], gforceX[idxPart]);
								fy.add(forceY[idxPart], gforceY[idxPart]);
								fz.add(forceZ[idxPart], gforceZ[idxPart]);
							}
						}
						if(potentialDiff.getRelativeInfNorm() > epsi || potentialDiff.getRelativeL2Norm() > epsi){
							std::cout << "[Up] potentialDiff is different at index " << gcell.getMortonIndex() << " is " << potentialDiff << "\n";
						}
						if(fx.getRelativeInfNorm() > epsi || fx.getRelativeL2Norm() > epsi){
							std::cout << "[Up] fx is different at index " << gcell.getMortonIndex() << " is " << fx << "\n";
						}
						if(fy.getRelativeInfNorm() > epsi || fy.getRelativeL2Norm() > epsi){
							std::cout << "[Up] fy is different at index " << gcell.getMortonIndex() << " is " << fy << "\n";
						}
						if(fz.getRelativeInfNorm() > epsi || fz.getRelativeL2Norm() > epsi){
							std::cout << "[Up] fz is different at index " << gcell.getMortonIndex() << " is " << fz << "\n";
						}
					}
				}
			}
		});

        //std::cout << "Comparing is over" << std::endl;
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
			sizeForEachGroup[treeHeight-1].push_back((int)size_last);
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

FSize getNbParticlesPerNode(FSize mpi_count, FSize mpi_rank, FSize total){
	if(mpi_rank < (total%mpi_count))
		return ((total - (total%mpi_count))/mpi_count)+1;
	return ((total - (total%mpi_count))/mpi_count);
}

