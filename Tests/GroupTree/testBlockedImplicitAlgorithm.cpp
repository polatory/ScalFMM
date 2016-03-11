
// Keep in private GIT
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/Core/FGroupTreeBerenger.hpp"

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
#include "../../Src/GroupTree/Core/FGroupTaskAlgorithm.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/GroupTree/TestKernel/FGroupTestParticleContainer.hpp"
#include "../../Src/GroupTree/TestKernel/FTestCellPOD.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"

std::vector<MortonIndex> getMortonIndex(const char* const mapping_filename);

int main(int argc, char* argv[]){
    setenv("STARPU_NCPU","1",1);
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
	const FParameterNames Mapping {
		{"-map"} ,
		"mapping  \\o/."
	};
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbParticles,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::InputFile, LocalOptionBlocSize, Mapping);
	//int provided = 0;
	//MPI_Init_thread(&argc,&argv, MPI_THREAD_SERIALIZED, &provided);

    typedef double FReal;

    // Initialize the types
    typedef FTestCellPODCore  GroupCellSymbClass;
    typedef FTestCellPODData  GroupCellUpClass;
    typedef FTestCellPODData  GroupCellDownClass;
    typedef FTestCellPOD      GroupCellClass;


    typedef FGroupTestParticleContainer<FReal>                                GroupContainerClass;
    typedef FGroupTreeBerenger< FReal, GroupCellClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass,
            GroupContainerClass, 0, 1, long long int>  GroupOctreeClass;
    typedef FStarPUAllCpuCapacities<FTestKernels< GroupCellClass, GroupContainerClass >>  GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUImplicitAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper > GroupAlgorithm;

    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer<FReal>      ContainerClass;
    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

    // FFmmAlgorithmTask FFmmAlgorithmThread
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 8);
    const char* const mapping_filename      = FParameters::getStr(argc,argv,Mapping.options, "mapping");
	std::vector<MortonIndex> distributedMortonIndex = getMortonIndex(mapping_filename);

#ifndef STARPU_USE_MPI
		cout << "Pas de mpi -_-\" " << endl;
#endif
#define LOAD_FILE
#ifndef LOAD_FILE
    const FSize NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, FSize(20));
    FRandomLoader<FReal> loader(NbParticles, 1.0, FPoint<FReal>(0,0,0), 0);
#else
    // Load the particles
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    FFmaGenericLoader<FReal> loader(filename);
#endif
    FAssertLF(loader.isOpen());

    // Usual octree
    OctreeClass tree(NbLevels, FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 2),
                     loader.getBoxWidth(), loader.getCenterOfBox());
    FTestParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
        loader.fillParticle(&particlePosition);//Same with file or not
        allParticles.push(particlePosition);
        tree.insert(particlePosition);
    }
    // Put the data into the tree
    //GroupOctreeClass groupedTree(NbLevels, groupSize, &tree);
	GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles, false, distributedMortonIndex);
	
    // Check tree structure at leaf level
    groupedTree.forEachCellLeaf<GroupContainerClass>([&](GroupCellClass gcell, GroupContainerClass* gleaf){
        const ContainerClass* src = tree.getLeafSrc(gcell.getMortonIndex());
        if(src == nullptr){
            std::cout << "[PartEmpty] Error cell should not exist " << gcell.getMortonIndex() << "\n";
        }
        else {
            if(src->getNbParticles() != gleaf->getNbParticles()){
                std::cout << "[Part] Nb particles is different at index " << gcell.getMortonIndex() << " is " << gleaf->getNbParticles() << " should be " << src->getNbParticles() << "\n";
            }
        }
    });

    // Run the algorithm
    GroupKernelClass groupkernel;
    GroupAlgorithm groupalgo(&groupedTree,&groupkernel, distributedMortonIndex);
    groupalgo.execute();
    // Usual algorithm
    KernelClass kernels;            // FTestKernels FBasicKernels
    FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
    algo.execute();
	int rank = groupalgo.getRank();
	for(int i = 0; i < groupedTree.getHeight(); ++i)
	{
		if(groupedTree.getNbCellGroupAtLevel(i) < groupalgo.getNProc())
			std::cout << "Error at level " << i << std::endl;
	}
	return 0;
    // Validate the result
	for(int idxLevel = 2 ; idxLevel < groupedTree.getHeight() ; ++idxLevel){
		for(int idxGroup = 0 ; idxGroup < groupedTree.getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
			//if(groupalgo.isDataOwned(idxGroup, groupedTree.getNbCellGroupAtLevel(idxLevel))){
			if(groupalgo.isDataOwnedBerenger(groupedTree.getCellGroup(idxLevel, idxGroup)->getStartingIndex(), idxLevel)){
				GroupOctreeClass::CellGroupClass* currentCells = groupedTree.getCellGroup(idxLevel, idxGroup);
				currentCells->forEachCell([&](GroupCellClass gcell){
						const CellClass* cell = tree.getCell(gcell.getMortonIndex(), idxLevel);
						if(cell == nullptr){
							std::cout << "[Empty node(" << rank << ")] Error cell should not exist " << gcell.getMortonIndex() << "\n";
						}
						else {
							if(gcell.getDataUp() != cell->getDataUp()){
								std::cout << "[Up node(" << rank << ")] Up is different at index " << gcell.getMortonIndex() << " level " << idxLevel << " is " << gcell.getDataUp() << " should be " << cell->getDataUp() << "\n";
							}
							if(gcell.getDataDown() != cell->getDataDown()){
								std::cout << "[Down node(" << rank << ")] Down is different at index " << gcell.getMortonIndex() << " level " << idxLevel << " is " << gcell.getDataDown() << " should be " << cell->getDataDown() << "\n";
							}
						}
				});
			}
		}
	}
	{
		int idxLevel = groupedTree.getHeight()-1;
		for(int idxGroup = 0 ; idxGroup < groupedTree.getNbCellGroupAtLevel(idxLevel) ; ++idxGroup){
			//if(groupalgo.isDataOwned(idxGroup, groupedTree.getNbCellGroupAtLevel(idxLevel))){
			if(groupalgo.isDataOwnedBerenger(groupedTree.getCellGroup(groupedTree.getHeight()-1, idxGroup)->getStartingIndex(), groupedTree.getHeight()-1)){
				GroupOctreeClass::ParticleGroupClass* particleGroup = groupedTree.getParticleGroup(idxGroup); 
				GroupOctreeClass::CellGroupClass* cellGroup = groupedTree.getCellGroup(idxLevel, idxGroup);
				cellGroup->forEachCell([&](GroupCellClass cell){
					MortonIndex midx = cell.getMortonIndex();
					const int leafIdx = particleGroup->getLeafIndex(midx);
					GroupContainerClass leaf = particleGroup->template getLeaf<GroupContainerClass>(leafIdx);
					const FSize nbPartsInLeaf = leaf.getNbParticles();
					if(cell.getDataUp() != nbPartsInLeaf){
						std::cout << "[P2M node(" << rank << ")] Error a Cell has " << cell.getDataUp() << " (it should be " << nbPartsInLeaf << ")\n";
					}
					const long long int* dataDown = leaf.getDataDown();
					for(FSize idxPart = 0 ; idxPart < nbPartsInLeaf ; ++idxPart){
						if(dataDown[idxPart] != loader.getNumberOfParticles()-1){
							std::cout << "[Full node(" << rank << ")] Error a particle has " << dataDown[idxPart] << " (it should be " << (loader.getNumberOfParticles()-1) << ") at index " << cell.getMortonIndex() << "\n";
						}
					}
				});
			}
		}
	}
    return 0;
}
std::vector<MortonIndex> getMortonIndex(const char* const mapping_filename)
{
	std::vector<MortonIndex> ret;

	std::ifstream fichier(mapping_filename, ios::in);  // on ouvre le fichier en lecture

	if(fichier)  // si l'ouverture a réussi
	{       
		int nbProcess;
		fichier >> nbProcess;
		for(int i = 0; i < nbProcess; ++i)
		{
			MortonIndex start, end;
			fichier >> start >> end;
			ret.push_back(start);
			ret.push_back(end);
		}
		// instructions
		fichier.close();  // on ferme le fichier
	}
	else  // sinon
		cerr << "Impossible d'ouvrir le fichier !" << endl;
	return ret;
}
