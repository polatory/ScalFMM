#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/GroupTree/FGroupSeqAlgorithm.hpp"
#ifdef ScalFMM_USE_OMP4
#include "../../Src/GroupTree/FGroupTaskDepAlgorithm.hpp"
#endif
#ifdef ScalFMM_USE_STARPU
#include "../../Src/GroupTree/FGroupTaskStarpuAlgorithm.hpp"
#endif
#include "../../Src/GroupTree/FP2PGroupParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../Src/GroupTree/FGroupTestParticleContainer.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbThreads,
                         FParameterDefinitions::NbParticles, LocalOptionBlocSize);
    // Initialize the types
    typedef FTestCell                                                       GroupCellClass;
    typedef FGroupTestParticleContainer                                     GroupContainerClass;
    typedef FGroupTree< GroupCellClass, GroupContainerClass, 2, long long int>  GroupOctreeClass;
    typedef FTestKernels< GroupCellClass, GroupContainerClass >                       GroupKernelClass;
#ifdef ScalFMM_USE_STARPU
    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#elif defined(ScalFMM_USE_OMP4)
    // Set the number of threads
    omp_set_num_threads(FParameters::getValue(argc,argv,FParameterDefinitions::NbThreads.options, omp_get_max_threads()));
    typedef FGroupTaskDepAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#else
    //typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
    typedef FGroupTaskAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#endif

    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer      ContainerClass;
    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

    // FFmmAlgorithmTask FFmmAlgorithmThread
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 20);
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

    // Load the particles
    FRandomLoader loader(NbParticles, 1.0, FPoint(0,0,0), 0);
    FAssertLF(loader.isOpen());

    // Usual octree
    OctreeClass tree(NbLevels, 2, loader.getBoxWidth(), loader.getCenterOfBox());

    FP2PParticleContainer<> allParticles;
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint particlePosition;
        loader.fillParticle(&particlePosition);
        allParticles.push(particlePosition);
        tree.insert(particlePosition);
    }

    // Put the data into the tree
    //GroupOctreeClass groupedTree(NbLevels, groupSize, &tree);
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    groupedTree.printInfoBlocks();

    // Check tree structure at leaf level
    groupedTree.forEachCellLeaf<FGroupTestParticleContainer>([&](GroupCellClass* gcell, FGroupTestParticleContainer* gleaf){
        const ContainerClass* src = tree.getLeafSrc(gcell->getMortonIndex());
        if(src == nullptr){
            std::cout << "[PartEmpty] Error cell should not exist " << gcell->getMortonIndex() << "\n";
        }
        else {
            if(src->getNbParticles() != gleaf->getNbParticles()){
                std::cout << "[Part] Nb particles is different at index " << gcell->getMortonIndex() << " is " << gleaf->getNbParticles() << " should be " << src->getNbParticles() << "\n";
            }
        }
    });

    // Run the algorithm
    GroupKernelClass groupkernel;
    GroupAlgorithm groupalgo(&groupedTree,&groupkernel);
    groupalgo.execute();

    // Usual algorithm
    KernelClass kernels;            // FTestKernels FBasicKernels
    FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
    algo.execute();

    // Validate the result
    groupedTree.forEachCellLeaf<FGroupTestParticleContainer>([&](GroupCellClass* cell, FGroupTestParticleContainer* leaf){
        const int nbPartsInLeaf = leaf->getNbParticles();
        if(cell->getDataUp() != nbPartsInLeaf){
            std::cout << "[P2M] Error a Cell has " << cell->getDataUp() << " (it should be " << nbPartsInLeaf << ")\n";
        }
    });
    groupedTree.forEachCellLeaf<FGroupTestParticleContainer>([&](GroupCellClass* cell, FGroupTestParticleContainer* leaf){
        const int nbPartsInLeaf = leaf->getNbParticles();
        const long long int* dataDown = leaf->getDataDown();
        for(int idxPart = 0 ; idxPart < nbPartsInLeaf ; ++idxPart){
            if(dataDown[idxPart] != NbParticles-1){
                std::cout << "[Full] Error a particle has " << dataDown[idxPart] << " (it should be " << (NbParticles-1) << ") at index " << cell->getMortonIndex() << "\n";
            }
        }
    });
    // Compare the results
    groupedTree.forEachCellWithLevel([&](GroupCellClass* gcell, const int level){
        const CellClass* cell = tree.getCell(gcell->getMortonIndex(), level);
        if(cell == nullptr){
            std::cout << "[Empty] Error cell should not exist " << gcell->getMortonIndex() << "\n";
        }
        else {
            if(gcell->getDataUp() != cell->getDataUp()){
                std::cout << "[Up] Up is different at index " << gcell->getMortonIndex() << " level " << level << " is " << gcell->getDataUp() << " should be " << cell->getDataUp() << "\n";
            }
            if(gcell->getDataDown() != cell->getDataDown()){
                std::cout << "[Down] Down is different at index " << gcell->getMortonIndex() << " level " << level << " is " << gcell->getDataDown() << " should be " << cell->getDataDown() << "\n";
            }
        }
    });

    return 0;
}
