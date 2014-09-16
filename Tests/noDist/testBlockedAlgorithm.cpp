#include "../../Src/GroupTree/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/GroupTree/FGroupSeqAlgorithm.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../Src/GroupTree/FGroupTestParticleContainer.hpp"

int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::NbParticles, LocalOptionBlocSize);
    // Initialize the types
    typedef FTestCell                                                       GroupCellClass;
    typedef FGroupTestParticleContainer                                     GroupContainerClass;
    typedef FSimpleLeaf< GroupContainerClass >                                   LeafClass;
    typedef FGroupTree< GroupCellClass, GroupContainerClass, 2, long long int>  GroupOctreeClass;
    typedef FTestKernels< GroupCellClass, GroupContainerClass >                       GroupKernelClass;
    typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 20);
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

    // Load the particles
    FRandomLoader loader(NbParticles, 1.0, FPoint(0,0,0), 0);
    FAssertLF(loader.isOpen());
    FP2PParticleContainer<> allParticles;
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint particlePosition;
        loader.fillParticle(&particlePosition);
        allParticles.push(particlePosition);
    }

    // Put the data into the tree
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    groupedTree.printInfoBlocks();

    // Run the algorithm
    GroupKernelClass kernel;
    GroupAlgorithm algo(&groupedTree,&kernel);
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

    return 0;
}
