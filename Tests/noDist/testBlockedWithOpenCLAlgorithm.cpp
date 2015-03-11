// Keep in private GIT
// @SCALFMM_PRIVATE

// @FUSE_STARPU
// @FUSE_OPENCL

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

#include "../../Src/GroupTree/FGroupTaskStarpuAlgorithm.hpp"

#include "../../Src/GroupTree/FP2PGroupParticleContainer.hpp"
#include "../../Src/GroupTree/FGroupTaskAlgorithm.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/GroupTree/FGroupTestParticleContainer.hpp"

#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/GroupTree/FStarPUKernelCapacities.hpp"
#include "../../Src/GroupTree/OpenCl/FTextReplacer.hpp"

struct FTestCell_Alignement{
    static const int dataUp;
    static const int dataDown;
    static const int mindex;
    static const int coord;
};

const int FTestCell_Alignement::dataUp = reinterpret_cast<std::size_t>(&((reinterpret_cast<FTestCell*>(0xF00))->dataUp)) - std::size_t(0xF00);
const int FTestCell_Alignement::dataDown = reinterpret_cast<std::size_t>(&((reinterpret_cast<FTestCell*>(0xF00))->dataDown)) - std::size_t(0xF00);
const int FTestCell_Alignement::mindex = reinterpret_cast<std::size_t>(&((reinterpret_cast<FTestCell*>(0xF00))->mortonIndex)) - std::size_t(0xF00);
const int FTestCell_Alignement::coord = reinterpret_cast<std::size_t>(&((reinterpret_cast<FTestCell*>(0xF00))->coordinate)) - std::size_t(0xF00);


// Initialize the types
class OpenCLSource{
    FTextReplacer kernelfile;
    size_t dim;

public:
    //OpenCLSource() : kernelfile("/home/berenger/Projets/ScalfmmGit/scalfmm/Src/GroupTree/OpenCl/FEmptyKernel.cl"){
    OpenCLSource() : kernelfile("/home/berenger/Projets/ScalfmmGit/scalfmm/Src/GroupTree/OpenCl/FTestKernel.cl"){
        kernelfile.replaceAll("___FReal___", "double");
        kernelfile.replaceAll("___FParticleValueClass___", "long long");
        kernelfile.replaceAll("___FCellClassSize___", sizeof(FTestCell));
        kernelfile.replaceAll("___NbAttributesPerParticle___", 2);
        kernelfile.replaceAll("___FCellUpOffset___", FTestCell_Alignement::dataUp);
        kernelfile.replaceAll("___FCellDownOffset___", FTestCell_Alignement::dataDown);
        kernelfile.replaceAll("___FCellMortonOffset___", FTestCell_Alignement::mindex);
        kernelfile.replaceAll("___FCellCoordinateOffset___", FTestCell_Alignement::coord);

        dim = 1;
    }

    const char* getKernelCode(const int /*inDevId*/){
        return kernelfile.getContent();
    }

    void releaseKernelCode(){
        kernelfile.clear();
    }

    size_t getNbDims() const {
        return 1;
    }

    const size_t* getDimSizes() const {
        return &dim;
    }
};


int main(int argc, char* argv[]){
    setenv("STARPU_NCPU","0",1);
    setenv("STARPU_NOPENCL","1",1);
    setenv("STARPU_OPENCL_ONLY_ON_CPUS","1",1);

    const FParameterNames LocalOptionBlocSize {
        {"-bs"},
        "The size of the block of the blocked tree"
    };
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles."
                         "Usually run with STARPU_NCPU=0 STARPU_NOPENCL=1 STARPU_OPENCL_ONLY_ON_CPUS=1 ./Tests/Release/testBlockedWithOpenCLAlgorithm",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::NbThreads,
                         FParameterDefinitions::NbParticles, LocalOptionBlocSize);

    typedef FTestCell                                                       GroupCellClass;
    typedef FGroupTestParticleContainer                                     GroupContainerClass;
    typedef FGroupTree< GroupCellClass, GroupContainerClass, 2, long long int>  GroupOctreeClass;
    typedef FStarPUAllCpuOpenCLCapacities<FTestKernels< GroupCellClass, GroupContainerClass >>  GroupKernelClass;
    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass
    #ifdef ScalFMM_ENABLE_CUDA_KERNEL
        , FCudaGroupOfCells<0>, FCudaGroupOfParticles<0, int>, FCudaGroupAttachedLeaf<0, int>, FCudaEmptyKernel<>
    #endif
        , FOpenCLDeviceWrapper<GroupKernelClass, OpenCLSource>
         > GroupAlgorithm;

    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer      ContainerClass;
    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

    // FFmmAlgorithmTask FFmmAlgorithmThread
    typedef FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int groupSize      = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

//#define LOAD_FILE
#ifndef LOAD_FILE
    const int NbParticles   = FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 20);
    FRandomLoader loader(NbParticles, 1.0, FPoint(0,0,0), 0);
#else
    // Load the particles
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    FFmaGenericLoader loader(filename);
#endif
    FAssertLF(loader.isOpen());

    // Usual octree
    OctreeClass tree(NbLevels, 2, loader.getBoxWidth(), loader.getCenterOfBox());

    FP2PParticleContainer<> allParticles;
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint particlePosition;
#ifndef LOAD_FILE
        loader.fillParticle(&particlePosition);
#else
        FReal ph;
        loader.fillParticle(&particlePosition, &ph);
#endif
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
    groupalgo.execute(); // FFmmP2M TODO

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
            if(dataDown[idxPart] != loader.getNumberOfParticles()-1){
                std::cout << "[Full] Error a particle has " << dataDown[idxPart] << " (it should be " << (loader.getNumberOfParticles()-1) << ") at index " << cell->getMortonIndex() << "\n";
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
