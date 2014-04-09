
#include "../../Src/GroupTree/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Files/FFmaBinLoader.hpp"

#include "../../Src/GroupTree/FGroupSeqAlgorithm.hpp"
#include "../../Src/GroupTree/FP2PGroupParticleContainer.hpp"

int main(int argc, char* argv[]){
    static const int P = 9;
    typedef FRotationCell<P>               CellClass;
    typedef FP2PParticleContainer<>          ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FGroupTree< CellClass, FP2PGroupParticleContainer<>, 4, FReal>  GroupOctreeClass;

    FTic counter;
    const int NbLevels      = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.bin.fma.double");

    FFmaBinLoader loader(filename);
    FAssertLF(loader.isOpen());

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    FP2PParticleContainer<> allParticles;

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint particlePosition;
        FReal physicalValue;
        loader.fillParticle(&particlePosition,&physicalValue);
        tree.insert(particlePosition, physicalValue );
        allParticles.push(particlePosition, physicalValue);
    }

    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.tacAndElapsed() << "s)." << std::endl;

    const int groupSize      = FParameters::getValue(argc,argv,"-bs", 250);

    counter.tic();
    GroupOctreeClass groupedTree2(NbLevels, groupSize, &tree);
    std::cout << "Done  " << "(@Converting the tree with all Octree = " << counter.tacAndElapsed() << "s). Group size is " << groupSize << "." << std::endl;

    counter.tic();
    GroupOctreeClass groupedTree3(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    std::cout << "Done  " << "(@Converting the tree with all Octree = " << counter.tacAndElapsed() << "s). Group size is " << groupSize << "." << std::endl;

    groupedTree2.printInfoBlocks();
    groupedTree3.printInfoBlocks();



    typedef FRotationKernel< CellClass, FP2PGroupParticleContainer<> , P>   KernelClass;
    typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, CellClass, KernelClass, typename GroupOctreeClass::ParticleGroupClass, FP2PGroupParticleContainer<> > GroupAlgorithm;

    KernelClass kernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
    GroupAlgorithm algo(&groupedTree2,&kernel);

    algo.execute();

    return 0;
}
