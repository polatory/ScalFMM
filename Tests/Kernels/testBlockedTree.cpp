
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



int main(int argc, char* argv[]){
    static const int P = 9;
    typedef FRotationCell<P>               CellClass;
    typedef FP2PParticleContainer<>          ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
    //typedef FRotationKernel< CellClass, ContainerClass , P>   KernelClass;
    typedef FGroupTree< CellClass, 4, FReal>  GroupOctreeClass;

    FTic counter;
    const int NbLevels      = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.bin.fma.double");

    FFmaBinLoader loader(filename);
    FAssertLF(loader.isOpen());

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint particlePosition;
        FReal physicalValue;
        loader.fillParticle(&particlePosition,&physicalValue);
        tree.insert(particlePosition, physicalValue );
    }

    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.tacAndElapsed() << "s)." << std::endl;

    const int groupSize      = FParameters::getValue(argc,argv,"-bs", 250);

    counter.tic();
    GroupOctreeClass groupedTree(NbLevels, groupSize, &tree);
    std::cout << "Done  " << "(@Converting the tree with leafs only = " << counter.tacAndElapsed() << "s). Group size is " << groupSize << "." << std::endl;

    counter.tic();
    //GroupOctreeClass groupedTree2(NbLevels, groupSize, &tree);
    std::cout << "Done  " << "(@Converting the tree with all Octree = " << counter.tacAndElapsed() << "s). Group size is " << groupSize << "." << std::endl;


    groupedTree.printInfoBlocks();
    //groupedTree2.printInfoBlocks();


    return 0;
}
