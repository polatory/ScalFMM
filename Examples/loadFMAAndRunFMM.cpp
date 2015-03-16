#include <fstream>
#include <memory>
#include <string>
#include <sys/ioctl.h>

#include "Utils/FMath.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"
#include "Files/FFmaGenericLoader.hpp"
#include "Core/FFmmAlgorithm.hpp"

#include "Containers/FOctree.hpp"
#include "Components/FBasicCell.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Components/FBasicParticleContainer.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "FChebBalanceSymKernel.hpp"
#include "loadFMAAndRunFMMArgs.hpp"
#include "CostZones.hpp"

typedef FCostCell                                       CellClass;
typedef FBasicParticleContainer<0>                      ContainerClass;
typedef FSimpleLeaf< ContainerClass >                   LeafClass;
typedef FOctree< CellClass, ContainerClass, LeafClass > OctreeClass;
typedef FInterpMatrixKernelR                            MatrixKernelClass;
typedef FChebBalanceSymKernel<CellClass, 
                              ContainerClass, 
                              MatrixKernelClass,
                              5,
                              OctreeClass>              KernelClass;

const FReal epsilon = 1e-4;

int main(int argc, char** argv)
{
    loadFMAAndRunFMMArgs args(argc, argv);

    FFmaGenericLoader loader(args.inFileName().c_str());
    OctreeClass tree(args.treeHeight(),
                     args.subTreeHeight(),
                     loader.getBoxWidth(),
                     loader.getCenterOfBox());

    FReal  physicalValue;
    FPoint particlePosition;
    // insertion
    for ( int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart ) {
        loader.fillParticle(&particlePosition, &physicalValue);
        tree.insert(particlePosition);
    }
    
    KernelClass kernel(&tree, epsilon);
    FFmmAlgorithm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >
        algo(&tree, &kernel);

    algo.execute();

    kernel.printResults(std::cout);

    OctreeClass::Iterator it(&tree);

    CostZones<OctreeClass, CellClass> costzones(&tree, args.zoneCount());
    costzones.run();
    auto zones = costzones.getZones();

    // GCC versions before 5.0 have not implemented move constructors to streams
    std::vector<std::unique_ptr<std::ofstream>> outfiles;
    for ( int i = 0; i < args.treeHeight(); i++ ) {
        std::unique_ptr<std::ofstream> out(
            new std::ofstream( args.outFileName()
                               + "_"
                               + std::to_string(args.zoneCount())
                               + "z"
                               + "."
                               + std::to_string(i)
                               + args.outFileExt()));
        *out << "x,y,z,zone" << std::endl;
        outfiles.push_back(std::move(out));
    }
    
    int zoneIdx = 0;
    for ( auto zone : zones) {
        for ( auto cell : zone) {
            *(outfiles[cell.first]) << cell.second->getCoordinate().getX() << ",";
            *(outfiles[cell.first]) << cell.second->getCoordinate().getY() << ",";
            *(outfiles[cell.first]) << cell.second->getCoordinate().getZ() << ",";
            *(outfiles[cell.first]) << zoneIdx << "," << cell.first << std::endl;
        }
        zoneIdx++;
    }

    return EXIT_SUCCESS;
}


