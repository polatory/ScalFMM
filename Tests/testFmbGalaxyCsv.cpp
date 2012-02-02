#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

#include "../Src/Kernels/FElecForcesKernels.hpp"
#include "../Src/Kernels/FComputeCell.hpp"
#include "../Src/Fmb/FFmbComponents.hpp"

#include "../Src/Extensions/FExtendVelocity.hpp"

#include "../Src/Files/FTreeCsvSaver.hpp"
#include "../Src/Files/FFmaLoader.hpp"
#include "../Src/Arranger/FOctreeArranger.hpp"


class FmbVeloParticle : public FmbParticle, public FExtendVelocity {
};

template <class ParticleClass>
class GalaxyLoader : public FFmaLoader<ParticleClass> {
public:
    GalaxyLoader(const char* const filename) : FFmaLoader<ParticleClass>(filename) {
    }

    void fillParticle(ParticleClass& inParticle){
        FReal x,y,z,data, vx, vy, vz;
        this->file >> x >> y >> z >> data >> vx >> vy >> vz;
        inParticle.setPosition(x,y,z);
        inParticle.setPhysicalValue(data);
        inParticle.setVelocity(vx,vy,vz);
    }
};


template <class OctreeClass, class ContainerClass , class ParticleClass>
class MassSaver : public FTreeCsvSaver<OctreeClass,ContainerClass, ParticleClass> {
public:
    MassSaver(const char inBasefile[], const bool inIncludeHeader = false)
        : FTreeCsvSaver<OctreeClass,ContainerClass, ParticleClass> (inBasefile,inIncludeHeader) {
    }

    virtual FReal getValue(ParticleClass*const part){
        return part->getPhysicalValue();
    }
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FmbVeloParticle         ParticleClass;
    typedef FComputeCell            CellClass;
    typedef FVector<ParticleClass>  ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FElecForcesKernels<ParticleClass, CellClass, ContainerClass >   KernelClass;

    typedef FFmmAlgorithmThread<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test fmb algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 6);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const FReal DT          = FParameters::getValue(argc,argv,"-dt", FReal(0.1));
    const int DevP          = FParameters::getValue(argc,argv,"-p", 5);

    FComputeCell::Init(DevP);

    GalaxyLoader<ParticleClass> loader("../Data/galaxy.fma.tmp");

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    {
        ParticleClass particleToFill;
        particleToFill.setPhysicalValue(FReal(0.10));

        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(particleToFill);
            tree.insert(particleToFill);
        }
    }

    // -----------------------------------------------------

    KernelClass kernels( DevP, NbLevels, loader.getBoxWidth());
    FmmClass algo( &tree, &kernels);
    FOctreeArranger<OctreeClass, ContainerClass, ParticleClass> arranger(&tree);
    MassSaver<OctreeClass, ContainerClass, ParticleClass> saver("./out/test%d.csv");

    for(int idx = 0; idx < 100 ; ++idx){
        algo.execute();
        { // update velocity and position
            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            do{
                typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
                while( iter.hasNotFinished() ){
                    kernels.computeVelocity(&iter.data(), DT);
                    kernels.updatePosition(&iter.data(), DT);
                    iter.gotoNext();
                }
            } while(octreeIterator.moveRight());
        }
        // update tree and vtk
        arranger.rearrange(true);
        saver.exportTree(&tree);
    }

    // -----------------------------------------------------

    return 0;
}
