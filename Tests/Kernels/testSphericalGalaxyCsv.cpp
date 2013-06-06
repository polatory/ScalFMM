// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#include "../../Src/Files/FTreeCsvSaver.hpp"
#include "../../Src/Files/FFmaLoader.hpp"
#include "../../Src/Arranger/FOctreeArranger.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"


class VelocityContainer : public FP2PParticleContainer {
    typedef FP2PParticleContainer Parent;

    FVector<FPoint> velocities;

public:
    template<typename... Args>
    void push(const FPoint& inParticlePosition, const FPoint& velocity, Args... args){
        Parent::push(inParticlePosition, args... );
        velocities.push(velocity);
    }

    const FVector<FPoint>& getVelocities() const{
        return velocities;
    }

    FVector<FPoint>& getVelocities() {
        return velocities;
    }

    void fillToCsv(const int partIdx, FReal values[4]) const {
        values[0] = Parent::getPositions()[0][partIdx];
        values[1] = Parent::getPositions()[1][partIdx];
        values[2] = Parent::getPositions()[2][partIdx];
        values[3] = Parent::getPotentials()[partIdx];
    }
};


class GalaxyLoader : public FFmaLoader {
public:
    GalaxyLoader(const char* const filename) : FFmaLoader(filename) {
    }

    void fillParticle(FPoint* position, FReal* physivalValue, FPoint* velocity){
        FReal x,y,z,data, vx, vy, vz;
        this->file >> x >> y >> z >> data >> vx >> vy >> vz;
        position->setPosition(x,y,z);
        *physivalValue = (data);
        velocity->setPosition(vx,vy,vz);
    }
};

struct TestParticle{
    FPoint position;
    FReal physicalValue;
    FReal forces[3];
    FReal potential;
    FPoint velocity;
    const FPoint& getPosition(){
        return position;
    }
};

template <class ParticleClass>
class Converter {
public:
    template <class ContainerClass>
    static ParticleClass GetParticle(ContainerClass* containers, const int idxExtract){
        const FReal*const positionsX = containers->getPositions()[0];
        const FReal*const positionsY = containers->getPositions()[1];
        const FReal*const positionsZ = containers->getPositions()[2];
        const FReal*const forcesX = containers->getForcesX();
        const FReal*const forcesY = containers->getForcesY();
        const FReal*const forcesZ = containers->getForcesZ();
        const FReal*const physicalValues = containers->getPhysicalValues();
        const FReal*const potentials = containers->getPotentials();
        FVector<FPoint> velocites = containers->getVelocities();

        TestParticle part;
        part.position.setPosition( positionsX[idxExtract],positionsY[idxExtract],positionsZ[idxExtract]);
        part.physicalValue = physicalValues[idxExtract];
        part.forces[0] = forcesX[idxExtract];
        part.forces[1] = forcesY[idxExtract];
        part.forces[2] = forcesZ[idxExtract];
        part.potential = potentials[idxExtract];
        part.velocity  = velocites[idxExtract];

        return part;
    }

    template <class OctreeClass>
    static void Insert(OctreeClass* tree, const ParticleClass& part){
        tree->insert(part.position , part.velocity, part.physicalValue, part.forces[0],
                part.forces[1],part.forces[2],part.potential);
    }
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FSphericalCell          CellClass;
    typedef VelocityContainer  ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel< CellClass, ContainerClass >   KernelClass;

    typedef FFmmAlgorithmThread<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 6);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const FReal DT          = FParameters::getValue(argc,argv,"-dt", FReal(0.1));
    const int DevP          = FParameters::getValue(argc,argv,"-p", 5);

    FSphericalCell::Init(DevP);

    GalaxyLoader loader(FParameters::getStr(argc,argv,"-f", "../Data/galaxy.fma.tmp"));

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    {
        FPoint position, velocity;
        FReal physicalValue;

        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&position, &physicalValue, &velocity);
            tree.insert(position, velocity, physicalValue);
        }
    }

    // -----------------------------------------------------

    KernelClass kernels( DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
    FmmClass algo( &tree, &kernels);
    FOctreeArranger<OctreeClass, ContainerClass, TestParticle, Converter<TestParticle> > arranger(&tree);
    FTreeCsvSaver<OctreeClass, ContainerClass> saver("./out/test%d.csv");

    for(int idx = 0; idx < 100 ; ++idx){
        algo.execute();
        { // update velocity and position
            OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            do{
                kernels.computeVelocity(octreeIterator.getCurrentListTargets(), DT);
                kernels.updatePosition(octreeIterator.getCurrentListTargets(), DT);
            } while(octreeIterator.moveRight());
        }
        // update tree and vtk
        arranger.rearrange(AllDirs);
        saver.exportTree(&tree);
    }

    // -----------------------------------------------------

    return 0;
}
