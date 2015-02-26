// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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
#include "../../Src/Files/FFmaGenericLoader.hpp"
#include "../../Src/Arranger/FOctreeArranger.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FParticleType.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"
#include "../../Src/Arranger/FAbstractMover.hpp"


class VelocityContainer : public FP2PParticleContainer<> {
    typedef FP2PParticleContainer<> Parent;
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


class GalaxyLoader : public FFmaGenericLoader {
public:
    GalaxyLoader(const char* const filename) : FFmaGenericLoader(filename) {
    }

    void fillParticle(FPoint* position, FReal* physivalValue, FPoint* velocity){
        FReal x,y,z,data, vx, vy, vz;
        *(this->file) >> x >> y >> z >> data >> vx >> vy >> vz;
        position->setPosition(x,y,z);
        *physivalValue = (data);
        velocity->setPosition(vx,vy,vz);
    }
};

template<class OctreeClass>
class GalaxyMover : public FAbstractMover<OctreeClass, VelocityContainer>{
private:
    VelocityContainer toStoreRemovedParts;

public:
    GalaxyMover() {
    }

    virtual ~GalaxyMover(){
    }

    /** To get the position of the particle at idx idxPart in leaf lf */
    void getParticlePosition(VelocityContainer* lf, const int idxPart, FPoint* particlePos){
        (*particlePos) = FPoint(lf->getPositions()[0][idxPart],lf->getPositions()[1][idxPart],lf->getPositions()[2][idxPart]);
    }

    /** Remove a particle but keep it to reinsert it later*/
    void removeFromLeafAndKeep(VelocityContainer* lf, const FPoint& particlePos, const int idxPart,FParticleType /*type*/){
        std::array<typename VelocityContainer::AttributesClass, VelocityContainer::NbAttributes> particleValues;
        for(int idxAttr = 0 ; idxAttr < VelocityContainer::NbAttributes ; ++idxAttr){
            particleValues[idxAttr] = lf->getAttribute(idxAttr)[idxPart];
        }

        toStoreRemovedParts.push(particlePos,lf->getVelocities()[idxPart],particleValues);

        lf->getVelocities().removeOne(idxPart);
        lf->removeParticles(&idxPart,1);
    }

    /** Reinsert the previously saved particles */
    void insertAllParticles(OctreeClass* tree){
        std::array<typename VelocityContainer::AttributesClass, VelocityContainer::NbAttributes> particleValues;

        for(int idxToInsert = 0; idxToInsert<toStoreRemovedParts.getNbParticles() ; ++idxToInsert){
            for(int idxAttr = 0 ; idxAttr < VelocityContainer::NbAttributes ; ++idxAttr){
                particleValues[idxAttr] = toStoreRemovedParts.getAttribute(idxAttr)[idxToInsert];
            }
            const FPoint particlePos(toStoreRemovedParts.getPositions()[0][idxToInsert],
                                     toStoreRemovedParts.getPositions()[1][idxToInsert],
                                     toStoreRemovedParts.getPositions()[2][idxToInsert]);

            tree->insert(particlePos, toStoreRemovedParts.getVelocities()[idxToInsert], particleValues);
        }

        toStoreRemovedParts.clear();
        toStoreRemovedParts.getVelocities().clear();
    }
};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Run a Spherical Harmonic (Old Implementation) FMM kernel with several time step.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::SHDevelopment,
                         FParameterDefinitions::DeltaT, FParameterDefinitions::OutputFile);

    typedef FSphericalCell          CellClass;
    typedef VelocityContainer  ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel< CellClass, ContainerClass >   KernelClass;

    typedef FFmmAlgorithmThread<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

    typedef GalaxyMover<OctreeClass> MoverClass;
    typedef FOctreeArranger<OctreeClass, ContainerClass, MoverClass> ArrangerClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 6);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const FReal DT          = FParameters::getValue(argc,argv,FParameterDefinitions::DeltaT.options, FReal(0.1));
    const int DevP          = FParameters::getValue(argc,argv,FParameterDefinitions::SHDevelopment.options, 5);

    FSphericalCell::Init(DevP);

    GalaxyLoader loader(FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/galaxy.fma.tmp"));

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
    ArrangerClass arranger(&tree);
    FTreeCsvSaver<OctreeClass, ContainerClass> saver(FParameters::getStr(argc,argv,FParameterDefinitions::OutputFile.options, "/tmp/test%d.csv"));

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
        arranger.rearrange();
        saver.exportTree(&tree);
    }

    // -----------------------------------------------------

    return 0;
}
