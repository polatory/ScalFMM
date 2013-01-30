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

// ==== CMAKE =====
// @FUSE_MPI
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMpi.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Extensions/FExtendVelocity.hpp"

#include "../../Src/Files/FTreeMpiCsvSaver.hpp"
#include "../../Src/Files/FFmaLoader.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"

class FmmVeloParticle : public FSphericalParticle, public FExtendVelocity {
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
        inParticle.setPhysicalValue(FReal(0.10));
    }
};


template <class OctreeClass, class ContainerClass , class ParticleClass>
class MassSaver : public FTreeMpiCsvSaver<OctreeClass,ContainerClass, ParticleClass> {
public:
    MassSaver(const char inBasefile[], const FMpi::FComm& inComm, const bool inIncludeHeader = false)
        : FTreeMpiCsvSaver<OctreeClass,ContainerClass, ParticleClass>(inBasefile, inComm, inIncludeHeader) {
    }

    virtual FReal getValue(ParticleClass*const part){
        return part->getPhysicalValue();
    }
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FmmVeloParticle         ParticleClass;
    typedef FBasicCell              CellClass;
    typedef FVector<ParticleClass>  ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////
    FMpi app( argc, argv);

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 6);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);

    GalaxyLoader<ParticleClass> loader(FParameters::getStr(argc,argv,"-f", "../Data/galaxy.fma.tmp"));

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    {
        ParticleClass particleToFill;
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(particleToFill);
            if( (idxPart+1) % (app.global().processId()+1) == 0) tree.insert(particleToFill);
        }
    }

    // -----------------------------------------------------

    {
        MassSaver<OctreeClass, ContainerClass, ParticleClass> saver("./out/test%d.csv", app.global() , false);
        saver.exportTree(&tree);
    }

    // -----------------------------------------------------

    {
        MassSaver<OctreeClass, ContainerClass, ParticleClass> saver("./out/htest%d.csv", app.global() , true);
        saver.exportTree(&tree);
    }

    return 0;
}
