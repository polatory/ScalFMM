// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Core/FFmmAlgorithmPeriodic.hpp"

#include "../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../Src/Kernels/Spherical/FSphericalParticle.hpp"

#include "../Src/Files/FEwalLoader.hpp"
#include "../Src/Components/FSimpleLeaf.hpp"

/** Ewal particle is used in the gadget program
  * here we try to make the same simulation
  */
class EwalParticle : public FSphericalParticle {
public:
    // Type of particle
    enum Type{
        OW,
        HW,
        Undefined,
    };

private:
    Type type; //< current type
    int index; //< current index in array

public:
    // Basic constructor
    EwalParticle() : type(Undefined), index(-1) {
    }

    Type getType() const{
        return type;
    }

    void setType(const Type inType) {
        type = inType;
    }

    int getIndex() const{
        return index;
    }

    void setIndex( const int inIndex ){
        index = inIndex;
    }
};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef EwalParticle            ParticleClass;
    typedef FSphericalCell          CellClass;
    typedef FVector<ParticleClass>  ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >   KernelClass;

    typedef FFmmAlgorithmPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels      = FParameters::getValue(argc,argv,"-h", 4);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 2);
    const int DevP          = FParameters::getValue(argc,argv,"-P", 5);
    const int BoundaryDeep  = FParameters::getValue(argc,argv,"-bd", 5);
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/testEwal417.dt");
    FTic counter;

    // -----------------------------------------------------

    std::cout << "Opening : " << filename << "\n";
    FEwalLoader<ParticleClass> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    CellClass::Init(DevP);
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tWidth : " << loader.getBoxWidth() << " \t center x : " << loader.getCenterOfBox().getX()
              << " y : " << loader.getCenterOfBox().getY() << " z : " << loader.getCenterOfBox().getZ() << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    counter.tic();

    ParticleClass * const particles = new ParticleClass[loader.getNumberOfParticles()];

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(particles[idxPart]);
        // reset forces and insert in the tree
        ParticleClass part = particles[idxPart];
        part.setForces(0,0,0);
        part.setPotential(0);
        tree.insert(part);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels(DevP, NbLevels,loader.getBoxWidth(),loader.getCenterOfBox(), BoundaryDeep);
    FmmClass algo(&tree,&kernels,BoundaryDeep);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------
    {
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){
                std::cout << ">> index " << iter.data().getIndex() << " type " << iter.data().getType() << std::endl;
                std::cout << "x " << iter.data().getPosition().getX() << " y " << iter.data().getPosition().getY() << " z " << iter.data().getPosition().getZ() << std::endl;
                std::cout << "fx " << iter.data().getForces().getX() << " fy " << iter.data().getForces().getY() << " fz " << iter.data().getForces().getZ() << std::endl;
                std::cout << "physical value " << iter.data().getPhysicalValue() << " potential " << iter.data().getPotential() << std::endl;

                const ParticleClass& part = particles[iter.data().getIndex()];
                std::cout << "x " << part.getPosition().getX() << " y " << part.getPosition().getY() << " z " << part.getPosition().getZ() << std::endl;
                std::cout << "fx " <<part.getForces().getX() << " fy " << part.getForces().getY() << " fz " << part.getForces().getZ() << std::endl;
                std::cout << "physical value " << part.getPhysicalValue() << " potential " << part.getPotential() << std::endl;

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());
    }

    // -----------------------------------------------------

    delete[] particles;

    return 0;
}



