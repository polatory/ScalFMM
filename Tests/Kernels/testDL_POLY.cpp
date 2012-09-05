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
#include <iomanip>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

//#include "../../Src/Files/FFmaLoader.hpp"
#include "../../Src/Files/FEwalLoader.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is little or longer
  * related that each other
  */

class DLPolyParticle : public FSphericalParticle {
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
    DLPolyParticle() : type(Undefined), index(-1) {
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

//
// ------------------------------------------------------------------
//
// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef DLPolyParticle             ParticleClass;
    typedef FSphericalCell                 CellClass;
    typedef FVector<ParticleClass>         ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >     KernelClass;

    typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    std::cout << ">> You can pass -sequential or -task (thread by default).\n";
    //////////////////////////////////////////////////////////////
    const int DevP = FParameters::getValue(argc,argv,"-p", 8);
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;

    //
    const FReal coeff_MD= FReal( 138935.4835 / 418.4 );
    //

    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");

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
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    loader.fillTree(tree);
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        ParticleClass particle;
        loader.fillParticle(particle);
        // reset forces and insert in the tree
        particle.setIndex(idxPart);
        particle.setForces(0,0,0);
        particle.setPotential(0);
        // insert in tree
        tree.insert(particle);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Create kernel ..." << std::endl;
    counter.tic();

    KernelClass kernels(DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    counter.tac();
    std::cout << "Done  " << " in " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;

    FmmClass algo(&tree,&kernels);
    counter.tic();
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    { // get sum forces&potential
        FReal  potential = 0;
        FPoint forces;
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){
                potential += iter.data().getPotential() * iter.data().getPhysicalValue();
                forces += iter.data().getForces();
                std::cout << " " << iter.data().getIndex()+1 << " \t "<<
                             std::setprecision(5)<< iter.data().getPosition().getX() << "  \t" <<
                             iter.data().getPosition().getY() << "  \t" <<
                             iter.data().getPosition().getZ() << "   Forces: \t"<<
                             std::setprecision(8) << iter.data().getForces().getX()*coeff_MD << "  \t " <<
                             iter.data().getForces().getY()*coeff_MD << "  \t " <<
                             iter.data().getForces().getZ()*coeff_MD << std::endl;
		
                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());

        std::cout << "Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
        std::cout << "Potential = " << potential*coeff_MD << std::endl;
    }
    std::cout << "Constante DL_POLY: " << coeff_MD << std::endl;
    // -----------------------------------------------------


    return 0;
}



