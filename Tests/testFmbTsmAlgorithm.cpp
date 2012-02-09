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

#include "../Src/Components/FTypedLeaf.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Core/FFmmAlgorithmTsm.hpp"

#include "../Src/Kernels/FSphericalKernel.hpp"
#include "../Src/Kernels/FSphericalCell.hpp"
#include "../Src/Kernels/FSphericalParticle.hpp"

#include "../Src/Files/FFmaTsmLoader.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FTypedSphericalParticle        ParticleClass;
    typedef FTypedSphericalCell            CellClass;
    typedef FVector<ParticleClass>         ContainerClass;

    typedef FTypedLeaf<ParticleClass, ContainerClass >                      LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

    typedef FFmmAlgorithmTsm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Fmb on a Tsm system.\n";
    //////////////////////////////////////////////////////////////
    const int DevP = FParameters::getValue(argc,argv,"-p", 8);
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;

    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.tsm.fma");
    std::cout << "Opening : " << filename << "\n";

    FFmaTsmLoader<ParticleClass> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------
    CellClass::Init(DevP);
    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    loader.fillTree(tree);

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels(DevP, NbLevels, loader.getBoxWidth());

    FmmClass algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    { // get sum forces&potential
        FReal potential = 0;
        F3DPosition forces;
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            FVector<ParticleClass>::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){
                potential += iter.data().getPotential() * iter.data().getPhysicalValue();
                forces += iter.data().getForces();

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());

        std::cout << "Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
        std::cout << "Potential = " << potential << std::endl;
    }


    // -----------------------------------------------------

    return 0;
}



