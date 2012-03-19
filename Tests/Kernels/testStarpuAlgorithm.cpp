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

// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_STARPU
// ================

#include <starpu.h>


#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/Components/FTestParticle.hpp"
#include "../../Src/Components/FTestCell.hpp"

#include "../../Src/Core/FFmmAlgorithmStarpu.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Components/FFmaParticle.hpp"
#include "../../Src/Extensions/FExtendForces.hpp"
#include "../../Src/Extensions/FExtendPotential.hpp"

#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Files/FFmaLoader.hpp"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass, class CellType>
KernelClass** FFmmAlgorithmStarpu<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass,CellType>::globalKernels = 0;

// export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
// Compile With openmp : g++ testFmbAlgorithm.cpp ../../Src/Utils/FDebug.cpp ../../Src/Utils/FTrace.cpp -lgomp -fopenmp -lstarpu -O2 -o testFmbAlgorithm.exe
//
// g++ -L../starpu/lib/ -I../starpu/include testFmbAlgorithmNoProc.cpp ../../Src/Utils/FDebug.cpp ../../Src/Utils/FTrace.cpp ../../Src/Utils/FMath.cpp ../../Src/Utils/F3DPosition.cpp -lgomp -fopenmp -lstarpu -O2 -o testFmbAlgorithm.exe

////////////////////////////////////////////////////////////////
// Define classes
////////////////////////////////////////////////////////////////

#define USE_TEST

template< class ParticleClass, class CellClass, class ContainerClass>
class TestKernels  : public FAbstractKernels<ParticleClass,CellClass,ContainerClass> {
public:
    /** Before upward */
    void P2M(CellClass* const pole, const ContainerClass* const particles) {
        // the pole represents all particles under
        *pole = (*pole) + particles->getSize();
    }

    /** During upward */
    void M2M(CellClass* const FRestrict pole, const CellClass *const FRestrict *const FRestrict child, const int /*level*/) {
        // A parent represents the sum of the child
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child[idx]){
                *pole = (*pole) + (*child[idx]);
            }
        }
    }

    /** Before Downward */
    void M2L(CellClass* const FRestrict pole, const CellClass* distantNeighbors[343], const int /*size*/, const int /*level*/) {
        // The pole is impacted by what represent other poles
        for(int idx = 0 ; idx < 343 ; ++idx){
            if(distantNeighbors[idx]){
                *pole = (*pole) + (*distantNeighbors[idx]);
            }
        }
    }

    /** During Downward */
    void L2L(const CellClass*const FRestrict local, CellClass* FRestrict *const FRestrict child, const int /*level*/) {
        // Each child is impacted by the father
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child[idx]){
                *child[idx] = (*local) + (*child[idx]);
            }
        }

    }

    /** After Downward */
    void L2P(const CellClass* const  local, ContainerClass*const particles){
        // The particles is impacted by the parent cell
        typename ContainerClass::BasicIterator iter(*particles);
        while( iter.hasNotFinished() ){
            iter.data().setDataDown(iter.data().getDataDown() + (*local));
            iter.gotoNext();
        }

    }


    /** After Downward */
    void P2P(const FTreeCoordinate& ,
                 ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
                 ContainerClass* const directNeighborsParticles[27], const int ){
        // Each particles targeted is impacted by the particles sources
        long long int inc = sources->getSize();
        if(targets == sources){
            inc -= 1;
        }
        for(int idx = 0 ; idx < 27 ; ++idx){
            if( directNeighborsParticles[idx] ){
                inc += directNeighborsParticles[idx]->getSize();
            }
        }

        typename ContainerClass::BasicIterator iter(*targets);
        while( iter.hasNotFinished() ){
            iter.data().setDataDown(iter.data().getDataDown() + inc);
            iter.gotoNext();
        }

    }

    /** After Downward */
    void P2PRemote(const FTreeCoordinate& ,
                 ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
                 ContainerClass* const directNeighborsParticles[27], const int ){

        // Each particles targeted is impacted by the particles sources
        long long int inc = 0;
        for(int idx = 0 ; idx < 27 ; ++idx){
            if( directNeighborsParticles[idx] ){
                inc += directNeighborsParticles[idx]->getSize();
            }
        }

        typename ContainerClass::BasicIterator iter(*targets);
        while( iter.hasNotFinished() ){
            iter.data().setDataDown(iter.data().getDataDown() + inc);
            iter.gotoNext();
        }

    }
};


class StarCell : public AbstractStarCell, public FSphericalCell {
public:
    void initHandle(){
        AbstractStarCell::handleUp.registerData(FSphericalCell::getMultipole(), sizeof(FComplexe) * FSphericalCell::GetPoleSize());
        AbstractStarCell::handleDown.registerData(FSphericalCell::getLocal(), sizeof(FComplexe) * FSphericalCell::GetLocalSize());
    }
};

class StarTestCell : public AbstractStarCell, public FTestCell {
public:
    void initHandle(){
        AbstractStarCell::handleUp.registerData( &dataUp, sizeof(dataUp));
        AbstractStarCell::handleDown.registerData( &dataDown, sizeof(dataDown));
    }
};

class TestParticle : public FTestParticle, public FExtendPhysicalValue{
};

////////////////////////////////////////////////////////////////
// Typedefs
////////////////////////////////////////////////////////////////
#ifdef USE_TEST
    typedef TestParticle             ParticleClass;
    typedef StarTestCell                 CellClass;
    typedef long long               CellType;
#else
    typedef FSphericalParticle        ParticleClass;
    typedef StarCell                  CellClass;
    typedef FComplexe                 CellType;
#endif

typedef StarVector<ParticleClass>      ContainerClass;

typedef FSimpleLeaf<ParticleClass, ContainerClass >                         LeafClass;
typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

#ifdef USE_TEST
typedef TestKernels<ParticleClass, CellType, DataVector<ParticleClass> >         KernelClass;
#else
typedef FSphericalKernel<ParticleClass, CellClass, DataVector<ParticleClass> >     KernelClass;
#endif

typedef FFmmAlgorithmStarpu<OctreeClass,ParticleClass,CellClass, ContainerClass,KernelClass,LeafClass,CellType>  AlgorithmClass;

////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test fmb algorithm.\n";
    //////////////////////////////////////////////////////////////
    const int DevP = FParameters::getValue(argc,argv,"-p", 8);
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");

    std::cout << "Opening : " << filename << "\n";

    FFmaLoader<ParticleClass> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------
#ifndef USE_TEST
    CellClass::Init(DevP);
#endif
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    loader.fillTree(tree);

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

#ifdef USE_TEST
    KernelClass kernel;
#else
    KernelClass kernel(DevP, NbLevels,loader.getBoxWidth(), loader.getCenterOfBox());
#endif
    AlgorithmClass algo( &tree, &kernel);
    std::cout << "There are " << starpu_worker_get_count() << " workers" << std::endl;
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;


    // Check result
#ifdef USE_TEST
    ValidateFMMAlgo<OctreeClass, ParticleClass, CellClass, ContainerClass, LeafClass>(&tree);
#else
    { // get sum forces&potential
        FReal potential = 0;
        F3DPosition forces;
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){
                potential += iter.data().getPotential() * iter.data().getPhysicalValue();
                forces += iter.data().getForces();

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());

        std::cout << "Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
        std::cout << "Potential = " << potential << std::endl;
    }
#endif

    return 0;
}
