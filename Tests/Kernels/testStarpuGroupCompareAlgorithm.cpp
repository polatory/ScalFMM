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
// @FUSE_STARPU
// ================

#include <starpu.h>


#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithmStarpuGroup.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Files/FFmaLoader.hpp"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>



static const FReal Epsilon = FReal(0.0005);

///////////////////////////////////////////////////////
// to test equality between good and potentialy bad solution
///////////////////////////////////////////////////////
/** To compare data */
template <class CellClass, class CellClass2>
bool isEqualPole(const CellClass& me, const CellClass2& other, FReal*const cumul){
    FMath::FAccurater accurate;
    for(int idx = 0; idx < CellClass::GetPoleSize(); ++idx){
        accurate.add(me.getMultipole()[idx].getImag(),other.getMultipole()[idx].getImag());
        accurate.add(me.getMultipole()[idx].getReal(),other.getMultipole()[idx].getReal());
    }
    *cumul = accurate.getInfNorm()+ accurate.getL2Norm();
    return accurate.getInfNorm() < Epsilon && accurate.getL2Norm() < Epsilon;//FMath::LookEqual(cumul,FReal(0.0));
}

/** To compare data */
template <class CellClass, class CellClass2>
bool isEqualLocal(const CellClass& me, const CellClass2& other, FReal*const cumul){
    FMath::FAccurater accurate;
    for(int idx = 0; idx < FSphericalCell::GetLocalSize(); ++idx){
        accurate.add(me.getLocal()[idx].getImag(),other.getLocal()[idx].getImag());
        accurate.add(me.getLocal()[idx].getReal(),other.getLocal()[idx].getReal());
    }
    *cumul = accurate.getInfNorm()+ accurate.getL2Norm();
    return accurate.getInfNorm() < Epsilon && accurate.getL2Norm() < Epsilon;//FMath::LookEqual(cumul,FReal(0.0));
}


template<class OctreeClass, class OctreeClass2, class ContainerClass, class ContainerClass2, class CellClass, class CellClass2>
void ValidateFMMAlgoProc(OctreeClass2* const badTree, OctreeClass* const valideTree){
    std::cout << "Check Result\n";
    {
        const int OctreeHeight = valideTree->getHeight();
        typename OctreeClass2::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level > 1 ; --level){
            while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
                octreeIteratorValide.moveRight();
            }

            do {
                if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << "Error index are not equal!" << std::endl;
                }
                else{
                    FReal cumul;
                    if( !isEqualPole(*octreeIterator.getCurrentCell(),*octreeIteratorValide.getCurrentCell(),&cumul) ){
                        std::cout << "Pole Data are different. Cumul " << cumul << " at level " << level << " index is " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }
                    if( !isEqualLocal(*octreeIterator.getCurrentCell(),*octreeIteratorValide.getCurrentCell(),&cumul) ){
                        std::cout << "Local Data are different. Cumul " << cumul << " at level " << level << " index is " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }
                }

            } while(octreeIterator.moveRight() && octreeIteratorValide.moveRight());

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            octreeIteratorValide.moveUp();
            octreeIteratorValide.gotoLeft();
        }
    }
    {
        // Check that each particle has been summed with all other
        typename OctreeClass2::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            octreeIteratorValide.moveRight();
        }

        do {
            if( octreeIterator.getCurrentListSrc()->getSize() != octreeIteratorValide.getCurrentListSrc()->getSize()){
                std::cout << " Particules numbers is different " << std::endl;
            }
            if( octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                std::cout << " Index are differents " << std::endl;
            }

            typename ContainerClass2::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            while( iter.hasNotFinished() ){

                typename ContainerClass::BasicIterator iterValide(*octreeIteratorValide.getCurrentListTargets());
                while( iterValide.hasNotFinished() ){
                    if( FMath::LookEqual(iterValide.data().getPosition().getX(),iter.data().getPosition().getX()) &&
                        FMath::LookEqual(iterValide.data().getPosition().getY(),iter.data().getPosition().getY()) &&
                        FMath::LookEqual(iterValide.data().getPosition().getZ(),iter.data().getPosition().getZ()) ){
                        break;
                    }
                    iterValide.gotoNext();
                }

                if( iterValide.hasNotFinished() ){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    bool error = false;
                    if( FMath::RelatifDiff(iterValide.data().getPotential() , iter.data().getPotential())  > Epsilon ){
                        std::cout << " Potential error : " << iterValide.data().getPotential()  << " " << iter.data().getPotential() << "\n";
                        error = true;
                    }
                    if( FMath::RelatifDiff(iterValide.data().getForces().getX(),iter.data().getForces().getX()) > Epsilon
                            || FMath::RelatifDiff(iterValide.data().getForces().getY(),iter.data().getForces().getY()) > Epsilon
                            || FMath::RelatifDiff(iterValide.data().getForces().getZ(),iter.data().getForces().getZ()) > Epsilon){
                        std::cout << " Forces error : x " << iterValide.data().getForces().getX() << " " << iter.data().getForces().getX()
                                  << " y " << iterValide.data().getForces().getY()  << " " << iter.data().getForces().getY()
                                  << " z " << iterValide.data().getForces().getZ()  << " " << iter.data().getForces().getZ() << "\n";
                        error = true;
                    }
                    if( error ){
                        std::cout << "At position " << iterValide.data().getPosition() << " == " << iter.data().getPosition() << std::endl;
                    }
                }
                else{
                    std::cout << "Particle not found " << iter.data().getPosition() << std::endl;
                }
                iter.gotoNext();
            }

        } while(octreeIterator.moveRight() && octreeIteratorValide.moveRight());
    }

    std::cout << "Done\n";
}


////////////////////////////////////////////////////////////////
// Typedefs
////////////////////////////////////////////////////////////////

class SphericalCell : public FSphericalCell {
public:
    void intialCopy(const SphericalCell*const other){
        FMemUtils::copyall(multipole_exp, other->multipole_exp, PoleSize);
        FMemUtils::copyall(local_exp, other->local_exp, LocalSize);
    }
    void copyUp(const SphericalCell*const other){
        FMemUtils::copyall(multipole_exp, other->multipole_exp, PoleSize);
    }
    void restoreCopy(const SphericalCell*const other){
        FMemUtils::copyall(multipole_exp, other->multipole_exp, PoleSize);
        FMemUtils::copyall(local_exp, other->local_exp, LocalSize);
    }
};

template <class ParticleClass>
class Container : public FVector<ParticleClass> {
public:
    void reduce(const Container*const other){
        for( int idx = 0 ; idx < FVector<ParticleClass>::getSize() ; ++idx){
            FVector<ParticleClass>::data()[idx].incForces(other->FVector<ParticleClass>::data()[idx].getForces());
            FVector<ParticleClass>::data()[idx].incPotential(other->FVector<ParticleClass>::data()[idx].getPotential());
        }
    }
};


typedef FSphericalParticle             ParticleClass;
typedef Container<ParticleClass>       ContainerClass;

typedef SphericalCell CellClass;


typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

typedef FFmmAlgorithmStarpuGroup<OctreeClass, ParticleClass, CellClass, ContainerClass,KernelClass,LeafClass>  AlgorithmClass;


typedef FSphericalParticle             ParticleClass2;
typedef FSphericalCell                 CellClass2;
typedef FVector<ParticleClass2>         ContainerClass2;

typedef FSimpleLeaf<ParticleClass2, ContainerClass2 >                     LeafClass2;
typedef FOctree<ParticleClass2, CellClass2, ContainerClass2 , LeafClass2 >  OctreeClass2;
typedef FSphericalKernel<ParticleClass2, CellClass2, ContainerClass2 >     KernelClass2;

typedef FFmmAlgorithmThread<OctreeClass2, ParticleClass2, CellClass2, ContainerClass2, KernelClass2, LeafClass2 > FmmClass2;
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
    const int BlockSize = FParameters::getValue(argc,argv,"-bs", 250);
    const int NbThread = FParameters::getValue(argc,argv,"-t", -1);
    const bool usePerfModel = (FParameters::findParameter(argc, argv, "-perf") != FParameters::NotFound);
    const bool useReductionPart = (FParameters::findParameter(argc, argv, "-reducepart") != FParameters::NotFound);
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");

    std::cout << "Opening : " << filename << "\n";

    FFmaLoader<ParticleClass> loader(filename);
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

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    KernelClass kernel(DevP, NbLevels,loader.getBoxWidth(), loader.getCenterOfBox());
    AlgorithmClass algo( &tree, &kernel, BlockSize, usePerfModel, useReductionPart);

    // -----------------------------------------------------

    std::cout << "Build gouped tree..." << std::endl;
    counter.tic();
    algo.buildGroups(NbThread);
    counter.tac();
    std::cout << "Done  in " << counter.elapsed() << "s." << std::endl;

    // -----------------------------------------------------

    std::cout << "Execute Fmm..." << std::endl;

    counter.tic();
    algo.execute();
    counter.tac();

    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Release gouped tree..." << std::endl;
    counter.tic();
    algo.releaseGroups();
    counter.tac();
    std::cout << "Done in " << counter.elapsed() << "s." << std::endl;

    // -----------------------------------------------------

    FFmaLoader<ParticleClass2> loader2(filename);
    if(!loader2.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------
    CellClass2::Init(DevP);
    OctreeClass2 tree2(NbLevels, SizeSubLevels, loader2.getBoxWidth(), loader2.getCenterOfBox());

    // -----------------------------------------------------

    loader2.fillTree(tree2);

    // -----------------------------------------------------

    KernelClass2 kernels2(DevP, NbLevels,loader2.getBoxWidth(), loader2.getCenterOfBox());
    FmmClass2 algo2(&tree2,&kernels2);

    counter.tic();
    algo2.execute();
    counter.tac();

    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    ValidateFMMAlgoProc<OctreeClass2, OctreeClass, ContainerClass2, ContainerClass, CellClass2, CellClass>(&tree, &tree2);

    return 0;
}
