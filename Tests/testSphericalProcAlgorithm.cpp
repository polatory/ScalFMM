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

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FMpi.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Kernels/FSphericalKernel.hpp"
#include "../Src/Kernels/FSphericalCell.hpp"
#include "../Src/Kernels/FSphericalParticle.hpp"

#include "../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Files/FMpiFmaLoader.hpp"
#include "../Src/Files/FMpiTreeBuilder.hpp"
#include "../Src/Files/FFmaBinLoader.hpp"

#include <iostream>

#include <cstdio>
#include <cstdlib>

// Uncoment to validate the FMM
//#define VALIDATE_FMM

/** This program show an example of use of
  * the fmm basic algo it also check that eachh particles is little or longer
  * related that each other
  */


///////////////////////////////////////////////////////
// to test equality between good and potentialy bad solution
///////////////////////////////////////////////////////
/** To compare data */
bool isEqualPole(const FSphericalCell& me, const FSphericalCell& other, FReal*const cumul){
    //return memcmp(multipole_exp, other.multipole_exp, sizeof(FComplexe)*MultipoleSize) == 0 &&
    //        memcmp(local_exp, other.local_exp, sizeof(FComplexe)*MultipoleSize) == 0;
    *cumul = 0.0;
    for(int idx = 0; idx < FSphericalCell::GetPoleSize(); ++idx){
        *cumul += FMath::Abs( me.getMultipole()[idx].getImag() - other.getMultipole()[idx].getImag() );
        *cumul += FMath::Abs( me.getMultipole()[idx].getReal() - other.getMultipole()[idx].getReal() );
    }

    return *cumul < 0.0001;//FMath::LookEqual(cumul,FReal(0.0));
}

/** To compare data */
bool isEqualLocal(const FSphericalCell& me, const FSphericalCell& other, FReal*const cumul){
    //return memcmp(multipole_exp, other.multipole_exp, sizeof(FComplexe)*MultipoleSize) == 0 &&
    //        memcmp(local_exp, other.local_exp, sizeof(FComplexe)*MultipoleSize) == 0;
    *cumul = 0.0;
    for(int idx = 0; idx < FSphericalCell::GetLocalSize(); ++idx){
        *cumul += FMath::Abs( me.getLocal()[idx].getImag() - other.getLocal()[idx].getImag() );
        *cumul += FMath::Abs( me.getLocal()[idx].getReal() - other.getLocal()[idx].getReal() );
    }

    return *cumul < 0.0001;//FMath::LookEqual(cumul,FReal(0.0));
}


#ifdef VALIDATE_FMM
template<class OctreeClass, class ContainerClass>
void ValidateFMMAlgoProc(OctreeClass* const badTree,
                         OctreeClass* const valideTree){
    std::cout << "Check Result\n";
    {
        const int OctreeHeight = valideTree->getHeight();
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level >= 1 ; --level){
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
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            octreeIteratorValide.moveRight();
        }

        do {
            typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            typename ContainerClass::BasicIterator iterValide(*octreeIteratorValide.getCurrentListTargets());

            if( octreeIterator.getCurrentListSrc()->getSize() != octreeIteratorValide.getCurrentListSrc()->getSize()){
                std::cout << " Particules numbers is different " << std::endl;
            }
            if( octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                std::cout << " Index are differents " << std::endl;
            }

            while( iter.hasNotFinished() && iterValide.hasNotFinished() ){
                // If a particles has been impacted by less than NbPart - 1 (the current particle)
                // there is a problem

                if( !FMath::LookEqual(iterValide.data().getPotential() , iter.data().getPotential()) ){
                    std::cout << " Potential error : " << iterValide.data().getPotential()  << " " << iter.data().getPotential() << "\n";
                }
                if( !FMath::LookEqual(iterValide.data().getForces().getX(),iter.data().getForces().getX())
                        || !FMath::LookEqual(iterValide.data().getForces().getY(),iter.data().getForces().getY())
                        || !FMath::LookEqual(iterValide.data().getForces().getZ(),iter.data().getForces().getZ()) ){
                    /*std::cout << idx << " Forces error : " << (iterValide.data().getForces().getX() - iter.data().getForces().getX())
                              << " " << (iterValide.data().getForces().getY() - iter.data().getForces().getY())
                              << " " << (iterValide.data().getForces().getZ() - iter.data().getForces().getZ()) << "\n";*/
                    std::cout << " Forces error : x " << iterValide.data().getForces().getX() << " " << iter.data().getForces().getX()
                              << " y " << iterValide.data().getForces().getY()  << " " << iter.data().getForces().getY()
                              << " z " << iterValide.data().getForces().getZ()  << " " << iter.data().getForces().getZ() << "\n";
                }
                iter.gotoNext();
                iterValide.gotoNext();
            }

        } while(octreeIterator.moveRight() && octreeIteratorValide.moveRight());
    }

    std::cout << "Done\n";
}
#endif


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FSphericalParticle     ParticleClass;
    typedef FSphericalCell         CellClass;
    typedef FVector<ParticleClass>         ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

    typedef FFmmAlgorithmThreadProc<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    typedef FFmmAlgorithmThread<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassNoProc;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app( argc, argv);

    const int DevP = FParameters::getValue(argc,argv,"-p", 8);
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.bin.fma");

    std::cout << "Opening : " << filename << "\n";

    FMpiFmaLoader<ParticleClass> loader(filename, app.global());
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

    if( app.global().processCount() != 1){
        //////////////////////////////////////////////////////////////////////////////////
        // Build tree from mpi loader
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Build Tree ..." << std::endl;
        counter.tic();

        FMpiTreeBuilder<ParticleClass>::LoaderToTree(app.global(), loader, tree);

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
    }
    else{
        ParticleClass partToInsert;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(partToInsert);
            tree.insert(partToInsert);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels(DevP, NbLevels,loader.getBoxWidth(), loader.getCenterOfBox());
    FmmClass algo(app.global(),&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    { // get sum forces&potential
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Sum Result" , __FILE__ , __LINE__) );

        FReal potential = 0;
        F3DPosition forces;

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
					ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished()){
                potential += iter.data().getPotential() * iter.data().getPhysicalValue();
                forces += iter.data().getForces();

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());

        std::cout << "My potential is " << potential << std::endl;

        potential = app.global().reduceSum(potential);
        forces.setX(app.global().reduceSum(forces.getX()));
        forces.setY(app.global().reduceSum(forces.getY()));
        forces.setZ(app.global().reduceSum(forces.getZ()));


        if(app.global().processId() == 0){
            std::cout << "Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
            std::cout << "Potential Sum = " << potential << std::endl;
        }
    }

#ifdef VALIDATE_FMM
    {
        OctreeClass treeValide(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
        {
            FFmaBinLoader<ParticleClass> loaderSeq(filename);
            ParticleClass partToInsert;
            for(FSize idxPart = 0 ; idxPart < loaderSeq.getNumberOfParticles() ; ++idxPart){
                loaderSeq.fillParticle(partToInsert);
                treeValide.insert(partToInsert);
            }
        }

        std::cout << "Working on particles ..." << std::endl;
        counter.tic();
        FmmClassNoProc algoValide(&treeValide,&kernels);
        algoValide.execute();
        counter.tac();
        std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

        FReal potentialValide = 0;
        F3DPosition forcesValide;

        OctreeClass::Iterator octreeIteratorValide(&treeValide);
        octreeIteratorValide.gotoBottomLeft();

        do{
            ContainerClass::ConstBasicIterator iterValide(*octreeIteratorValide.getCurrentListTargets());
            while( iterValide.hasNotFinished()){
                potentialValide += iterValide.data().getPotential() * iterValide.data().getPhysicalValue();
                forcesValide += iterValide.data().getForces();

                iterValide.gotoNext();
            }

        } while(octreeIteratorValide.moveRight());

        std::cout << "Valide Foces Sum  x = " << forcesValide.getX() << " y = " << forcesValide.getY() << " z = " << forcesValide.getZ() << std::endl;
        std::cout << "Valide Potential = " << potentialValide << std::endl;

        ValidateFMMAlgoProc<OctreeClass,ContainerClass>(&tree,&treeValide);
    }
#endif


    // -----------------------------------------------------

    return 0;
}


