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
// @FUSE_MPI
// @FUSE_STARPU
// ================

#include "./Src/Utils/FMpi.hpp"

#include "./Src/Utils/FTic.hpp"
#include "./Src/Utils/FParameters.hpp"

#include "./Src/Containers/FOctree.hpp"
#include "./Src/Containers/FVector.hpp"

#include "./Src/Components/FTestParticle.hpp"
#include "./Src/Components/FTestCell.hpp"

#include "./Src/Core/FFmmAlgorithmStarpuProc.hpp"
#include "./Src/Core/FFmmAlgorithmThread.hpp"

#include "./Src/Components/FSimpleLeaf.hpp"

#include "./Src/Components/FTestKernels.hpp"

#include "./Src/Fmb/FFmbKernels.hpp"

#include "./Src/Files/FMpiFmaLoader.hpp"
#include "./Src/Files/FMpiTreeBuilder.hpp"
#include "./Src/Files/FFmaBinLoader.hpp"

#include "./Src/Components/FFmaParticle.hpp"
#include "./Src/Extensions/FExtendForces.hpp"
#include "./Src/Extensions/FExtendPotential.hpp"

#include "./Src/Components/FBasicCell.hpp"
#include "./Src/Fmb/FExtendFmbCell.hpp"

#include <starpu.h>

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


template<class OctreeClass, class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
KernelClass** FFmmAlgorithmStarpuProc<OctreeClass,ParticleClass,CellClass,ContainerClass,KernelClass,LeafClass>::globalKernels = 0;

// export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
// ./configure --enable-maxbuffers=190 --prefix=/home/bramas/starpu/starpu
// Compile With openmp : g++ testFmbAlgorithm.cpp ./Src/Utils/FDebug.cpp ./Src/Utils/FTrace.cpp -lgomp -fopenmp -lstarpu -lstarpumpi -O2 -o testFmbAlgorithm.exe
// compile with intel : mpicxx -I/home/bramas/starpu/starpu-0.9.2/include -I/home/bramas/starpu/starpu-0.9.2/mpi -openmp -lstarpu -lstarpumpi testFmbAlgorithm.cpp ./Src/Utils/FDebug.cpp ./Src/Utils/FTrace.cpp ./Src/Utils/FMath.cpp ./Src/Utils/F3DPosition.cpp -g -o testFmbAlgorithm.exe

////////////////////////////////////////////////////////////////
// Define classes
////////////////////////////////////////////////////////////////

//#define USE_TEST

#ifdef USE_TEST

template<class OctreeClass>
void ValidateTree(OctreeClass& realTree,
                        OctreeClass& treeValide, const FMpi::FComm& comm){
    FSize totalNbLeafs = 0;
    {

        typename OctreeClass::Iterator octreeIterator(&treeValide);
        octreeIterator.gotoBottomLeft();
        do {
            ++totalNbLeafs;
        } while(octreeIterator.moveRight());
    }

    const FSize myLeftLeaf = comm.getLeft(totalNbLeafs);
    const FSize myRightLeaf = comm.getRight(totalNbLeafs);

    //printf("%d should go from %d to %d leaf (on %d total leafs)\n", comm.processId(), myLeftLeaf, myRightLeaf, totalNbLeafs);

    typename OctreeClass::Iterator octreeIteratorValide(&treeValide);
    octreeIteratorValide.gotoBottomLeft();
    for(FSize idxLeaf = 0 ; idxLeaf < myLeftLeaf ; ++idxLeaf){
        if(!octreeIteratorValide.moveRight()){
            printf("Error cannot access to the left leaf %lld in the valide tree\n", myLeftLeaf);
        }
    }

    typename OctreeClass::Iterator octreeIterator(&realTree);
    octreeIterator.gotoBottomLeft();

    for(FSize idxLeaf = myLeftLeaf ; idxLeaf < myRightLeaf ; ++idxLeaf){
        if(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            printf("Error index are different, valide %lld invalid %lld\n",octreeIteratorValide.getCurrentGlobalIndex(),
                   octreeIterator.getCurrentGlobalIndex());
            break;
        }
        if(octreeIteratorValide.getCurrentListSrc()->getSize() != octreeIterator.getCurrentListSrc()->getSize()){
            printf("Error leafs do not have the same number of particles, valide %d, invalide %d\n",
                   octreeIteratorValide.getCurrentListSrc()->getSize(), octreeIterator.getCurrentListSrc()->getSize() );
        }

        //printf("index %lld with %d particles\n", octreeIteratorValide.getCurrentGlobalIndex(), octreeIteratorValide.getCurrentListSrc()->getSize());

        if(!octreeIteratorValide.moveRight() && idxLeaf != myRightLeaf - 1){
            printf("Error cannot valide tree end to early, idxLeaf %lld myRightLeaf %lld\n", idxLeaf, myRightLeaf);
            break;
        }

        if(!octreeIterator.moveRight() && idxLeaf != myRightLeaf - 1){
            printf("Error cannot test tree end to early, idxLeaf %lld myRightLeaf %lld\n", idxLeaf, myRightLeaf);
            break;
        }
    }

}

template<class OctreeClass, class ContainerClass>
void ValidateFMMAlgoProc(OctreeClass* const badTree,
                         OctreeClass* const valideTree){
    const int OctreeHeight = badTree->getHeight();
    {
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level > 0 ; --level){

            while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()) {
                octreeIteratorValide.moveRight();
            }

            FSize countCheck = 0;
            do{
                if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << "Error index are not equal!" << std::endl;
                }
                else{
                    if(octreeIterator.getCurrentCell()->getDataUp() != octreeIteratorValide.getCurrentCell()->getDataUp()){
                        std::cout << "M2M error at level " << level << " up bad " << octreeIterator.getCurrentCell()->getDataUp()
                                << " good " << octreeIteratorValide.getCurrentCell()->getDataUp() << " index " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }
                    if(octreeIterator.getCurrentCell()->getDataDown() != octreeIteratorValide.getCurrentCell()->getDataDown()){
                        std::cout << "L2L error at level " << level << " down bad " << octreeIterator.getCurrentCell()->getDataDown()
                                << " good " << octreeIteratorValide.getCurrentCell()->getDataDown() << " index " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }
                }
                ++countCheck;
            } while(octreeIteratorValide.moveRight() && octreeIterator.moveRight());

            // Check that each particle has been summed with all other

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            octreeIteratorValide.moveUp();
            octreeIteratorValide.gotoLeft();
        }
    }
    {
        {
            // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            do {
                if(octreeIterator.getCurrentListSrc()->getSize() != octreeIterator.getCurrentCell()->getDataUp()){
                    printf("P2M problem nb part %d data up %ld \n",
                           octreeIterator.getCurrentListSrc()->getSize(), octreeIterator.getCurrentCell()->getDataUp());
                }
            } while( octreeIterator.moveRight() );
        }
    }
    {
        FSize NbPart = 0;
        FSize NbLeafs = 0;
        { // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(valideTree);
            octreeIterator.gotoBottomLeft();
            do{
                NbPart += octreeIterator.getCurrentListSrc()->getSize();
                ++NbLeafs;
            } while(octreeIterator.moveRight());
        }
        {
            // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            do {
                typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
                const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());
                for(FSize idxPart = 0 ; idxPart < octreeIterator.getCurrentListTargets()->getSize() ; ++idxPart){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( (!isUsingTsm && iter.data().getDataDown() != NbPart - 1) ||
                        (isUsingTsm && iter.data().getDataDown() != NbPart) ){
                        std::cout << "Problem L2P + P2P, value on particle is : " << iter.data().getDataDown() <<
                                     " at pos " << idxPart << " index is " << octreeIterator.getCurrentGlobalIndex() << "\n";
                    }
                    iter.gotoNext();
                }
            } while( octreeIterator.moveRight());
        }
    }

    {
        // Check that each particle has been summed with all other
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator valideOctreeIterator(valideTree);
        valideOctreeIterator.gotoBottomLeft();
        while(valideOctreeIterator.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            valideOctreeIterator.moveRight();
        }

        do {
            if(valideOctreeIterator.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
                printf("Do not have the same index valide %lld invalide %lld \n",
                       valideOctreeIterator.getCurrentGlobalIndex(), octreeIterator.getCurrentGlobalIndex());
                break;
            }

            if(octreeIterator.getCurrentListTargets()->getSize() != valideOctreeIterator.getCurrentListTargets()->getSize()){
                printf("Do not have the same number of particle at leaf id %lld, valide %d invalide %d \n",
                       octreeIterator.getCurrentGlobalIndex(), valideOctreeIterator.getCurrentListTargets()->getSize(), octreeIterator.getCurrentListTargets()->getSize());
            }
            else {
                typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
                typename ContainerClass::BasicIterator iterValide(*valideOctreeIterator.getCurrentListTargets());

                for(int idxPart = 0 ; idxPart < octreeIterator.getCurrentListTargets()->getSize() ; ++idxPart){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( iter.data().getDataDown() != iterValide.data().getDataDown()){
                        std::cout << "Problem on leaf " << octreeIterator.getCurrentGlobalIndex() <<
                                     " part " << idxPart << " valide data down " << iterValide.data().getDataDown() <<
                                     " invalide " << iter.data().getDataDown() << "\n";
                        std::cout << "Data down for leaf is: valide " << valideOctreeIterator.getCurrentCell()->getDataDown()
                                  << " invalide " << octreeIterator.getCurrentCell()->getDataDown()
                                  << " size is: valide " <<  valideOctreeIterator.getCurrentListTargets()->getSize()
                                  << " invalide " << octreeIterator.getCurrentListTargets()->getSize() << std::endl;
                    }
                    iter.gotoNext();
                    iterValide.gotoNext();
                }
            }

        }while( octreeIterator.moveRight() && valideOctreeIterator.moveRight());
    }

}

#endif

/** Fmb class has to extend {FExtendForces,FExtendPotential,FExtendPhysicalValue}
  * Because we use fma loader it needs {FFmaParticle}
  */
class FmbParticle : public FExtendForces, public FFmaParticle, public FExtendPotential {
public:
};

/** Custom cell
  *
  */
class FmbCell : public FBasicCell, public FExtendFmbCell {
public:
};

class TestParticlePhys : public FTestParticle , public FExtendPhysicalValue {};

////////////////////////////////////////////////////////////////
// Typedefs
////////////////////////////////////////////////////////////////
#ifdef USE_TEST
    typedef TestParticlePhys             ParticleClass;
    typedef FTestCell                 CellClass;
#else
    typedef FmbParticle               ParticleClass;
    typedef FmbCell                   CellClass;
#endif

typedef FVector<ParticleClass >        ContainerClass;
typedef StarContainer<ContainerClass > StarContainerClass;
typedef StarVector<ParticleClass>      StarVectorClass;
typedef StarCell<CellClass>            StarCellClass;

typedef FSimpleLeaf<ParticleClass, StarContainerClass >                         LeafClass;
typedef FOctree<ParticleClass, StarCellClass, StarContainerClass , LeafClass >  OctreeClass;

#ifdef USE_TEST
typedef FTestKernels<ParticleClass, StarCellClass, StarVectorClass >                  KernelClass;
typedef FTestKernels<ParticleClass, StarCellClass, StarContainerClass >                  ValideKernelClass;
#else
typedef FFmbKernels<ParticleClass, StarCellClass, StarVectorClass >                   KernelClass;
typedef FFmbKernels<ParticleClass, StarCellClass, StarContainerClass >                   ValideKernelClass;
#endif


////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test fmb algorithm.\n";
    //////////////////////////////////////////////////////////////
    FMpi app(argc, argv);

    std::cout << "I am " << app.global().processId() << std::endl;

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 4);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 2);
    FTic counterTime;

    // -----------------------------------------------------

    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counterTime.tic();


    const char defaultFilename[] = "testLoaderFMA.fma"; //../../Data/ "testLoaderFMA.fma" "testFMAlgorithm.fma" Sphere.fma
    const char* filename;

    if(argc == 1){
        std::cout << "You have to give a .fma file in argument.\n";
        std::cout << "The program will try a default file : " << defaultFilename << "\n";
        filename = defaultFilename;
    }
    else{
        filename = argv[1];
        std::cout << "Opening : " << filename << "\n";
    }

    FMpiFmaLoader<ParticleClass> loader(filename, app.global());
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }


    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    counterTime.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counterTime.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    if( app.global().processCount() != 1){
        //////////////////////////////////////////////////////////////////////////////////
        // Build tree from mpi loader
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Build Tree ..." << std::endl;

        FMpiTreeBuilder<ParticleClass>::LoaderToTree(app.global(), loader, tree);

        //////////////////////////////////////////////////////////////////////////////////
    }
    else{
        ParticleClass partToInsert;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(partToInsert);
            tree.insert(partToInsert);
        }
    }

    // -----------------------------------------------------



    std::cout << "Working on particles ..." << std::endl;
    counterTime.tic();

#ifdef USE_TEST
    KernelClass kernel;
#else
    KernelClass kernel(tree.getHeight(), 1.0);
#endif
    FFmmAlgorithmStarpuProc<OctreeClass,ParticleClass,StarCellClass,StarContainerClass,KernelClass,LeafClass> algo(app.global(),  &tree, &kernel);
    algo.execute();

    counterTime.tac();
    std::cout << "Done  " << "(@Algorithm = " << counterTime.elapsed() << "s)." << std::endl;


    // Check result
#ifdef USE_TEST

    OctreeClass treeValide(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
    {
        FFmaBinLoader<ParticleClass> loaderSeq(filename);
        ParticleClass partToInsert;
        for(FSize idxPart = 0 ; idxPart < loaderSeq.getNumberOfParticles() ; ++idxPart){
            loaderSeq.fillParticle(partToInsert);
            treeValide.insert(partToInsert);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Check particles in tree
    //////////////////////////////////////////////////////////////////////////////////
    std::cout << "Validate tree ..." << std::endl;

    ValidateTree(tree, treeValide, app.global());


    ValideKernelClass validekernel;

    FFmmAlgorithmThread<OctreeClass,ParticleClass,StarCellClass,StarContainerClass,ValideKernelClass,LeafClass> algoValide(&treeValide,&validekernel);
    algoValide.execute();

    ValidateFMMAlgoProc<OctreeClass,ContainerClass>(&tree,&treeValide);

#else
    { // get sum forces&potential
        FReal potential = 0;
        F3DPosition forces;
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
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

