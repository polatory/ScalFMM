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
// @FUSE_MPI
// ================

#include "../../Src/Utils/FMpi.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"
#include "../../Src/Components/FAbstractSendable.hpp"

#include "../../Src/Components/FFmaParticle.hpp"
#include "../../Src/Components/FTestParticle.hpp"
#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"
#include "../../Src/Extensions/FExtendPhysicalValue.hpp"


#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Files/FFmaBinLoader.hpp"
#include "../../Src/Files/FMpiFmaLoader.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Components/FBasicKernels.hpp"

#include <iostream>
#include <cstdio>
#include <cstdlib>


/** This program show an example of use of the fmm threaded + mpi algo
  * it also check that each particles is impacted each other particles
  */

/////////////////////////////////////////////////////////////////////////////
// Test function
/////////////////////////////////////////////////////////////////////////////

// Check if tree is built correctly
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



/** This function tests the octree to be sure that the fmm algorithm
  * has worked completly.
  */
template<class OctreeClass, class ContainerClass, class FmmClassProc>
void ValidateFMMAlgoProc(OctreeClass* const badTree,
                         OctreeClass* const valideTree,
                         FmmClassProc* const fmm){
    const int OctreeHeight = badTree->getHeight();
    {
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level > 0 && fmm->hasWorkAtLevel(level) ; --level){

            while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()) {
                octreeIteratorValide.moveRight();
            }

            while(octreeIteratorValide.getCurrentGlobalIndex() != fmm->getWorkingInterval(level).min){
                octreeIteratorValide.moveRight();
                octreeIterator.moveRight();
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
        {
            // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            do {
                if(octreeIterator.getCurrentListSrc()->getSize() != octreeIterator.getCurrentCell()->getDataUp()){
                    printf("P2M problem nb part %d data up %lld \n",
                           octreeIterator.getCurrentListSrc()->getSize(), octreeIterator.getCurrentCell()->getDataUp());
                }
            } while( octreeIterator.moveRight() );
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


/** To print an octree
  * used to debug and understand how the values were passed
  */
template<class OctreeClass>
void print(OctreeClass* const valideTree){
    typename OctreeClass::Iterator octreeIterator(valideTree);
    for(int idxLevel = valideTree->getHeight() - 1 ; idxLevel > 1 ; --idxLevel ){
        do{
            std::cout << "[" << octreeIterator.getCurrentGlobalIndex() << "] up:" << octreeIterator.getCurrentCell()->getDataUp() << " down:" << octreeIterator.getCurrentCell()->getDataDown() << "\t";
        } while(octreeIterator.moveRight());
        std::cout << "\n";
        octreeIterator.gotoLeft();
        octreeIterator.moveDown();
    }
}

/////////////////////////////////////////////////////////////////////
// Types
/////////////////////////////////////////////////////////////////////



/** class has to extend {FExtendForces,FExtendPotential,FExtendPhysicalValue}
  * Because we use fma loader it needs {FExtendPhysicalValue}
  */
class TestParticle : public FTestParticle, public FExtendPhysicalValue {
public:
    /** Save current object */
    void save(FBufferWriter& buffer) const {
        FTestParticle::save(buffer);
        FExtendPhysicalValue::save(buffer);
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        FTestParticle::restore(buffer);
        FExtendPhysicalValue::restore(buffer);
    }
};

/////////////////////////////////////////////////////////////////////
// Define the classes to use
/////////////////////////////////////////////////////////////////////

typedef TestParticle               ParticleClass;
typedef FTestCell                  CellClass;
typedef FVector<ParticleClass>     ContainerClass;

typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
typedef FTestKernels<ParticleClass, CellClass, ContainerClass >         KernelClass;

typedef FFmmAlgorithmThread<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
typedef FFmmAlgorithmThreadProc<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClassProc;

/////////////////////////////////////////////////////////////////////
// Main
/////////////////////////////////////////////////////////////////////

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app( argc, argv);

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const char* const defaultFilename = (sizeof(FReal) == sizeof(float))?
                                    "../../Data/test20k.bin.fma.single":
                                    "../../Data/test20k.bin.fma.double";
    const char* const filename = FParameters::getStr(argc,argv,"-f", defaultFilename);
    std::cout << "Opening : " << filename << "\n";

    FMpiFmaLoader<ParticleClass> loader(filename,app.global());
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // The real tree to work on
    OctreeClass realTree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    if( app.global().processCount() != 1){
        //////////////////////////////////////////////////////////////////////////////////
        // Build tree from mpi loader
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Build Tree ..." << std::endl;
        counter.tic();

        FMpiTreeBuilder<ParticleClass>::LoaderToTree(app.global(), loader, realTree);

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
    }    
    else{
        loader.fillTree(realTree);
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Create real tree
    //////////////////////////////////////////////////////////////////////////////////

    OctreeClass treeValide(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
    {
        FFmaBinLoader<ParticleClass> loaderSeq(filename);
        loaderSeq.fillTree(treeValide);
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Check particles in tree
    //////////////////////////////////////////////////////////////////////////////////
    std::cout << "Validate tree ..." << std::endl;
    counter.tic();

    ValidateTree(realTree, treeValide, app.global());

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working parallel particles ..." << std::endl;
    counter.tic();

    KernelClass kernels;

    FmmClassProc algo(app.global(),&realTree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working sequential particles ..." << std::endl;
    counter.tic();

    FmmClass algoValide(&treeValide,&kernels);
    algoValide.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Checking data ..." << std::endl;
    counter.tic();

    ValidateFMMAlgoProc<OctreeClass,ContainerClass, FmmClassProc>(&realTree,&treeValide,&algo);

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}



