// /!\ Please, you must read the license at the bottom of this page

#include "../Src/Utils/FMpi.hpp"
#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"
#include "../Src/Utils/FParameters.hpp"
#include "../Src/Utils/FGlobal.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/F3DPosition.hpp"
#include "../Src/Utils/FAbstractSendable.hpp"

#include "../Src/Components/FFmaParticle.hpp"
#include "../Src/Components/FTestParticle.hpp"
#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestKernels.hpp"
#include "../Src/Extensions/FExtendPhysicalValue.hpp"

#include "../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FMpiFmaLoader.hpp"
#include "../Src/Files/FMpiTreeBuilder.hpp"

#include "../Src/Components/FBasicKernels.hpp"

#include <iostream>

#include <stdio.h>
#include <stdlib.h>


// Compile by : g++ testFmmAlgorithmProc.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmmAlgorithmProc.exe

/** This program show an example of use of
  * the fmm threaded + mpi algo
  * it also check that each particles is impacted each other particles
  */

/////////////////////////////////////////////////////////////////////////////
// Test function
/////////////////////////////////////////////////////////////////////////////

// Check if tree is built correctly
template<class OctreeClass>
void ValidateTree(OctreeClass& realTree,
                        OctreeClass& treeValide, FMpi& app){
    int totalNbLeafs = 0;
    {

        typename OctreeClass::Iterator octreeIterator(&treeValide);
        octreeIterator.gotoBottomLeft();
        do {
            ++totalNbLeafs;
        } while(octreeIterator.moveRight());
    }

    const int myLeftLeaf = app.getLeft(totalNbLeafs);
    const int myRightLeaf = app.getRight(totalNbLeafs);

    //printf("%d should go from %d to %d leaf (on %d total leafs)\n", app.processId(), myLeftLeaf, myRightLeaf, totalNbLeafs);

    typename OctreeClass::Iterator octreeIteratorValide(&treeValide);
    octreeIteratorValide.gotoBottomLeft();
    for(int idxLeaf = 0 ; idxLeaf < myLeftLeaf ; ++idxLeaf){
        if(!octreeIteratorValide.moveRight()){
            printf("Error cannot access to the left leaf %d in the valide tree\n", myLeftLeaf);
        }
    }

    typename OctreeClass::Iterator octreeIterator(&realTree);
    octreeIterator.gotoBottomLeft();

    for(int idxLeaf = myLeftLeaf ; idxLeaf < myRightLeaf ; ++idxLeaf){
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
            printf("Error cannot valide tree end to early, idxLeaf %d myRightLeaf %d\n", idxLeaf, myRightLeaf);
            break;
        }

        if(!octreeIterator.moveRight() && idxLeaf != myRightLeaf - 1){
            printf("Error cannot test tree end to early, idxLeaf %d myRightLeaf %d\n", idxLeaf, myRightLeaf);
            break;
        }
    }

}



/** This function test the octree to be sure that the fmm algorithm
  * has worked completly.
  */
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

            int countCheck = 0;
            do{
                if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << "Error index are not equal!" << std::endl;
                }
                else{
                    /*if(octreeIterator.getCurrentCell()->getDataUp() != octreeIteratorValide.getCurrentCell()->getDataUp()){
                        std::cout << "M2M error at level " << level << " up bad " << octreeIterator.getCurrentCell()->getDataUp()
                                << " good " << octreeIteratorValide.getCurrentCell()->getDataUp() << " index " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                    }*/
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
        int NbPart = 0;
        int NbLeafs = 0;
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
                for(int idxPart = 0 ; idxPart < octreeIterator.getCurrentListTargets()->getSize() ; ++idxPart){
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
                    printf("P2M problem nb part %d data up %ld \n",
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

/** Fmb class has to extend {FExtendForces,FExtendPotential,FExtendPhysicalValue}
  * Because we use fma loader it needs {FExtendPhysicalValue}
  */
class TestParticle : public FTestParticle, public FExtendPhysicalValue {
};


class TestCell : public FTestCell , public FAbstractSendable {
public:
    static const int SerializedSizeUp = sizeof(long);
    void serializeUp(void* const buffer) const {
        *(long*)buffer = this->dataUp;
    }
    void deserializeUp(const void* const buffer){
        this->dataUp = *(long*)buffer;
    }

    static const int SerializedSizeDown = sizeof(long);
    void serializeDown(void* const buffer) const {
        *(long*)buffer = this->dataDown;
    }
    void deserializeDown(const void* const buffer){
        this->dataDown = *(long*)buffer;
    }
};

typedef TestParticle               ParticleClass;
typedef TestCell                   CellClass;
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

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 9);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    char defaultFilename[] = "testLoaderFMA.fma"; //../../Data/ "testLoaderFMA.fma" "testFMAlgorithm.fma" Sphere.fma
    char* filename;
    FTic counter;

    if(argc == 1){
        std::cout << "You have to give a .fma file in argument.\n";
        std::cout << "The program will try a default file : " << defaultFilename << "\n";
        filename = defaultFilename;
    }
    else{
        filename = argv[1];
        std::cout << "Opening : " << filename << "\n";
    }

    FMpiFmaLoader<ParticleClass> loader(filename,app);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // The real tree to work on
    OctreeClass realTree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    if( app.processCount() != 1){
        //////////////////////////////////////////////////////////////////////////////////
        // We sort our particles
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Create intervals ..." << std::endl;
        counter.tic();

        FMpiTreeBuilder<ParticleClass> builder;

        builder.splitAndSortFile(loader, NbLevels);

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
        // Now we can build the real tree
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Create real tree ..." << std::endl;
        counter.tic();

        builder.intervalsToTree(realTree);

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
    }    
    else{
        FFmaBinLoader<ParticleClass> loaderSeq(filename);
        ParticleClass partToInsert;
        for(int idxPart = 0 ; idxPart < loaderSeq.getNumberOfParticles() ; ++idxPart){
            loaderSeq.fillParticle(partToInsert);
            realTree.insert(partToInsert);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Create real tree
    //////////////////////////////////////////////////////////////////////////////////

    OctreeClass treeValide(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
    {
        FFmaBinLoader<ParticleClass> loaderSeq(filename);
        ParticleClass partToInsert;
        for(int idxPart = 0 ; idxPart < loaderSeq.getNumberOfParticles() ; ++idxPart){
            loaderSeq.fillParticle(partToInsert);
            treeValide.insert(partToInsert);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Check particles in tree
    //////////////////////////////////////////////////////////////////////////////////
    std::cout << "Validate tree ..." << std::endl;
    counter.tic();

    ValidateTree(realTree, treeValide, app);

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working parallel particles ..." << std::endl;
    counter.tic();

    KernelClass kernels;

    FmmClassProc algo(app,&realTree,&kernels);
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

    ValidateFMMAlgoProc<OctreeClass,ContainerClass>(&realTree,&treeValide);

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}


// [--LICENSE--]
