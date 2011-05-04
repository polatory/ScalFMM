// /!\ Please, you must read the license at the bottom of this page


#include "../Src/Utils/FMpi.hpp"
#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FList.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/F3DPosition.hpp"
#include "../Src/Utils/FAbstractSendable.hpp"

#include "../Src/Components/FFmaParticle.hpp"
#include "../Src/Components/FTestParticle.hpp"
#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestKernels.hpp"
#include "../Src/Extenssions/FExtendPhysicalValue.hpp"

#include "../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

#include "../Src/Files/FFmaLoader.hpp"

#include "../Src/Components/FBasicKernels.hpp"

#include <iostream>

#include <stdio.h>
#include <stdlib.h>


// Compile by : g++ testFmmAlgorithmProc.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmmAlgorithmProc.exe

/** This program show an example of use of
  * the fmm threaded + mpi algo
  * it also check that each particles is impacted each other particles
  */


/** Fmb class has to extend {FExtendForces,FExtendPotential,FExtendPhysicalValue}
  * Because we use fma loader it needs {FExtendPhysicalValue}
  */
class TestParticle : public FTestParticle, public FExtendPhysicalValue {
public:
};

class FTestCellPar : public FTestCell, public FAbstractSendable{
public :
    int bytesToSendUp() const{
        return sizeof(long);
    }
    int writeUp(void* const buffer, const int) const {
        const long tmpUp = getDataUp();
        memcpy(buffer,&tmpUp,bytesToSendUp());
        return bytesToSendUp();
    }
    int bytesToReceiveUp() const{
        return sizeof(long);
    }
    int readUp(void* const buffer, const int) {
        long tmpUp;
        memcpy(&tmpUp,buffer,bytesToSendUp());
        setDataUp(tmpUp);
        return bytesToReceiveUp();
    }

    int bytesToSendDown() const{
        return sizeof(long);
    }
    int writeDown(void* const buffer, const int) const {
        const long tmpDown = getDataDown();
        memcpy(buffer,&tmpDown,bytesToSendDown());
        return bytesToSendDown();
    }
    int bytesToReceiveDown() const{
        return sizeof(long);
    }
    int readDown(void* const buffer, const int) {
        long tmpDown;
        memcpy(&tmpDown,buffer,bytesToSendDown());
        setDataDown(tmpDown + getDataDown());
        return bytesToReceiveDown();
    }
};


/////////////////////////////////////////////////////////////////////////////
// Test function
/////////////////////////////////////////////////////////////////////////////

/** This function test the octree to be sure that the fmm algorithm
  * has worked completly.
  */
template<template< class ParticleClass, class CellClass, int OctreeHeight> class KernelClass,
        class ParticleClass, class CellClass,
        template<class ParticleClass> class LeafClass,
        int OctreeHeight, int SubtreeHeight>
void ValidateFMMAlgoProc(FOctree<ParticleClass, CellClass, LeafClass, OctreeHeight, SubtreeHeight>* const badTree,
                         FOctree<ParticleClass, CellClass, LeafClass, OctreeHeight, SubtreeHeight>* const valideTree,
                         FFmmAlgorithmThreadProc<FTestKernels, ParticleClass, CellClass, LeafClass, OctreeHeight, SubtreeHeight>*const fmm){
    std::cout << "Check Result\n";
    {
        typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level > 0 ; --level){
            int NbLeafs = 0;
            do{
                ++NbLeafs;
            } while(octreeIterator.moveRight());
            octreeIterator.gotoLeft();

            const int startIdx = fmm->getLeft(NbLeafs);
            const int endIdx = fmm->getRight(NbLeafs);
            // Check that each particle has been summed with all other

            for(int idx = 0 ; idx < startIdx ; ++idx){
                octreeIterator.moveRight();
                octreeIteratorValide.moveRight();
            }

            for(int idx = startIdx ; idx < endIdx ; ++idx){
                if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << "Error index are not equal!" << std::endl;
                }
                else{
                    if(octreeIterator.getCurrentCell()->getDataUp() != octreeIteratorValide.getCurrentCell()->getDataUp()){
                        std::cout << "M2M error at level " << level << " up bad " << octreeIterator.getCurrentCell()->getDataUp()
                                << " good " << octreeIteratorValide.getCurrentCell()->getDataUp() << " idx " << idx << std::endl;
                    }
                    if(octreeIterator.getCurrentCell()->getDataDown() != octreeIteratorValide.getCurrentCell()->getDataDown()){
                        std::cout << "L2L error at level " << level << " down bad " << octreeIterator.getCurrentCell()->getDataDown()
                                << " good " << octreeIteratorValide.getCurrentCell()->getDataDown() << " idx " << idx << std::endl;
                    }
                }

                octreeIterator.moveRight();
                octreeIteratorValide.moveRight();
            }

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
            typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();
            do{
                NbPart += octreeIterator.getCurrentListSrc()->getSize();
                ++NbLeafs;
            } while(octreeIterator.moveRight());
            std::cout << "There is " << NbPart << " particles on " << NbLeafs << " Leafs" << std::endl;
        }
        {
            const int startIdx = fmm->getLeft(NbLeafs);
            const int endIdx = fmm->getRight(NbLeafs);
            // Check that each particle has been summed with all other
            typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            for(int idx = 0 ; idx < startIdx ; ++idx){
                octreeIterator.moveRight();
            }

            for(int idx = startIdx ; idx < endIdx ; ++idx){
                typename FList<ParticleClass*>::BasicIterator iter(*octreeIterator.getCurrentListTargets());

                const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());

                while( iter.isValide() ){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem
                    if( (!isUsingTsm && iter.value()->getDataDown() != NbPart - 1) ||
                        (isUsingTsm && iter.value()->getDataDown() != NbPart) ){
                        std::cout << "Problem L2P + P2P, value on particle is : " << iter.value()->getDataDown() << "\n";
                    }
                    iter.progress();
                }
                octreeIterator.moveRight();
            }
        }
    }

    std::cout << "Done\n";
}

/** To print an octree
  * used to debug and understand how the values were passed
  */
template<template< class ParticleClass, class CellClass, int OctreeHeight> class KernelClass,
        class ParticleClass, class CellClass,
        template<class ParticleClass> class LeafClass,
        int OctreeHeight, int SubtreeHeight>
void print(FOctree<ParticleClass, CellClass, LeafClass, OctreeHeight, SubtreeHeight>* const valideTree){
    typename FOctree<ParticleClass, CellClass,LeafClass, OctreeHeight, SubtreeHeight>::Iterator octreeIterator(valideTree);
    for(int idxLevel = OctreeHeight - 1 ; idxLevel > 1 ; --idxLevel ){
        do{
            std::cout << "[" << octreeIterator.getCurrentGlobalIndex() << "] up:" << octreeIterator.getCurrentCell()->getDataUp() << " down:" << octreeIterator.getCurrentCell()->getDataDown() << "\t";
        } while(octreeIterator.moveRight());
        std::cout << "\n";
        octreeIterator.gotoLeft();
        octreeIterator.moveDown();
    }
}


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app( argc, argv);

    const int NbLevels = 10;//10;
    const int SizeSubLevels = 3;//3
    const char* const defaultFilename = "testLoaderFMA.fma"; //../../Data/ "testLoaderFMA.fma" "testFMAlgorithm.fma" Sphere.fma
    const char* filename;
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

    FFmaLoader<TestParticle> loader(filename);
    if(!loader.isValide()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    counter.tic();

    TestParticle*const particles = new TestParticle[loader.getNumberOfParticles()];
    TestParticle*const particlesValide = new TestParticle[loader.getNumberOfParticles()];

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&particles[idxPart]);
        particlesValide[idxPart] = particles[idxPart];
    }

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FOctree<TestParticle, FTestCellPar, FSimpleLeaf, NbLevels, SizeSubLevels> tree(loader.getBoxWidth(),loader.getCenterOfBox());

    FOctree<TestParticle, FTestCellPar, FSimpleLeaf, NbLevels, SizeSubLevels> treeValide(loader.getBoxWidth(),loader.getCenterOfBox());

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Inserting particles ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        tree.insert(&particles[idxPart]);
        treeValide.insert(&particlesValide[idxPart]);
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    FTestKernels<TestParticle, FTestCellPar, NbLevels> kernels;

    FFmmAlgorithmThreadProc<FTestKernels, TestParticle, FTestCellPar, FSimpleLeaf, NbLevels, SizeSubLevels> algo(app,&tree,&kernels);
    algo.execute();

    FFmmAlgorithmThread<FTestKernels, TestParticle, FTestCellPar, FSimpleLeaf, NbLevels, SizeSubLevels> algoValide(&treeValide,&kernels);
    algoValide.execute();

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    ValidateFMMAlgoProc<FTestKernels, TestParticle, FTestCellPar, FSimpleLeaf, NbLevels, SizeSubLevels>(&tree,&treeValide,&algo);

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    std::cout << "Deleting particles ..." << std::endl;
    counter.tic();
    delete [] particles;
    delete [] particlesValide;
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}


// [--LICENSE--]
