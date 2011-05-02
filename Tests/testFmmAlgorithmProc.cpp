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


// Compile by : g++ testFmmAlgorithmProc.cpp ../Src/Utils/FAssertable.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmmAlgorithmProc.exe

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

    TestParticle* particles = new TestParticle[loader.getNumberOfParticles()];
    TestParticle* particlesValide = new TestParticle[loader.getNumberOfParticles()];

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

    algo.ValidateFMMAlgoProc(&treeValide);

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
