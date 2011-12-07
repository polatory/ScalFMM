
// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Src/Utils/FParameters.hpp"
#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Components/FTestParticle.hpp"
#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestKernels.hpp"

#include "../Src/Core/FFmmAlgorithmThreadProcPeriodic.hpp"

#include "../Src/Components/FBasicKernels.hpp"
#include "../Src/Files/FMpiTreeBuilder.hpp"

#include "../Src/Utils/FAbstractSendable.hpp"

// Compile by : g++ testFmmAlgorithm.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmmAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is impacted each other particles
  */



template< class ParticleClass, class CellClass, class ContainerClass>
class FTestPeriodicKernels : public FTestKernels<ParticleClass,CellClass,ContainerClass> {
public:

    /** Before Downward */
    void M2L(CellClass* const FRestrict pole, const CellClass* distantNeighbors[189], FTreeCoordinate [189], const int size, const int ) {
        // The pole is impacted by what represent other poles
        for(int idx = 0 ; idx < size ; ++idx){
            pole->setDataDown(pole->getDataDown() + distantNeighbors[idx]->getDataUp());
        }
    }


    /** After Downward */
    void P2P(const MortonIndex ,
             ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
             ContainerClass* const directNeighborsParticles[26], const FTreeCoordinate [26], const int size) {

        // Each particles targeted is impacted by the particles sources
        long long int inc = sources->getSize();
        if(targets == sources){
            inc -= 1;
        }
        for(int idx = 0 ; idx < size ; ++idx){
            inc += directNeighborsParticles[idx]->getSize();
        }

        typename ContainerClass::BasicIterator iter(*targets);
        while( iter.hasNotFinished() ){
            iter.data().setDataDown(iter.data().getDataDown() + inc);
            iter.gotoNext();
        }

    }
};

class TestCell : public FTestCell , public FAbstractSendable {
public:
    static const int SerializedSizeUp = sizeof(long long int);
    void serializeUp(void* const buffer) const {
        *(long long int*)buffer = this->dataUp;
    }
    void deserializeUp(const void* const buffer){
        this->dataUp = *(long long int*)buffer;
    }

    static const int SerializedSizeDown = sizeof(long);
    void serializeDown(void* const buffer) const {
        *(long long int*)buffer = this->dataDown;
    }
    void deserializeDown(const void* const buffer){
        this->dataDown = *(long long int*)buffer;
    }
};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FTestParticle               ParticleClass;
    typedef TestCell                    CellClass;
    typedef FVector<ParticleClass>      ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestPeriodicKernels<ParticleClass, CellClass, ContainerClass >         KernelClass;

    typedef FFmmAlgorithmThreadProcPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels          = FParameters::getValue(argc,argv,"-h", 7);
    const int SizeSubLevels     = FParameters::getValue(argc,argv,"-sh", 3);
    const long NbPartPerBoxesPerProc   = FParameters::getValue(argc,argv,"-nb", 1);

    FMpi app(argc, argv);

    const long NbPartPerBoxes = NbPartPerBoxesPerProc * app.global().processCount();

    FTic counter;

    srand ( 1 ); // volontary set seed to constant

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    const FReal BoxWidth = 1.0;
    const FReal BoxCenter = 0.5;

    OctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, F3DPosition(BoxCenter,BoxCenter,BoxCenter));

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPartPerBoxes << " particles per boxes ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    long NbPart = 0;
    {
        const long NbSmallBoxesPerSide = (1 << (NbLevels-1));
        const FReal SmallBoxWidth = BoxWidth / FReal(NbSmallBoxesPerSide);
        const FReal SmallBoxWidthDiv2 = SmallBoxWidth / 2;

        NbPart = NbSmallBoxesPerSide * NbSmallBoxesPerSide * NbSmallBoxesPerSide * NbPartPerBoxesPerProc;

        FTestParticle*const particles = new FTestParticle[NbPart];
        int indexPart = 0;

        for(int idxX = 0 ; idxX < NbSmallBoxesPerSide ; ++idxX){
            for(int idxY = 0 ; idxY < NbSmallBoxesPerSide ; ++idxY){
                for(int idxZ = 0 ; idxZ < NbSmallBoxesPerSide ; ++idxZ){
                    for(int idxPart = 0 ; idxPart < NbPartPerBoxesPerProc ; ++idxPart ){
                        particles[indexPart++].setPosition(FReal(idxX)*SmallBoxWidth + SmallBoxWidthDiv2,
                                                   FReal(idxY)*SmallBoxWidth + SmallBoxWidthDiv2,
                                                   FReal(idxZ)*SmallBoxWidth + SmallBoxWidthDiv2);
                    }
                }
            }
        }

        FMpiTreeBuilder<ParticleClass>::ArrayToTree(app.global(), particles, NbPart, F3DPosition(BoxCenter,BoxCenter,BoxCenter), BoxWidth, tree);
        delete[] particles;
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    const int PeriodicDeep = 3;
    // FTestKernels FBasicKernels
    KernelClass kernels;
    //FFmmAlgorithm FFmmAlgorithmThread
    FmmClass algo( app.global(), &tree, &kernels, PeriodicDeep);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    { // Check that each particle has been summed with all other
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            if( octreeIterator.getCurrentListTargets()->getSize() != NbPartPerBoxes){
                std::cout << "Incorrect number of particles per leaf, should be " << NbPartPerBoxes <<
                             " iter.data().getDataDown() "<< octreeIterator.getCurrentListTargets()->getSize() << std::endl;
            }
        } while(octreeIterator.moveRight());
    }
    { // Check that each particle has been summed with all other
        long long int counterNbPart = 0;

        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            // on each leaf we should have the same number of particles
            if(octreeIterator.getCurrentCell()->getDataUp() != octreeIterator.getCurrentListSrc()->getSize()
                    || NbPartPerBoxes != octreeIterator.getCurrentCell()->getDataUp() ){
                    std::cout << "Problem P2M NbPartPerBoxes = " << NbPartPerBoxes <<
                                 " Data up = " << octreeIterator.getCurrentCell()->getDataUp() <<
                                 " Size = " << octreeIterator.getCurrentListSrc()->getSize() << "\n";
            }
            // we also count the number of particles.
            counterNbPart += octreeIterator.getCurrentListSrc()->getSize();
        } while(octreeIterator.moveRight());

        const long long int allPart = app.global().reduceSum( counterNbPart );
        if( app.global().processId() == 0 && allPart != NbPart * app.global().processCount()){
            std::cout << "Problem global nb part, counter = " << allPart << " created = " << NbPart * app.global().processCount() << std::endl;
        }
    }
    { // Ceck if there is number of NbPart summed at level 1
        long long particlesPerBox = NbPartPerBoxes;

        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        for(int idxLevel = NbLevels - 1 ; idxLevel >= 1 ; --idxLevel ){
            if(algo.hasWorkAtLevel(idxLevel)){
                while(octreeIterator.getCurrentGlobalIndex() != algo.getWorkingInterval(idxLevel).min){
                    octreeIterator.moveRight();
                }

                do{
                    if(octreeIterator.getCurrentCell()->getDataUp() != particlesPerBox){
                        std::cout << "Problem M2M particlesPerBox = " << particlesPerBox <<
                                     " Data up = " << octreeIterator.getCurrentCell()->getDataUp() <<
                                     " level = " << idxLevel << "\n";
                    }
                } while(octreeIterator.moveRight());
            }

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            particlesPerBox *= 8;
        }
    }
    {
        long long int counterL2L = 0;
        if( PeriodicDeep ){
            long long int particlesPerBox = NbPartPerBoxes * FMath::pow(8,NbLevels-1);

            long long int counterUp[PeriodicDeep];
            for( int idxLevel = 0 ; idxLevel < PeriodicDeep ; ++idxLevel ){
                counterUp[idxLevel] = particlesPerBox;
                particlesPerBox *= 8;
            }

            long long int counterM2L[PeriodicDeep];
            for( int idxLevel = 0 ; idxLevel < PeriodicDeep ; ++idxLevel ){
                counterM2L[idxLevel] = counterUp[idxLevel] * 189;
            }

            for( int idxLevel = PeriodicDeep - 1 ; idxLevel >= 0 ; --idxLevel ){
                counterL2L += counterM2L[idxLevel];
            }
        }
        {
            long long int particlesPerBox = NbPartPerBoxes * FMath::pow(8,NbLevels-2);

            typename OctreeClass::Iterator octreeIterator(&tree);
            for(int idxLevel = 1 ; idxLevel < NbLevels ; ++idxLevel ){
                counterL2L = particlesPerBox * 189 + counterL2L;

                if(algo.hasWorkAtLevel(idxLevel)){
                    while(octreeIterator.getCurrentGlobalIndex() != algo.getWorkingInterval(idxLevel).min){
                        octreeIterator.moveRight();
                    }

                    do{
                        if(octreeIterator.getCurrentCell()->getDataDown() != counterL2L){
                            std::cout << app.global().processId() <<  " >>Problem L2L counterL2L = " << counterL2L <<
                                         " Data Down = " << octreeIterator.getCurrentCell()->getDataDown() <<
                                         " level = " << idxLevel << "\n";
                        }
                    } while(octreeIterator.moveRight());
                }

                octreeIterator.gotoLeft();
                octreeIterator.moveDown();

                particlesPerBox /= 8;
            }
        }
    }
    { // Check that each particle has been summed with all other
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            const long long int sumParticles = octreeIterator.getCurrentCell()->getDataDown() + (26 * NbPartPerBoxes) + (NbPartPerBoxes - 1);

            while( iter.hasNotFinished() ){
                if( sumParticles != iter.data().getDataDown()){
                    std::cout << "P2P probleme, should be " << sumParticles <<
                                 " iter.data().getDataDown() "<< iter.data().getDataDown() << std::endl;
                }

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());
    }

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}


// [--LICENSE--]
