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

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"


#include "../Src/Components/FBasicKernels.hpp"

// Compile by : g++ testFmmDemonstration.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmmDemonstration.exe

// My cell is actually a basic cell => minimum of data
class MyCell : public FBasicCell {
};

// My fack particle is a basic particle and! a index
// that correspond to the index in the user particles array
class MyFackParticle : public FBasicParticle {
    int myIndex;
public:
    MyFackParticle() : myIndex(0) {
    }
    void setIndex(const int inIndex){
        this->myIndex = inIndex;
    }
    int getIndex() const{
        return this->myIndex;
    }
};

// My leaf store the indexes of the particles it receives
// in a vector
class MyLeaf : public FAbstractLeaf<MyFackParticle, FVector<int> > {
    FVector<int> indexes;
public:
    void push(const MyFackParticle& particle){
        indexes.push( particle.getIndex() );
    }
    FVector<int>* getSrc(){
        return &indexes;
    }
    FVector<int>* getTargets(){
        return &indexes;
    }
};

// My kernel actually does nothing except showing how to retreive data from an
// array from the indexes vector giving by the leaf in the P2M
template< class ParticleClass, class CellClass, class ContainerClass>
class MyKernel : public FAbstractKernels<ParticleClass,CellClass,ContainerClass> {
    FBasicParticle*const realParticles;
public:
    MyKernel(FBasicParticle*const inRealParticles): realParticles(inRealParticles) {
    }

    void P2M(CellClass* const , const ContainerClass* const particlesIndexes) {
        typename ContainerClass::ConstBasicIterator iter(*particlesIndexes);
        while( iter.hasNotFinished() ){
            realParticles[iter.data()].getPosition();
            iter.gotoNext();
        }
    }

    void M2M(CellClass* const FRestrict , const CellClass*const FRestrict *const FRestrict , const int ) {
    }

    void M2L(CellClass* const FRestrict , const CellClass* [], const int , const int ) {
    }

    void L2L(const CellClass* const FRestrict , CellClass* FRestrict *const FRestrict  , const int ) {
    }

    void L2P(const CellClass* const , ContainerClass* const ){
    }

    void P2P(ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
                     const ContainerClass* const [26], const int ) {
    }

    void P2P(const MortonIndex ,
                     ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
                     ContainerClass* const [26], const MortonIndex [26], const int ){
    }
};


int main(int argc, char ** argv){
    typedef MyFackParticle    ParticleClass;
    typedef MyCell            CellClass;
    typedef FVector<int>      ContainerClass;
    typedef MyLeaf            LeafClass;

    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef MyKernel<ParticleClass, CellClass, ContainerClass >         KernelClass;
    typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;

    //////////////////////////////////////////////////////////////
    ///////////////////////What we do/////////////////////////////

    std::cout << ">> This executable has to be used to demonstrate the use of scalfmm.\n";

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const long NbPart = FParameters::getValue(argc,argv,"-pn", 20L);

    FTic counter;

    srand ( 1 );

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    OctreeClass tree(NbLevels, SizeSubLevels,1.0,F3DPosition(0.5,0.5,0.5));
    FBasicParticle*const realsParticles = new FBasicParticle[NbPart];

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    {
        ParticleClass particleToFill;
        for(int idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particleToFill.setPosition(FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX);
            particleToFill.setIndex(idxPart);

            tree.insert(particleToFill);

            realsParticles[idxPart].setPosition(particleToFill.getPosition());
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Which particles are in wich leafs ..." << std::endl;
    counter.tic();

    typename OctreeClass::Iterator octreeIterator(&tree);
    octreeIterator.gotoBottomLeft();
    do{
        typename ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
        const MortonIndex indexAtThisLeaf = octreeIterator.getCurrentGlobalIndex();

        while( iter.hasNotFinished() ){
            std::cout << "Particles with index " << iter.data() << " has a morton index of " << indexAtThisLeaf << std::endl;

            const F3DPosition& particlePosition = realsParticles[iter.data()].getPosition();
            std::cout << "\t The real position of this particle is (" << particlePosition.getX() << ";" << particlePosition.getY() << ";" << particlePosition.getZ() << ")" << std::endl;

            iter.gotoNext();
        }
    } while(octreeIterator.moveRight());

    counter.tac();
    std::cout << "Done  " << "(@Counting = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels(realsParticles);
    FmmClass algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    delete [] realsParticles;

    return 0;
}


// [--LICENSE--]
