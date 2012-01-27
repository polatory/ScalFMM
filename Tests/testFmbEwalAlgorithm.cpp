// [--License--]

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Core/FFmmAlgorithmPeriodic.hpp"

#include "../Src/Components/FBasicKernels.hpp"
#include "../Src/Fmb/FFmbKernels.hpp"
#include "../Src/Fmb/FFmbComponents.hpp"

#include "../Src/Files/FEwalLoader.hpp"

/** Ewal particle is used in the gadget program
  * here we try to make the same simulation
  */
class EwalParticle : public FmbParticle {
public:
    // Type of particle
    enum Type{
        OW,
        HW,
        Undefined,
    };

private:
    Type type; //< current type
    int index; //< current index in array

public:
    // Basic constructor
    EwalParticle() : type(Undefined), index(-1) {
    }

    Type getType() const{
        return type;
    }

    void setType(const Type inType) {
        type = inType;
    }

    int getIndex() const{
        return index;
    }

    void setIndex( const int inIndex ){
        index = inIndex;
    }
};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef EwalParticle            ParticleClass;
    typedef FmbCell                 CellClass;
    typedef FVector<ParticleClass>  ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FBasicKernels<ParticleClass, CellClass, ContainerClass >        KernelClass;

    typedef FFmmAlgorithmPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test fmb algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 3);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 2);
    FTic counter;
    const char* const defaultFilename = "../Data/testEwal417.dt";
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

    FEwalLoader<ParticleClass> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tWidth : " << loader.getBoxWidth() << " \t center x : " << loader.getCenterOfBox().getX()
              << " y : " << loader.getCenterOfBox().getY() << " z : " << loader.getCenterOfBox().getZ()
              << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    counter.tic();

    ParticleClass * const particles = new ParticleClass[loader.getNumberOfParticles()];

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(particles[idxPart]);
        tree.insert(particles[idxPart]);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    //KernelClass kernels(NbLevels,loader.getBoxWidth());
    KernelClass kernels;
    FmmClass algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------
    {
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){
                ///potential += iter.data().getPotential() * iter.data().getPhysicalValue();
                ///forces += iter.data().getForces();

                std::cout << "index " << iter.data().getIndex() << " type " << iter.data().getType() << std::endl;
                std::cout << "x " << iter.data().getPosition().getX() << " y " << iter.data().getPosition().getY() << " z " << iter.data().getPosition().getZ() << std::endl;
                std::cout << "fx " << iter.data().getForces().getX() << " fy " << iter.data().getForces().getY() << " fz " << iter.data().getForces().getZ() << std::endl;
                std::cout << "potential " << iter.data().getPhysicalValue() << std::endl;

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());
    }

    // -----------------------------------------------------

    delete[] particles;

    return 0;
}


// [--END--]
