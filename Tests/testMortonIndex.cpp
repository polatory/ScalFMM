// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Components/FBasicParticle.hpp"

#include "../Src/Components/FBasicCell.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"
#include "../Src/Core/FFmmAlgorithmThreadUs.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Fmb/FFmbKernels.hpp"


#include "../Src/Files/FHLoader.hpp"


/** Basic function to convert a morton index in decimal string */
std::string MortonToBinary(MortonIndex index, int level){
    std::string str;
    int bits = 1 << ((level * 3) - 1);
    int dim = 0;
    while(bits){
        if(index & bits) str.append("1");
        else str.append("0");
        bits >>= 1;
        // we put a dot each 3 values
        if(++dim == 3){
            str.append(".");
            dim = 0;
        }
    }
    return str;
}

/** This program show an example of use of
  * the moront indexing from a file
  */

class Particle : public FBasicParticle {
    char data;

public:
    char getData() const{
        return data;
    }
    void setData(const char inData){
        this->data = inData;
    }
};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to know the box used.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = 6;
    const int SizeSubLevels = 2;
    FTic counter;
    const char* const defaultFilename = "../Data/testH.txt";
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

    FHLoader<Particle> loader(filename);
    if(!loader.hasNotFinished()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    FOctree<Particle, FBasicCell, FVector, FSimpleLeaf> tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating and inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "Width of box is " << loader.getBoxWidth() << std::endl;
    std::cout << "Center of box is x " << loader.getCenterOfBox().getX()  << " y " << loader.getCenterOfBox().getY() << " z " << loader.getCenterOfBox().getZ()<< std::endl;
    counter.tic();

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        Particle particle;
        loader.fillParticle(particle);

        tree.insert(particle);
    }

    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;


    // -----------------------------------------------------


    { // print indexes
        FOctree<Particle, FBasicCell, FVector, FSimpleLeaf>::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            const MortonIndex currentIndex = octreeIterator.getCurrentGlobalIndex();
            std::cout << "Current Morton Index : " << currentIndex << " or in binary " << MortonToBinary(currentIndex,NbLevels-1) << std::endl;
            std::cout << "Particles :" << std::endl;

            FVector<Particle>::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){

                printf("\tx = %e y = %e z = %e data = %c\n",iter.data().getPosition().getX(),iter.data().getPosition().getY(),iter.data().getPosition().getZ(),iter.data().getData());

                iter.gotoNext();
            }

        } while(octreeIterator.moveRight());
    }


    // -----------------------------------------------------

    return 0;
}


// [--LICENSE--]
