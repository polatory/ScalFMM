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

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Components/FBasicParticle.hpp"

#include "../Src/Components/FBasicCell.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Utils/FGlobal.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Files/FHLoader.hpp"

static const int PrintPrecision = 6;

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
struct AtomIndex  {
    int          globalIndex ;
    MortonIndex  indexMorton;
};

bool compareIndex(const AtomIndex& a, const AtomIndex& b){
    return (a.indexMorton < b.indexMorton);
}

class Particle : public FBasicParticle {
    char         data   ;
    AtomIndex    index ;
    //
public:
    char getData() const{
        return data;
    }
    void setData(const char inData){
        this->data = inData;
    }

    MortonIndex mortonIndex() const{
        return index.indexMorton;
    }
    void setMortonIndex(const MortonIndex inIndex) {
        index.indexMorton = inIndex;
    }

    int globalIndex() const{
        return index.globalIndex;
    }
    void setGlobalIndex(const int inIndex) {
        index.globalIndex = inIndex;
    }

};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FVector<Particle>      ContainerClass;
    typedef FSimpleLeaf<Particle, ContainerClass >                     LeafClass;
    typedef FOctree<Particle, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to know the box used.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/testMortonIndex.txt");
    std::cout << "Opening : " << filename << "\n";

    FHLoader<Particle> loader(filename);
    if(!loader.isOpen() ){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating and inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "Width of box is " << loader.getBoxWidth() << std::endl;
    std::cout << "Center of box is x " << loader.getCenterOfBox().getX()  << " y " << loader.getCenterOfBox().getY() << " z " << loader.getCenterOfBox().getZ()<< std::endl;
    counter.tic();
    {
        Particle particle;
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(particle);
            particle.setGlobalIndex( idxPart+1 );
            particle.setMortonIndex( -1 );
            tree.insert(particle);
        }
    }
    counter.tac();
    std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    counter.tic();
    int idx = 1 ;
    std::vector<AtomIndex> permutation( int(loader.getNumberOfParticles()) ) ;

    { // print indexes

        // -----------------------------------------------------

        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        idx = 0;
        do{
            const MortonIndex currentIndex = octreeIterator.getCurrentGlobalIndex();
            // std::cout << "Current Morton Index : " << currentIndex << " or in binary " << MortonToBinary(currentIndex,NbLevels-1) << std::endl;
            // std::cout << "Particles :" << std::endl;

            typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){
                iter.data().setMortonIndex(currentIndex);
                //	  iter.data().mortonIndex()       = currentIndex ;
                permutation[idx].indexMorton   = currentIndex ;
                permutation[idx].globalIndex   = idx + 1;
                ++idx;
                //	  printf("\tx = %e y = %e z = %e data = %c\n",iter.value()->getPosition().getX(),iter.value()->getPosition().getY(),iter.value()->getPosition().getZ(),iter.value()->getData());
                //	  std::cout << "  " <<iter.value()->globalIndex()  << "   " << (int) iter.value()->mortonIndex() << std::endl;
                std::cout << "  "  << iter.data().getData()
                          << std::setprecision(PrintPrecision)
                          << std::setw(15) << iter.data().getPosition().getX()
                          << std::setw(15) << iter.data().getPosition().getY()
                          << std::setw(15) <<iter.data().getPosition().getZ()
                          << std::endl;
                iter.gotoNext();
            }

        } while(octreeIterator.moveRight());
    }
    


    return 0;
}



