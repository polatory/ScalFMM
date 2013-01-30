// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <time.h>


#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FAssertable.hpp"
#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FBasicParticle.hpp"
#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Utils/FTic.hpp"

/**
* In this file we show how to use octree with iteration
* This is a good example to understand FOctree::Iterator.
*/

int main(int argc, char ** argv){
    typedef FVector<FBasicParticle>      ContainerClass;
    typedef FSimpleLeaf<FBasicParticle, ContainerClass >                     LeafClass;
    typedef FOctree<FBasicParticle, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use octree iterator.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 9);
    const int NbSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    const int NbPart = FParameters::getValue(argc,argv,"-nb", 2000000);
    const FReal FRandMax = FReal(RAND_MAX);

    FTic counterTime;

    srand ( 1 ); // volontary set seed to constant
    // -----------------------------------------------------

    OctreeClass tree(NbLevels, NbSubLevels, 1.0, FPoint(0.5,0.5,0.5));

    // -----------------------------------------------------
    std::cout << "Creating and inserting " << NbPart << " particles ..." << std::endl;
    counterTime.tic();
    {
        FBasicParticle particle;
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particle.setPosition(FReal(rand())/FRandMax,FReal(rand())/FRandMax,FReal(rand())/FRandMax);
            tree.insert(particle);
        }
    }
    counterTime.tac();
    std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------
    {
        std::cout << "Itering on particles ..." << std::endl;
        counterTime.tic();

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        for(int idxLevel = NbLevels - 1 ; idxLevel >= 1 ; --idxLevel ){
            int counter = 0;
            do{
                ++counter;
                //counter += octreeIterator.getCurrentList()->getSize();
            } while(octreeIterator.moveRight());
            octreeIterator.moveUp();
            octreeIterator.gotoLeft();
            std::cout << "Cells at this level " << counter << " ...\n";
        }
        counterTime.tac();
        std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;
    }
    // -----------------------------------------------------
    {
        std::cout << "Itering on particles fast ..." << std::endl;
        counterTime.tic();

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        OctreeClass::Iterator avoidGoLeft(octreeIterator);

        for(int idx = 0 ; idx < NbLevels - 1; ++idx ){
            do{
            } while(octreeIterator.moveRight());
            avoidGoLeft.moveUp();
            octreeIterator = avoidGoLeft;
        }
        counterTime.tac();
        std::cout << "Done  " << "(" << counterTime.elapsed() << "s)." << std::endl;
    }
    // -----------------------------------------------------

    return 0;
}



