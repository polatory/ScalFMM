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

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FBasicParticleContainer.hpp"
#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Arranger/FOctreeArranger.hpp"

struct TestParticle{
    FPoint position;
    const FPoint& getPosition(){
        return position;
    }
};

template <class ParticleClass>
class Converter {
public:
    template <class ContainerClass>
    static ParticleClass GetParticle(ContainerClass* container, const int idxExtract){
        TestParticle part;
        part.position.setPosition(
                    container->getPositions()[0][idxExtract],
                    container->getPositions()[1][idxExtract],
                    container->getPositions()[2][idxExtract]);
        return part;
    }

    template <class OctreeClass>
    static void Insert(OctreeClass* tree, const ParticleClass& part){
        tree->insert(part.position);
    }
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FBasicCell                      CellClass;
    typedef FBasicParticleContainer<0>      ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels          = FParameters::getValue(argc,argv,"-h", 7);
    const int SizeSubLevels     = FParameters::getValue(argc,argv,"-sh", 3);
    const int NbPart           = FParameters::getValue(argc,argv,"-nb", 2000000);
    const FReal FRandMax        = FReal(RAND_MAX);

    FTic counter;

    srand ( 1 ); // volontary set seed to constant

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    const FReal BoxWidth = 1.0;
    const FReal BoxCenter = 0.5;

    OctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, FPoint(BoxCenter,BoxCenter,BoxCenter));

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating & Inserting " << NbPart << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();


    {

        FPoint particleToFill;
        for(int idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particleToFill.setPosition(
                        (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)),
                        (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)),
                        (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2)));
            tree.insert(particleToFill);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    { // Create new position for each particles
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass* particles = octreeIterator.getCurrentListTargets();
            for(int idxPart = 0; idxPart < particles->getNbParticles() ; ++idxPart){
                particles->getWPositions()[0][idxPart] = (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2));
                particles->getWPositions()[1][idxPart] = (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2));
                particles->getWPositions()[2][idxPart] = (BoxWidth*FReal(rand())/FRandMax) + (BoxCenter-(BoxWidth/2));
            }
        } while(octreeIterator.moveRight());
    }

    counter.tac();
    std::cout << "Done  " << "(@Moving = " << counter.elapsed() << "s)." << std::endl;


    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Arrange ..." << std::endl;
    counter.tic();

    FOctreeArranger<OctreeClass, ContainerClass, TestParticle, Converter<TestParticle> > arrange(&tree);
    arrange.rearrange();

    counter.tac();
    std::cout << "Done  " << "(@Arrange = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Test ..." << std::endl;
    counter.tic();

    { // Check that each particle has been put into the right leaf
        long counterPart = 0;

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            const MortonIndex leafIndex = octreeIterator.getCurrentGlobalIndex();

            ContainerClass* particles = octreeIterator.getCurrentListTargets();
            for(int idxPart = 0; idxPart < particles->getNbParticles() ; ++idxPart){
                const FPoint particlePosition( particles->getWPositions()[0][idxPart],
                                               particles->getWPositions()[1][idxPart],
                                               particles->getWPositions()[2][idxPart]);

                const MortonIndex particleIndex = tree.getMortonFromPosition( particlePosition );
                if( leafIndex != particleIndex){
                    std::cout << "Index problem, should be " << leafIndex <<
                                 " particleIndex "<< particleIndex << std::endl;
                }

            }

            counterPart += octreeIterator.getCurrentListTargets()->getNbParticles();
            if(octreeIterator.getCurrentListTargets()->getNbParticles() == 0){
                std::cout << "Problem, leaf is empty at index " << leafIndex << std::endl;
            }
        } while(octreeIterator.moveRight());

        if( counterPart != NbPart ){
            std::cout <<"Wrong particles number, should be " << NbPart << " but is " << counterPart << std::endl;
        }
    }

    { // Check that each particle has been summed with all other
        OctreeClass::Iterator octreeIterator(&tree);
        OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        const int heightMinusOne = NbLevels - 1;
        for(int idxLevel = 1 ; idxLevel < heightMinusOne ; ++idxLevel ){
            // for each cells
            do{
                int countChild = 0;
                CellClass** const child = octreeIterator.getCurrentChild();
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild ){
                    if( child[idxChild] ){
                        countChild += 1;
                    }
                }

                if(countChild == 0){
                    std::cout << "Problem at level " << idxLevel << " cell has no child " << octreeIterator.getCurrentGlobalIndex() << std::endl;
                }

            } while(octreeIterator.moveRight());

            avoidGotoLeftIterator.moveDown();
            octreeIterator = avoidGotoLeftIterator;
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Test = " << counter.elapsed() << "s)." << std::endl;
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}



