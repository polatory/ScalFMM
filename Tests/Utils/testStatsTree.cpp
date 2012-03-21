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

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FFmaParticle.hpp"
#include "../../Src/Extensions/FExtendForces.hpp"
#include "../../Src/Extensions/FExtendPotential.hpp"

#include "../../Src/Components/FBasicCell.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"


#include "../../Src/Utils/FMath.hpp"

#include "../../Src/Files/FFmaLoader.hpp"


/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FVector<FFmaParticle>      ContainerClass;
    typedef FSimpleLeaf<FFmaParticle, ContainerClass >                     LeafClass;
    typedef FOctree<FFmaParticle, FBasicCell, ContainerClass , LeafClass >  OctreeClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to show some stat about the tree.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;

    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
    std::cout << "Opening : " << filename << "\n";

    FFmaLoader<FFmaParticle> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating and Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();


    loader.fillTree(tree);

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;


    // -----------------------------------------------------

    { // get stats
        {
            std::cout << "[STAT] level is " << NbLevels << std::endl;
            std::cout << "[STAT] potentials leafs number is " << (1 << (3* (NbLevels-1) )) << std::endl;

            FReal averageParticles = 0;
            {
                int nbLeafs = 0;
                OctreeClass::Iterator octreeIterator(&tree);
                octreeIterator.gotoBottomLeft();
                do{
                    averageParticles += FReal(octreeIterator.getCurrentListTargets()->getSize());
                    ++nbLeafs;
                } while(octreeIterator.moveRight());
                averageParticles /= FReal(nbLeafs);

                std::cout << "[STAT] Nb leafs : " << nbLeafs << std::endl;
            }
            std::cout << "[STAT] Average particles on leafs = " << averageParticles << std::endl;

            FReal varianceParticles = 0;
            {
                int nbLeafs = 0;
                OctreeClass::Iterator octreeIterator(&tree);
                octreeIterator.gotoBottomLeft();
                do{
                    varianceParticles += FReal(octreeIterator.getCurrentListTargets()->getSize() * octreeIterator.getCurrentListTargets()->getSize());
                    ++nbLeafs;
                } while(octreeIterator.moveRight());
                varianceParticles /= FReal(nbLeafs);
            }
            std::cout << "[STAT] Variances of particles on leafs is = " << (varianceParticles - (averageParticles*averageParticles)) << std::endl;
        }

        {
            FReal averageNeighbors = 0;
            int nbLeafs = 0;
            OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                ContainerClass* neighbors[27];
                // need the current particles and neighbors particles
                averageNeighbors += FReal(tree.getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),NbLevels-1));
                ++nbLeafs;
            } while(octreeIterator.moveRight());
            std::cout << "[STAT] Average neighbors for each leafs = " << (averageNeighbors/FReal(nbLeafs)) << std::endl;
        }

        {
            long long int totalCells = 0;
            long long int totalM2L = 0;
            long long int totalM2ML2L = 0;

            int nbCellsAtTop = 0;
            int nbCellsAtBottom = 0;

            OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel){

                int nbCellsAtLevel = 0;
                int nbChildAtLevel = 0;
                int nbNeighborsAtLevel = 0;

                do{
                    ++nbCellsAtLevel;

                    if( idxLevel != NbLevels - 1 ){
                        FBasicCell** child = octreeIterator.getCurrentChild();
                        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                            if(child[idxChild]) ++nbChildAtLevel;
                        }
                    }

                    const FBasicCell* neighbors[343];
                    nbNeighborsAtLevel += tree.getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),idxLevel);

                } while(octreeIterator.moveRight());

                octreeIterator.moveUp();
                octreeIterator.gotoLeft();

                std::cout << "[STAT] Level = " << idxLevel << "\n";
                std::cout << "[STAT] >> Nb Cells = " << nbCellsAtLevel << "\n";
                std::cout << "[STAT] >> Nb M2M/L2L interactions = " << nbChildAtLevel << "\n";
                std::cout << "[STAT] >> Average M2M/L2L interactions = " << FReal(nbChildAtLevel)/FReal(nbCellsAtLevel) << "\n";
                std::cout << "[STAT] >> Nb M2L interactions = " << nbNeighborsAtLevel << "\n";
                std::cout << "[STAT] >> Average M2L interactions = " << FReal(nbNeighborsAtLevel)/FReal(nbCellsAtLevel) << "\n";

                totalCells += (long long int)(nbCellsAtLevel);
                totalM2L += (long long int)(nbNeighborsAtLevel);
                totalM2ML2L += (long long int)(nbChildAtLevel);
                nbCellsAtTop = nbCellsAtLevel;
                if( idxLevel == NbLevels - 1 ) nbCellsAtBottom = nbCellsAtLevel;
            }

            std::cout << "[STAT] For all the tree\n";
            std::cout << "[STAT] >> Total Nb Cells = " << totalCells-nbCellsAtTop << "\n";
            std::cout << "[STAT] >> Total Nb M2M/L2L interactions = " << totalM2ML2L << "\n";
            std::cout << "[STAT] >> Total Average M2M/L2L interactions = " << FReal(totalM2ML2L-nbCellsAtTop)/FReal(totalCells-nbCellsAtBottom) << "\n";
            std::cout << "[STAT] >> Total Nb M2L interactions per cell = " << totalM2L << "\n";
            std::cout << "[STAT] >> Total Average M2L interactions per cell = " << FReal(totalM2L)/FReal(totalCells) << "\n";

       }
    }

    // -----------------------------------------------------


    return 0;
}


