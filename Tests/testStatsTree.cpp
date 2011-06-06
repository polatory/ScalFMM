// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FList.hpp"

#include "../Src/Components/FFmaParticle.hpp"
#include "../Src/Extensions/FExtendForces.hpp"
#include "../Src/Extensions/FExtendPotential.hpp"

#include "../Src/Components/FBasicCell.hpp"
#include "../Src/Fmb/FExtendFmbCell.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"
#include "../Src/Core/FFmmAlgorithmThreadUs.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Fmb/FFmbKernels.hpp"

#include "../Src/Utils/FMath.hpp"

#include "../Src/Files/FFmaLoader.hpp"

// With openmp : g++ testFmbAlgorithm.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmbAlgorithm.exe
// icpc -openmp -openmp-lib=compat testFmbAlgorithm.cpp ../Src/Utils/FAssertable.cpp ../Src/Utils/FDebug.cpp -O2 -o testFmbAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to show some stat about the tree.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = 7;//10;
    const int SizeSubLevels = 3;//3
    FTic counter;
    const char* const defaultFilename = "testLoaderFMA.fma"; //../../Data/ "testLoaderFMA.fma" "testFMAlgorithm.fma" Sphere.fma
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

    FFmaLoader<FFmaParticle> loader(filename);
    if(!loader.hasNotFinished()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    FOctree<FFmaParticle, FBasicCell, FSimpleLeaf>
            tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating and Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    counter.tic();


    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FFmaParticle particle;

        loader.fillParticle(particle);

        tree.insert(particle);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;


    // -----------------------------------------------------

    { // get stats
        {
            std::cout << "[STAT] level is " << NbLevels << std::endl;
            std::cout << "[STAT] potentials leafs number is " << (1 << (3*NbLevels)) << std::endl;

            FReal averageParticles = 0;
            {
                long nbLeafs = 0;
                FOctree<FFmaParticle, FBasicCell, FSimpleLeaf>::Iterator octreeIterator(&tree);
                octreeIterator.gotoBottomLeft();
                do{
                    averageParticles += octreeIterator.getCurrentListTargets()->getSize();
                    ++nbLeafs;
                } while(octreeIterator.moveRight());
                averageParticles /= nbLeafs;

                std::cout << "[STAT] Nb leafs : " << nbLeafs << std::endl;
            }
            std::cout << "[STAT] Average particles on leafs = " << averageParticles << std::endl;

            FReal varianceParticles = 0;
            {
                long nbLeafs = 0;
                FOctree<FFmaParticle, FBasicCell, FSimpleLeaf>::Iterator octreeIterator(&tree);
                octreeIterator.gotoBottomLeft();
                do{
                    varianceParticles += (averageParticles - octreeIterator.getCurrentListTargets()->getSize()) - (averageParticles - octreeIterator.getCurrentListTargets()->getSize());
                    ++nbLeafs;
                } while(octreeIterator.moveRight());
                varianceParticles /= nbLeafs;
            }
            std::cout << "[STAT] Variances of particles on leafs is = " << FMath::Sqrt(varianceParticles) << std::endl;
        }


        {
            FReal averageReduction = 0;
            FOctree<FFmaParticle, FBasicCell, FSimpleLeaf>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            long previousCells = 0;

            do{
                ++previousCells;
            } while(octreeIterator.moveRight());

            for(long idxLevel = NbLevels - 2 ; idxLevel > 1 ; --idxLevel){
                octreeIterator.moveUp();
                octreeIterator.gotoLeft();

                long nbCells = 0;
                do{
                    ++nbCells;
                } while(octreeIterator.moveRight());
                averageReduction += nbCells/float(previousCells);
                previousCells = nbCells;
            }
            std::cout << "[STAT] Average reduction by level = " << (averageReduction/(NbLevels-2)) << std::endl;
        }

        {
            FReal averageNeighbors = 0;
            long nbLeafs = 0;
            FOctree<FFmaParticle, FBasicCell, FSimpleLeaf>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                FList<FFmaParticle>* neighbors[26];
                // need the current particles and neighbors particles
                const long counterNeighbors = tree.getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalIndex(),NbLevels-1);
                ++nbLeafs;
                for(long idx = 0 ; idx < counterNeighbors ; ++idx){
                    averageNeighbors += neighbors[idx]->getSize();
                }
            } while(octreeIterator.moveRight());
            std::cout << "[STAT] Average neighbors for each leafs = " << (averageNeighbors/nbLeafs) << std::endl;
        }

        {
            FOctree<FFmaParticle, FBasicCell, FSimpleLeaf>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            long fullNbCells = 0;
            long fullNbLeafs = 0;
            long fullNbChild = 0;
            long fullM2LCalculus = 0;

            long nbCellsAtLevel[NbLevels - 3];
            long nbChildAtLevel[NbLevels - 3];
            long M2LCalculusAtLevel[NbLevels - 3];

            // Get data
            do{
                ++fullNbLeafs;
            } while(octreeIterator.moveRight());

            for(long idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel){
                octreeIterator.moveUp();
                octreeIterator.gotoLeft();

                const int idxArray = idxLevel - 2;
                nbChildAtLevel[idxArray] = 0;
                nbCellsAtLevel[idxArray] = 0;
                M2LCalculusAtLevel[idxArray] = 0;

                do{
                    ++nbCellsAtLevel[idxArray];

                    FBasicCell** child = octreeIterator.getCurrentChild();
                    for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                        if(child[idxChild]) ++nbChildAtLevel[idxArray];
                    }

                    const FBasicCell* neighbors[208];
                    MortonIndex neighborsIndexes[208];
                    FTreeCoordinate currentPosition;
                    FTreeCoordinate neighborsPosition[208];

                    M2LCalculusAtLevel[idxArray] += tree.getDistantNeighborsWithIndex(neighbors, neighborsIndexes, currentPosition, neighborsPosition, octreeIterator.getCurrentGlobalIndex(),idxLevel);

                } while(octreeIterator.moveRight());

                fullNbCells += nbCellsAtLevel[idxArray];
                fullNbChild += nbChildAtLevel[idxArray];
                fullM2LCalculus += M2LCalculusAtLevel[idxArray];
            }


            // compute basic stats
            {
                const double averageChildPerCell = (fullNbChild/double(fullNbCells));
                std::cout << "[STAT] Average number of children = " << averageChildPerCell << std::endl;

                const double averageM2LCalculusPerCell = (fullM2LCalculus/double(fullNbCells));
                std::cout << "[STAT] Average M2L operations = " << averageM2LCalculusPerCell << std::endl;
            }
            // compute by scoring level by level
            double averageChildByLevel = 0;
            double averageM2LByLevel = 0;
            for(long idxLevel = NbLevels - 2 ; idxLevel > 1 ; --idxLevel){
                const int idxArray = idxLevel - 2;
                averageChildByLevel += nbChildAtLevel[idxArray]/double(nbCellsAtLevel[idxArray]);
                averageM2LByLevel += M2LCalculusAtLevel[idxArray]/double(nbCellsAtLevel[idxArray]);
            }
            averageChildByLevel /= NbLevels - 3;
            averageM2LByLevel /= NbLevels - 3;
            std::cout << "[STAT] Average of children per level = " << averageChildByLevel << std::endl;
            std::cout << "[STAT] Average M2L per level = " << averageM2LByLevel << std::endl;

            // doing a variance for theses data
            double varianceChildByLevel = 0;
            double varianceM2LByLevel = 0;
            for(long idxLevel = NbLevels - 2 ; idxLevel > 1 ; --idxLevel){
                const int idxArray = idxLevel - 2;

                const double averageAtLevel = nbChildAtLevel[idxArray]/double(nbCellsAtLevel[idxArray]);
                varianceChildByLevel += (averageAtLevel - averageChildByLevel) * (averageAtLevel - averageChildByLevel);

                const double averageAtLevelM2L = M2LCalculusAtLevel[idxArray]/double(nbCellsAtLevel[idxArray]);
                varianceM2LByLevel += (averageAtLevelM2L - averageM2LByLevel) * (averageAtLevelM2L - averageM2LByLevel);
            }
            std::cout << "[STAT] Variance number of children per level = " << FMath::Sqrt(varianceChildByLevel/double(NbLevels - 3)) << std::endl;
            std::cout << "[STAT] Variance M2L per level = " << FMath::Sqrt(varianceM2LByLevel/double(NbLevels - 3)) << std::endl;
        }
    }

    // -----------------------------------------------------


    return 0;
}


// [--LICENSE--]
