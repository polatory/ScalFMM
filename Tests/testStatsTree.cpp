// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Components/FFmaParticle.hpp"
#include "../Src/Extensions/FExtendForces.hpp"
#include "../Src/Extensions/FExtendPotential.hpp"

#include "../Src/Components/FBasicCell.hpp"
#include "../Src/Fmb/FExtendFmbCell.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Fmb/FFmbKernels.hpp"

#include "../Src/Utils/FMath.hpp"

#include "../Src/Files/FFmaLoader.hpp"


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


    FFmaParticle particle;
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(particle);
        tree.insert(particle);
    }

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
                typename OctreeClass::Iterator octreeIterator(&tree);
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
                typename OctreeClass::Iterator octreeIterator(&tree);
                octreeIterator.gotoBottomLeft();
                do{
                    varianceParticles += (averageParticles - FReal(octreeIterator.getCurrentListTargets()->getSize())) - (averageParticles - FReal(octreeIterator.getCurrentListTargets()->getSize()));
                    ++nbLeafs;
                } while(octreeIterator.moveRight());
                varianceParticles /= FReal(nbLeafs);
            }
            std::cout << "[STAT] Variances of particles on leafs is = " << FMath::Sqrt(varianceParticles) << std::endl;
        }


        {
            FReal averageReduction = 0;
            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            int previousCells = 0;

            do{
                ++previousCells;
            } while(octreeIterator.moveRight());

            for(int idxLevel = NbLevels - 2 ; idxLevel > 1 ; --idxLevel){
                octreeIterator.moveUp();
                octreeIterator.gotoLeft();

                int nbCells = 0;
                do{
                    ++nbCells;
                } while(octreeIterator.moveRight());
                averageReduction += FReal(nbCells)/FReal(previousCells);
                previousCells = nbCells;
            }
            std::cout << "[STAT] Average reduction by level = " << (averageReduction/FReal(NbLevels-2)) << std::endl;
        }

        {
            FReal averageNeighbors = 0;
            int nbLeafs = 0;
            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                ContainerClass* neighbors[26];
                // need the current particles and neighbors particles
                averageNeighbors += FReal(tree.getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalIndex(),NbLevels-1));
                ++nbLeafs;
            } while(octreeIterator.moveRight());
            std::cout << "[STAT] Average neighbors for each leafs = " << (averageNeighbors/FReal(nbLeafs)) << std::endl;
        }

        {
            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            int fullNbCells = 0;
            int fullNbLeafs = 0;
            int fullNbChild = 0;
            int fullM2LCalculus = 0;

            int nbCellsAtLevel[NbLevels - 3];
            int nbChildAtLevel[NbLevels - 3];
            int M2LCalculusAtLevel[NbLevels - 3];

            // Get data
            do{
                ++fullNbLeafs;
            } while(octreeIterator.moveRight());

            for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel){
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

                    const FBasicCell* neighbors[189];

                    M2LCalculusAtLevel[idxArray] += tree.getDistantNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),idxLevel);

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
            for(int idxLevel = NbLevels - 2 ; idxLevel > 1 ; --idxLevel){
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
            for(int idxLevel = NbLevels - 2 ; idxLevel > 1 ; --idxLevel){
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



