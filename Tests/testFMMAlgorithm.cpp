// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Sources/Utils/FTic.hpp"

#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Utils/F3DPosition.hpp"

#include "../Sources/Core/FTestParticule.hpp"
#include "../Sources/Core/FTestCell.hpp"
#include "../Sources/Core/FTestKernels.hpp"

#include "../Sources/Core/FFMMAlgorithm.hpp"
#include "../Sources/Core/FFMMAlgorithmThreaded.hpp"

// Compile by : g++ testFMMAlgorithm.cpp ../Sources/Utils/FAssertable.cpp ../Sources/Utils/FDebug.cpp -lgomp -fopenmp -O2 -o testFMMAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particules is impacted each other particules
  */


// Simply create particules and try the kernels
int main(int , char ** ){
        const int NbLevels = 10;//10;
        const int SizeSubLevels = 3;//3
        const long NbPart = 2000000;//2000000
        FTestParticule* particules = new FTestParticule[NbPart];
        FTic counter;

        srand ( 1 ); // volontary set seed to constant

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        std::cout << "Creating " << NbPart << " particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particules[idxPart].setPosition(FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX,FReal(rand())/RAND_MAX);
        }

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        FOctree<FTestParticule, FTestCell, NbLevels, SizeSubLevels> tree(1.0,F3DPosition(0.5,0.5,0.5));

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        std::cout << "Inserting particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            tree.insert(&particules[idxPart]);
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        std::cout << "Working on particules ..." << std::endl;
        counter.tic();

        // FTestKernels FBasicKernels
        FTestKernels<FTestParticule, FTestCell> kernels;
        //FFMMAlgorithm FFMMAlgorithmThreaded
        FFMMAlgorithmThreaded<FTestKernels, FTestParticule, FTestCell, NbLevels, SizeSubLevels> algo(&tree,&kernels);
        algo.execute();

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        std::cout << "Check Result\n";
        { // Check that each particule has been summed with all other
            FOctree<FTestParticule, FTestCell, NbLevels, SizeSubLevels>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            do{
                if(octreeIterator.getCurrentCell()->getDataUp() != octreeIterator.getCurrentList()->getSize() ){
                        std::cout << "Problem P2M : " << (octreeIterator.getCurrentCell()->getDataUp() - octreeIterator.getCurrentList()->getSize()) << "\n";
                }
            } while(octreeIterator.moveRight());
        }
        { // Ceck if there is number of NbPart summed at level 1
            FOctree<FTestParticule, FTestCell, NbLevels, SizeSubLevels>::Iterator octreeIterator(&tree);
            octreeIterator.moveDown();
            long res = 0;
            do{
                res += octreeIterator.getCurrentCell()->getDataUp();
            } while(octreeIterator.moveRight());
            if(res != NbPart){
                std::cout << "Problem M2M at level 1 : " << res << "\n";
            }
        }
        { // Ceck if there is number of NbPart summed at level 1
            FOctree<FTestParticule, FTestCell, NbLevels, SizeSubLevels>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel ){
                long res = 0;
                do{
                    res += octreeIterator.getCurrentCell()->getDataUp();
                } while(octreeIterator.moveRight());
                if(res != NbPart){
                    std::cout << "Problem M2M at level " << idxLevel << " : " << res << "\n";
                }
                octreeIterator.moveUp();
                octreeIterator.gotoLeft();
            }
        }
        { // Check that each particule has been summed with all other
            FOctree<FTestParticule, FTestCell, NbLevels, SizeSubLevels>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            do{
                FList<FTestParticule*>::BasicIterator iter(*octreeIterator.getCurrentList());
                while( iter.isValide() ){
                    // If a particules has been impacted by less than NbPart - 1 (the current particule)
                    // there is a problem
                    if(iter.value()->getDataDown() != NbPart - 1){
                        std::cout << "Problem L2P + P2P : " << iter.value()->getDataDown() << "\n";
                    }
                    iter.progress();
                }
            } while(octreeIterator.moveRight());
        }
        std::cout << "Done\n";

        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Deleting particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particules[idxPart].~FTestParticule();
        }
        delete [] particules;
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        //////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////

        return 0;
}


// [--LICENSE--]
