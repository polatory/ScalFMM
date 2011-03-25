// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Sources/Utils/FTic.hpp"

#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Utils/FAssertable.hpp"
#include "../Sources/Utils/F3DPosition.hpp"

#include "../Sources/Core/FBasicParticule.hpp"
#include "../Sources/Core/FBasicCell.hpp"

#include "../Sources/Core/FFMMAlgorithm.hpp"
#include "../Sources/Core/FAbstractKernels.hpp"

// Compile by : g++ testFMMAlgorithm.cpp ../Sources/Utils/FAssertable.cpp ../Sources/Utils/FDebug.cpp -O2 -o testFMMAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particules is little or longer
  * related that each other
  */


/** Custom particule class */
class MyTestParticule : public FBasicParticule{
protected:
    // To store data during downard pass
    long dataDown;
public:
    MyTestParticule(): dataDown(0){
    }
    long getDataDown() const {
        return this->dataDown;
    }
    void setDataDown(const long inData){
        this->dataDown = inData;
    }
};

/** Custom cell */
class MyTestCell : public FBasicCell {
    // To store data during upward and downward pass
    long dataUp, dataDown;
public:
    MyTestCell(): dataUp(0) , dataDown(0){
    }
    long getDataUp() const {
        return this->dataUp;
    }
    void setDataUp(const long inData){
        this->dataUp = inData;
    }
    long getDataDown() const {
        return this->dataDown;
    }
    void setDataDown(const long inData){
        this->dataDown = inData;
    }
};

/** Custom Kernel
  * This kernel simply store in each element the number of elements
  * it represents
  */
template< class ParticuleClass, class CellClass>
class MyTestKernels : public FAbstractKernels<ParticuleClass,CellClass> {
public:
    // Before upward
    virtual void P2M(CellClass* const pole, FList<ParticuleClass*>* const particules) {
        // the pole represents all particules under
        pole->setDataUp(particules->getSize());
    }
    // During upward
    virtual void M2M(CellClass* const pole, CellClass** const child, const int inLevel) {
        // A parent represents the sum of the child
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child[idx]){
                pole->setDataUp(pole->getDataUp() + child[idx]->getDataUp());
            }
        }
    }
    // Before Downward
    virtual void M2L(CellClass* const pole, CellClass** const distantNeighbors, const int size, const int inLevel) {
        // The pole is impacted by what represent other poles
        for(int idx = 0 ; idx < size ; ++idx){
            pole->setDataDown(pole->getDataDown() + distantNeighbors[idx]->getDataUp());
        }
    }
    // During Downward
    virtual void L2L(CellClass* const local, CellClass** const child, const int inLevel) {
        // Each child is impacted by the father
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child[idx]){
                child[idx]->setDataDown(local->getDataDown() + child[idx]->getDataDown());
            }
        }
    }
    // After Downward
    virtual void L2P(CellClass* const local, FList<ParticuleClass*>* const particules){
        // The particules is impacted by the parent cell
        typename FList<ParticuleClass*>::BasicIterator iter(*particules);
        while( iter.isValide() ){
            iter.value()->setDataDown(iter.value()->getDataDown() + local->getDataDown());
            iter.progress();
        }
    }
    // After Downward
    virtual void P2P(FList<ParticuleClass*>* const currentBox, FList<ParticuleClass*>** directNeighbors, const int size) {
        // Each particules targeted is impacted by the particules sources
        long inc = currentBox->getSize() - 1;
        for(int idx = 0 ; idx < size ; ++idx){
            inc += directNeighbors[idx]->getSize();
        }

        typename FList<ParticuleClass*>::BasicIterator iter(*currentBox);
        while( iter.isValide() ){
            iter.value()->setDataDown(iter.value()->getDataDown() + inc);
            iter.progress();
        }
    }
};

// Simply create particules and try the kernels
int main(int , char ** ){
        const int NbLevels = 10;//10;
        const int SizeSubLevels = 3;//3
        const long NbPart = 2000000;//2E6;
        MyTestParticule* particules = new MyTestParticule[NbPart];
        FTic counter;

        srand ( 1 ); // volontary set seed to constant

        // -----------------------------------------------------
        std::cout << "Creating " << NbPart << " particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particules[idxPart].setPosition(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
        }

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

        FOctree<MyTestParticule, MyTestCell, NbLevels, SizeSubLevels> tree(1.0,F3DPosition(0.5,0.5,0.5));

        // -----------------------------------------------------
        std::cout << "Inserting particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            tree.insert(&particules[idxPart]);
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

        std::cout << "Working on particules ..." << std::endl;
        counter.tic();

        MyTestKernels<MyTestParticule, MyTestCell> kernels;
        FFMMAlgorithm<MyTestKernels, MyTestParticule, MyTestCell, NbLevels, SizeSubLevels> algo(&tree,&kernels);
        algo.execute();

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        // -----------------------------------------------------

        std::cout << "Check Result\n";
        { // Check that each particule has been summed with all other
            FOctree<MyTestParticule, MyTestCell, NbLevels, SizeSubLevels>::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            do{
                FList<MyTestParticule*>::BasicIterator iter(*octreeIterator.getCurrentList());
                while( iter.isValide() ){
                    // If a particules has been impacted by less than NbPart - 1 (the current particule)
                    // there is a problem
                    if(iter.value()->getDataDown() != NbPart - 1){
                        std::cout << "Problem : " << iter.value()->getDataDown() << "\n";
                    }
                    iter.progress();
                }
            } while(octreeIterator.moveRight());
        }
        { // Ceck if there is number of NbPart summed at level 1
            FOctree<MyTestParticule, MyTestCell, NbLevels, SizeSubLevels>::Iterator octreeIterator(&tree);
            octreeIterator.moveDown();
            long res = 0;
            do{
                res += octreeIterator.getCurrentCell()->getDataUp();
            } while(octreeIterator.moveRight());
            if(res != NbPart){
                std::cout << "Problem at level 1 : " << res << "\n";
            }
        }
        std::cout << "Done\n";

        // -----------------------------------------------------
        std::cout << "Deleting particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particules[idxPart].~MyTestParticule();
        }
        delete [] particules;
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

	return 0;
}


// [--LICENSE--]
