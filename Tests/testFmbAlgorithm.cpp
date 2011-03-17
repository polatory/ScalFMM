// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include "../Sources/Utils/FTic.hpp"

#include "../Sources/Containers/FOctree.hpp"
#include "../Sources/Containers/FList.hpp"

#include "../Sources/Utils/FAssertable.hpp"
#include "../Sources/Utils/F3DPosition.hpp"
#include "../Sources/Utils/FComplexe.hpp"

#include "../Sources/Core/FBasicParticule.hpp"
#include "../Sources/Core/FBasicCell.hpp"

#include "../Sources/Core/FFMMAlgorithm.hpp"
#include "../Sources/Core/FFmbKernels.hpp"

// Compile by : g++ testFMMAlgorithm.cpp ../Sources/Utils/FAssertable.cpp ../Sources/Utils/FDebug.cpp -O2 -o testFMMAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particules is little or longer
  * related that each other
  */


/** Custom particule class */
class FmbParticule : public FBasicParticule{
protected:
    double value;
public:
    FmbParticule(const double inX, const double inY, const double inZ) : FBasicParticule(inX,inY,inZ) {
        setValue(0);
    }
    FmbParticule(const F3DPosition& inPos) : FBasicParticule(inPos) {
        setValue(0);
    }
    FmbParticule(){
        setValue( 0 );
    }
    double getValue() const {
        return this->value;
    }
    void setValue(const double inValue){
        this->value = inValue;
    }
};

/** Custom cell */
class FmbCell : public FBasicCell {
    static const int FMB_Info_P = 2;
    static const int MultipoleSize = int(((FMB_Info_P)+1) * ((FMB_Info_P)+2) * 0.5);
    FComplexe multipole_exp[MultipoleSize];
public:
    FmbCell() {
    }
    FComplexe* getMultipole() {
        return this->multipole_exp;
    }
};

/*
 Before putting position in cells :
Creating 2000000 particules ...
Inserting particules ...
Done  (1.82029s).

 After putting position in cells :
Creating 2000000 particules ...
Inserting particules ...
Done  (1.94356s).
 */

// Simply create particules and try the kernels
int main(int , char ** ){
        const int NbLevels = 10;//10;
        const int NbSubLevels = 3;//3
        const long NbPart = 2;//2E6;
        const int BoxWidth = 2;//1
        FmbParticule* particules = new FmbParticule[NbPart];
        FTic counter;

        srand ( 1 ); // volontary set seed to constant

        // -----------------------------------------------------
        std::cout << "Creating " << NbPart << " particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            new (&particules[idxPart]) FmbParticule(double(rand())/RAND_MAX,double(rand())/RAND_MAX,double(rand())/RAND_MAX);
        }

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

        FOctree<FmbParticule, FmbCell, NbLevels, NbSubLevels> tree(BoxWidth,F3DPosition(0.0,0.0,0.0));

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

        FFmbKernels<FmbParticule, FmbCell> kernels(NbLevels,BoxWidth);
        FFMMAlgorithm<FmbParticule, FmbCell, NbLevels, NbSubLevels> algo(&tree,&kernels);
        algo.execute();

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        // -----------------------------------------------------
        std::cout << "Deleting particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
            particules[idxPart].~FmbParticule();
        }
        delete [] particules;
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

	return 0;
}


// [--LICENSE--]
