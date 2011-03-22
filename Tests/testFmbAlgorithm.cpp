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

#include "../Sources/Core/FFmaParticule.hpp"
#include "../Sources/Core/FBasicCell.hpp"

#include "../Sources/Core/FFMMAlgorithm.hpp"
#include "../Sources/Core/FFmbKernels.hpp"

#include "../Sources/Files/FFMALoader.hpp"

// Compile by : g++ testFmbAlgorithm.cpp ../Sources/Utils/FAssertable.cpp ../Sources/Utils/FDebug.cpp -O2 -o testFmbAlgorithm.exe

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particules is little or longer
  * related that each other
  */


/** Custom particule class */
class FmbParticule : public FFmaParticule {
public:
    FmbParticule(const double inX, const double inY, const double inZ, const double inValue)
        : FFmaParticule(inX,inY,inZ,inValue) {
    }
    FmbParticule(const F3DPosition& inPos) : FFmaParticule(inPos) {
    }
    FmbParticule(){
    }
};

/** Custom cell */
class FmbCell : public FBasicCell {
    static const int FMB_Info_P = 2;
    static const int MultipoleSize = int(((FMB_Info_P)+1) * ((FMB_Info_P)+2) * 0.5);
    FComplexe multipole_exp[MultipoleSize];
    FComplexe local_exp[MultipoleSize];
public:
    FmbCell() {
        for(int idxPole = 0 ; idxPole < MultipoleSize ; ++idxPole){
            this->multipole_exp[idxPole].setImag(0);
            this->multipole_exp[idxPole].setReal(0);
            this->local_exp[idxPole].setImag(0);
            this->local_exp[idxPole].setReal(0);
        }
    }
    FComplexe* getMultipole() {
        return this->multipole_exp;
    }
    FComplexe* getLocal() {
        return this->local_exp;
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
        const int NbLevels = 9;//10;
        const int NbSubLevels = 3;//3
        FTic counter;
        const char* const filename = "testFMAlgorithm.fma";

        FFMALoader<FmbParticule> loader(filename);
        if(!loader.isValide()){
            std::cout << "Loader Error, " << filename << "is missing\n";
            return 1;
        }
        // -----------------------------------------------------

        FOctree<FmbParticule, FmbCell, NbLevels, NbSubLevels> tree(loader.getBoxWidth(),loader.getCenterOfBox());

        // -----------------------------------------------------
        std::cout << "Creating " << loader.getNumberOfParticules() << " particules ..." << std::endl;
        counter.tic();

        FmbParticule* particules = new FmbParticule[loader.getNumberOfParticules()];

        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticules() ; ++idxPart){
            loader.fillParticule(&particules[idxPart]);
        }

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        // -----------------------------------------------------
        std::cout << "Inserting particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < loader.getNumberOfParticules() ; ++idxPart){
            tree.insert(&particules[idxPart]);
        }
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

        std::cout << "Working on particules ..." << std::endl;
        counter.tic();

        FFmbKernels<FmbParticule, FmbCell> kernels(NbLevels,loader.getBoxWidth());
        FFMMAlgorithm<FmbParticule, FmbCell, NbLevels, NbSubLevels> algo(&tree,&kernels);
        algo.execute();

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        // -----------------------------------------------------
        std::cout << "Deleting particules ..." << std::endl;
        counter.tic();
        for(long idxPart = 0 ; idxPart < loader.getNumberOfParticules() ; ++idxPart){
            particules[idxPart].~FmbParticule();
        }
        delete [] particules;
        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;
        // -----------------------------------------------------

	return 0;
}


// [--LICENSE--]
