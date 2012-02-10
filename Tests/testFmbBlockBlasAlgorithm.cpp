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

#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"
#include "../Src/Components/FBasicCell.hpp"

#include "../Src/Kernels/FSphericalBlockBlasKernel.hpp"
#include "../Src/Kernels/FSphericalParticle.hpp"

#include "../Src/Files/FFmaScanfLoader.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


class PointingSphericalCell : public FBasicCell {
protected:
    static int DevP;
    static int LocalSize;
    static int PoleSize;
    static bool UseBlas;

    FComplexe* multipole_exp; //< For multipole extenssion
    FComplexe* local_exp;     //< For local extenssion

public:
    static void Init(const int inDevP, const bool inUseBlas = false){
        DevP  = inDevP;
        const int ExpP  = int((inDevP+1) * (inDevP+2) * 0.5);
        const int NExpP = (inDevP+1) * (inDevP+1);

        LocalSize = ExpP;
        if(inUseBlas) {
            PoleSize = NExpP;
        }
        else{
            PoleSize = ExpP;
        }
    }

    static int GetLocalSize(){
        return LocalSize;
    }

    static int GetPoleSize(){
        return PoleSize;
    }

    /** Default constructor */
    PointingSphericalCell()
        : multipole_exp(0), local_exp(0){
        //multipole_exp = new FComplexe[PoleSize];
        local_exp = new FComplexe[LocalSize];
    }

    /** Default destructor */
    virtual ~PointingSphericalCell(){
        //delete[] multipole_exp;
        delete[] local_exp;
    }

    /** Copy constructor */
    PointingSphericalCell& operator=(const PointingSphericalCell& other) {
        FMemUtils::copyall(multipole_exp, other.multipole_exp, PoleSize);
        FMemUtils::copyall(local_exp, other.local_exp, LocalSize);
        return *this;
    }

    /** Set Multipole adresse */
    const void setMultipole(FComplexe*const inMultipole){
        multipole_exp = inMultipole;
    }

    /** Get Multipole */
    const FComplexe* getMultipole() const {
        return multipole_exp;
    }
    /** Get Local */
    const FComplexe* getLocal() const {
        return local_exp;
    }

    /** Get Multipole */
    FComplexe* getMultipole() {
        return multipole_exp;
    }
    /** Get Local */
    FComplexe* getLocal() {
        return local_exp;
    }
};

int PointingSphericalCell::DevP(-1);
int PointingSphericalCell::LocalSize(-1);
int PointingSphericalCell::PoleSize(-1);



// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FSphericalParticle             ParticleClass;
    typedef PointingSphericalCell                 CellClass;
    typedef FVector<ParticleClass>         ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalBlockBlasKernel<ParticleClass, CellClass, ContainerClass > KernelClass;

    typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test fmb algorithm.\n";
    //////////////////////////////////////////////////////////////
    const int DevP = FParameters::getValue(argc,argv,"-p", 8);
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;

    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
    std::cout << "Opening : " << filename << "\n";

    FFmaScanfLoader<ParticleClass> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------
    CellClass::Init(DevP);
    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    loader.fillTree(tree);

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // ---------------------------------------------

    int cellsPerLevel[NbLevels];
    tree.getNbCellsPerLevel(cellsPerLevel);
    FComplexe* matrix[NbLevels];
    memset(matrix, 0, sizeof(FComplexe*) * NbLevels);

    {
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        OctreeClass::Iterator avoidGoLeft(octreeIterator);

        for(int idxLevel = NbLevels - 1 ; idxLevel > 0; --idxLevel ){
            matrix[idxLevel] = new FComplexe[CellClass::GetPoleSize() * cellsPerLevel[idxLevel]];
            int poleIndex = 0;
            do{
                octreeIterator.getCurrentCell()->setMultipole(&matrix[idxLevel][poleIndex]);
                poleIndex += CellClass::GetPoleSize();
            } while(octreeIterator.moveRight());
            avoidGoLeft.moveUp();
            octreeIterator = avoidGoLeft;
        }
    }

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels(DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    FmmClass algo(&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    // ---------------------------------------------

    FMemUtils::DeleteAll(matrix, NbLevels);

    // ---------------------------------------------

    { // get sum forces&potential
        FReal potential = 0;
        F3DPosition forces;
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){
                potential += iter.data().getPotential() * iter.data().getPhysicalValue();
                forces += iter.data().getForces();

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());

        std::cout << "Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
        std::cout << "Potential = " << potential << std::endl;
    }

    return 0;
}



