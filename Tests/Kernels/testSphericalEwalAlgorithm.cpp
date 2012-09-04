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
#include <iomanip>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithmPeriodic.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"

#include "../../Src/Files/FEwalLoader.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

/** Ewal particle is used in the gadget program
  * here we try to make the same simulation
  */
class EwalParticle : public FSphericalParticle {
public:
    // Type of particle
    enum Type{
        OW,
        HW,
        Undefined,
    };

private:
    Type type; //< current type
    int index; //< current index in array

public:
    // Basic constructor
    EwalParticle() : type(Undefined), index(-1) {
    }

    Type getType() const{
        return type;
    }

    void setType(const Type inType) {
        type = inType;
    }

    int getIndex() const{
        return index;
    }

    void setIndex( const int inIndex ){
        index = inIndex;
    }
};


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef EwalParticle            ParticleClass;
    typedef FSphericalCell          CellClass;
    typedef FVector<ParticleClass>  ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >   KernelClass;

    typedef FFmmAlgorithmPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels      = FParameters::getValue(argc,argv,"-h", 4);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 2);
    const int DevP          = FParameters::getValue(argc,argv,"-P", 9);
    const int PeriodicDeep      = FParameters::getValue(argc,argv,"-per", 2);
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/testEwal417.dt");
    FTic counter;
    const FReal coeff_MD= 138935.4835 / 418.4 ;

    // -----------------------------------------------------

    std::cout << "Opening : " << filename << "\n";
    FEwalLoader<ParticleClass> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    CellClass::Init(DevP);
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tWidth : " << loader.getBoxWidth() << " \t center x : " << loader.getCenterOfBox().getX()
              << " y : " << loader.getCenterOfBox().getY() << " z : " << loader.getCenterOfBox().getZ() << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    counter.tic();

    ParticleClass * const particles = new ParticleClass[loader.getNumberOfParticles()];

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(particles[idxPart]);
        // reset forces and insert in the tree
        particles[idxPart].setIndex(idxPart);
        ParticleClass part = particles[idxPart];
        part.setForces(0,0,0);
        part.setPotential(0);
        tree.insert(part);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Create kernel & run simu ..." << std::endl;
    counter.tic();

    FmmClass algo(&tree,PeriodicDeep);
    KernelClass kernels( DevP, algo.extendedTreeHeight(), algo.extendedBoxWidth(),algo.extendedBoxCenter());
    algo.setKernel(&kernels);
    algo.execute();

    counter.tac();

    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------
    {
        FReal  potential = 0;
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){

                const ParticleClass& part = particles[iter.data().getIndex()];
                // std::cout << ">> index " << iter.data().getIndex() << " type " << iter.data().getType() << std::endl;
                // std::cout << "Good x " << part.getPosition().getX() << " y " << part.getPosition().getY() << " z " << part.getPosition().getZ() << std::endl;
                // std::cout << "FMM  x " << iter.data().getPosition().getX() << " y " << iter.data().getPosition().getY() << " z " << iter.data().getPosition().getZ() << std::endl;
                // std::cout << "Good fx " <<part.getForces().getX() << " fy " << part.getForces().getY() << " fz " << part.getForces().getZ() << std::endl;
                // std::cout << "FMM  fx " << iter.data().getForces().getX() << " fy " << iter.data().getForces().getY() << " fz " << iter.data().getForces().getZ() << std::endl;
                // std::cout << "GOOD physical value " << part.getPhysicalValue() << " potential " << part.getPotential() << std::endl;
                // std::cout << "FMM  physical value " << iter.data().getPhysicalValue() << " potential " << iter.data().getPotential() << std::endl;
                // std::cout << "\n";
                potential += iter.data().getPotential() * iter.data().getPhysicalValue();

		std::cout << " " << iter.data().getIndex()+1 << " \t "<<  std::setprecision(5)<< iter.data().getPosition().getX() << "  \t" << iter.data().getPosition().getY() << "  \t" << iter.data().getPosition().getZ() << "   Forces: \t"<< std::setprecision(8)<< iter.data().getForces().getX() << "  \t " << iter.data().getForces().getY() << "  \t " << iter.data().getForces().getZ()<< std::endl;

                potentialDiff.add(part.getPotential(),iter.data().getPotential());
                fx.add(part.getForces().getX(),iter.data().getForces().getX());
                fy.add(part.getForces().getY(),iter.data().getForces().getY());
                fz.add(part.getForces().getZ(),iter.data().getForces().getZ());

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());
	
        printf("Potential diff is = \n");
        printf("%f\n",potentialDiff.getL2Norm());
        printf("%f\n",potentialDiff.getInfNorm());
        printf("Fx diff is = \n");
        printf("%f\n",fx.getL2Norm());
        printf("%f\n",fx.getInfNorm());
        printf("Fy diff is = \n");
        printf("%f\n",fy.getL2Norm());
        printf("%f\n",fy.getInfNorm());
        printf("Fz diff is = \n");
        printf("%f\n",fz.getL2Norm());
        printf("%f\n",fz.getInfNorm());
        std::cout << std::endl<< std::endl<< "Potential = " << potential*coeff_MD << std::endl;

    }

    // -----------------------------------------------------

    delete[] particles;

    return 0;
}



