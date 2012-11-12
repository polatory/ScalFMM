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

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmTsm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmThreadTsm.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FTypedLeaf.hpp"


/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FSphericalParticle             ParticleClass;
    typedef FSphericalCell                 CellClass;
    typedef FVector<ParticleClass>         ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                      LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

    typedef FFmmAlgorithmThread<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

    typedef FTypedSphericalParticle             ParticleClassTyped;
    typedef FTypedSphericalCell                 CellClassTyped;
    typedef FVector<ParticleClassTyped>  ContainerClassTyped;

    typedef FTypedLeaf<ParticleClassTyped, ContainerClassTyped >                      LeafClassTyped;
    typedef FOctree<ParticleClassTyped, CellClassTyped, ContainerClassTyped , LeafClassTyped >  OctreeClassTyped;
    typedef FSphericalKernel<ParticleClassTyped, CellClassTyped, ContainerClassTyped >          KernelClassTyped;

    typedef FFmmAlgorithmThreadTsm<OctreeClassTyped, ParticleClassTyped, CellClassTyped, ContainerClassTyped, KernelClassTyped, LeafClassTyped > FmmClassTyped;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical on a Tsm system.\n";
    std::cout << ">> It compares the results between Tms and no Tms (except P2P & L2P).\n";
    //////////////////////////////////////////////////////////////
    const int DevP = FParameters::getValue(argc,argv,"-p", 8);
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const int NbPart = 200000;//2000000
    const double BoxWidth = 1.0;
    const FPoint CenterOfBox(0.5,0.5,0.5);
    const FReal FRandMax = FReal(RAND_MAX);


    // -----------------------------------------------------
    CellClass::Init(DevP);

    OctreeClass tree(NbLevels, SizeSubLevels,BoxWidth,CenterOfBox);
    OctreeClassTyped treeTyped(NbLevels, SizeSubLevels,BoxWidth,CenterOfBox);


    std::cout << "Inserting particles ..." << std::endl;
    counter.tic();
    for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
        const FReal x = FReal(rand())/FRandMax;
        const FReal y = FReal(rand())/FRandMax;
        const FReal z = FReal(rand())/FRandMax;

        ParticleClass particles;
        ParticleClassTyped particlesTyped;
        ParticleClassTyped particlesTyped2;

        // Particle for standart model
        particles.setPosition(x,y,z);
        particles.setPhysicalValue(1);

        // Create a clone for typed (Tsm) version
        particlesTyped.setPosition(x,y,z);
        particlesTyped2.setPosition(x,y,z);

        particlesTyped.setPhysicalValue(1);
        particlesTyped2.setPhysicalValue(1);

        particlesTyped.setAsSource();
        particlesTyped2.setAsTarget();

        tree.insert(particles);
        treeTyped.insert(particlesTyped);
        treeTyped.insert(particlesTyped2);
    }
    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------
    std::cout << "Create kernels ..." << std::endl;

    KernelClass kernels(DevP, NbLevels, BoxWidth, CenterOfBox);
    KernelClassTyped kernelsTyped(DevP, NbLevels, BoxWidth, CenterOfBox);


    std::cout << "Working on particles ..." << std::endl;

    FmmClass algo(&tree,&kernels);
    FmmClassTyped algoTyped(&treeTyped,&kernelsTyped);

    counter.tic();
    algo.execute();
    algoTyped.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    // Here we compare the cells of the trees that must contains the same values

    std::cout << "Start checking ..." << std::endl;
    {
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        OctreeClassTyped::Iterator octreeIteratorTyped(&treeTyped);
        octreeIteratorTyped.gotoBottomLeft();

        for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel ){
            std::cout << "\t test level " << idxLevel << "\n";

            do{
                bool poleDiff = false;
                bool localDiff = false;
                for(int idxValues = 0 ; idxValues < FSphericalCell::GetPoleSize() && !(poleDiff && localDiff); ++idxValues){
                    const FComplexe pole = octreeIterator.getCurrentCell()->getMultipole()[idxValues];
                    const FComplexe poleTyped = octreeIteratorTyped.getCurrentCell()->getMultipole()[idxValues];
                    if(!FMath::LookEqual(pole.getImag(),poleTyped.getImag()) || !FMath::LookEqual(pole.getReal(),poleTyped.getReal())){
                        poleDiff = true;
                        printf("Pole diff imag( %.15e , %.15e ) real( %.15e , %.15e)\n",
                               pole.getImag(),poleTyped.getImag(),pole.getReal(),poleTyped.getReal());
                    }
                }
                for(int idxValues = 0 ; idxValues < FSphericalCell::GetPoleSize() && !(poleDiff && localDiff); ++idxValues){
                    const FComplexe local = octreeIterator.getCurrentCell()->getLocal()[idxValues];
                    const FComplexe localTyped = octreeIteratorTyped.getCurrentCell()->getLocal()[idxValues];
                    if(!FMath::LookEqual(local.getImag(),localTyped.getImag()) || !FMath::LookEqual(local.getReal(),localTyped.getReal())){
                        localDiff = true;
                        printf("Pole diff imag( %.15e , %.15e ) real( %.15e , %.15e)\n",
                               local.getImag(),localTyped.getImag(),local.getReal(),localTyped.getReal());
                    }
                }
                if(poleDiff){
                    std::cout << "Multipole error at level " << idxLevel << "\n";
                }
                if(localDiff){
                    std::cout << "Locale error at level " << idxLevel << "\n";
                }
            } while(octreeIterator.moveRight() && octreeIteratorTyped.moveRight());

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            octreeIteratorTyped.moveUp();
            octreeIteratorTyped.gotoLeft();
        }
    }

    std::cout << "Done ..." << std::endl;

    return 0;
}



