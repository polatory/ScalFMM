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

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/Files/FPerLeafLoader.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Components/FTestParticleContainer.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Components/FTestCell.hpp"
#include "../../Src/Components/FTestKernels.hpp"

#include "../../Src/Core/FFmmAlgorithmPeriodic.hpp"

/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is impacted each other particles
  */


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FTestCell                   CellClass;
    typedef FTestParticleContainer      ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FTestKernels< CellClass, ContainerClass >         KernelClass;

    typedef FFmmAlgorithmPeriodic<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass >     FmmClass;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels          = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels     = FParameters::getValue(argc,argv,"-sh", 3);
    const long NbParticles      = FParameters::getValue(argc,argv,"-nb", 1000);
    const int PeriodicDeep      = FParameters::getValue(argc,argv,"-per", 2);
    // choose in +x dir or -/+x dir or all dirs
    int PeriodicDirs          = (FParameters::existParameter(argc,argv,"-x")?DirMinusX:0) |
                                (FParameters::existParameter(argc,argv,"+x")?DirPlusX:0) |
                                (FParameters::existParameter(argc,argv,"-y")?DirMinusY:0) |
                                (FParameters::existParameter(argc,argv,"+y")?DirPlusY:0) |
                                (FParameters::existParameter(argc,argv,"-z")?DirMinusZ:0) |
                                (FParameters::existParameter(argc,argv,"+z")?DirPlusZ:0);
    if( PeriodicDirs == 0 ) PeriodicDirs =  AllDirs;

    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    FRandomLoader loader(NbParticles);
    //FPerLeafLoader loader(NbLevels);

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    {
        FPoint particlePosition;
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particlePosition);
            tree.insert(particlePosition);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels;
    FmmClass algo( &tree, PeriodicDeep, PeriodicDirs);
    algo.setKernel(&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    { // Check that each particle has been summed with all other
        long long int counterNbPart = 0;

        tree.forEachCellLeaf([&](CellClass* cell, LeafClass* leaf){
            if(cell->getDataUp() != leaf->getSrc()->getNbParticles() ){
                    std::cout << "Problem P2M Data up = " << cell->getDataUp() <<
                                 " Size = " << leaf->getSrc()->getNbParticles() << "\n";
            }
            // we also count the number of particles.
            counterNbPart += leaf->getSrc()->getNbParticles();
        });

        if( counterNbPart != loader.getNumberOfParticles()){
            std::cout << "Problem global nb part, counter = " << counterNbPart << " created = " <<
                         loader.getNumberOfParticles() << std::endl;
        }
    }
    {
        const FTreeCoordinate repetitions = algo.repetitions();
        const int totalRepeatedBox = repetitions.getX() * repetitions.getY() * repetitions.getZ();
        std::cout << "The box is repeated " << repetitions.getX() <<" "<< repetitions.getY()<<" "<<
                     repetitions.getZ() << " there are " << totalRepeatedBox << " boxes in total\n";
        const long long NbParticlesEntireSystem = loader.getNumberOfParticles() * totalRepeatedBox;
        std::cout << "The total number of particles is "  << NbParticlesEntireSystem << "\n";

        FTreeCoordinate min, max;
        algo.repetitionsIntervals(&min, &max);
        std::cout << "Min is " << min << " Max is " << max << std::endl;

        tree.forEachLeaf([&](LeafClass* leaf){
            for(int idxPart = 0 ; idxPart < leaf->getSrc()->getNbParticles() ; ++idxPart ){
                if( NbParticlesEntireSystem - 1 != leaf->getSrc()->getDataDown()[idxPart]){
                    std::cout << "P2P probleme, should be " << NbParticlesEntireSystem - 1 <<
                                 " iter.data().getDataDown() "<< leaf->getSrc()->getDataDown()[idxPart] << std::endl;
                }
            }
        });
    }

    std::cout << "Test done..." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    return 0;
}



