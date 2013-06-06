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


#include "../../Src/Files/FFmaLoader.hpp"
#include "../../Src/Files/FTreeIO.hpp"

#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FSphericalCell                 CellClass;
    typedef FP2PParticleContainer         ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to load or retrieve an entier tree.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
    std::cout << "Opening : " << filename << "\n";

    FFmaLoader loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------
    CellClass::Init(5);
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint particlePosition;
        FReal physicalValue = 0.0;
        loader.fillParticle(&particlePosition,&physicalValue);
        tree.insert(particlePosition, physicalValue );
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Save tree ..." << std::endl;

    FTreeIO::Save<OctreeClass, CellClass, LeafClass, ContainerClass >("/tmp/tree.data", tree);

    // -----------------------------------------------------

    std::cout << "Load tree ..." << std::endl;

    FTreeIO::Load<OctreeClass, CellClass, LeafClass, ContainerClass >("/tmp/tree.data", tree);

    return 0;
}



