// ===================================================================================
// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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
#include <iomanip>

#include <cstdio>  //printf
#include <cstdlib>
#include <cstring>  //memset

#include <cmath>
#include <algorithm>
#include <string>

#include  "ScalFmmConfig.h"
#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FParameterNames.hpp"
#include "../../Src/Utils/FIOVtk.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/Kernels/P2P/FP2PR.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         ">> This executable test the efficiency of the computation of the P2P",
                         FParameterDefinitions::NbParticles);

    const int nbParticles = FParameters::getValue(argc, argv, FParameterDefinitions::NbParticles.options, 1000);
    std::cout << "Test with " << nbParticles << " particles." << std::endl;

    //////////////////////////////////////////////////////////

    FRandomLoader loader(nbParticles*2);

    FTic timer;
    FP2PParticleContainer<> leaf1;
    for(int idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
        FPoint pos;
        loader.fillParticle(&pos);
        leaf1.push(pos, 1.0);
    }

    FP2PParticleContainer<> leaf2;
    for(int idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
        FPoint pos;
        loader.fillParticle(&pos);
        leaf2.push(pos, 1.0);
    }
    FP2PParticleContainer<> * const pleaf2 = &leaf2;
    std::cout << "Timer taken to create and insert the particles = " << timer.tacAndElapsed() << "s" << std::endl;

    //////////////////////////////////////////////////////////

    timer.tic();
    FP2PR::FullMutual<FP2PParticleContainer<>>( &leaf1, &pleaf2, 1);
    timer.tac();
    std::cout << "Timer taken by FullMutual = " << timer.elapsed() << "s" << std::endl;

    timer.tic();
    FP2PR::FullRemote<FP2PParticleContainer<>>( &leaf1, &pleaf2, 1);
    timer.tac();
    std::cout << "Timer taken by FullRemote = " << timer.elapsed() << "s" << std::endl;

    //////////////////////////////////////////////////////////

    return 0;
}
