// ===================================================================================
// Copyright ScalFmm 2016 INRIA, Olivier Coulaud, Bérenger Bramas,
// Matthias Messner olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the
// FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
// An extension to the license is given to allow static linking of scalfmm
// inside a proprietary application (no matter its license).
// See the main license file for more details.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================

// Keep in private GIT
// @SCALFMM_PRIVATE

// @FUSE_FFT
#include <iostream>
#include <cstdio>


#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Components/FSimpleIndexedLeaf.hpp"

#include "../../Src/Utils/FPoint.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Adaptive/FAdaptChebSymKernel.hpp"
#include "Kernels/Uniform/FUnifCell.hpp"
#include "Adaptive/FAdaptUnifKernel.hpp"

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Components/FBasicKernels.hpp"

#include "../../Src/Files/FRandomLoader.hpp"

#include "../../Src/Adaptive/FAdaptiveCell.hpp"
#include "../../Src/Adaptive/FAdaptiveKernelWrapper.hpp"
#include "../../Src/Adaptive/FAbstractAdaptiveKernel.hpp"
#include "../../Src/Adaptive/FAdaptiveTestKernel.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include "../../Src/Kernels/P2P/FP2PR.hpp"


/** This program show an example of use of the fmm basic algo
  * it also check that each particles is impacted each other particles
  */

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    const FParameterNames LocalOptionMinMultipoleThreshod {
        {"-sM"},
        " s_min^M threshold for Multipole (l+1)^2 for Spherical harmonic."
    };
    const FParameterNames LocalOptionMinLocalThreshod {
        {"-SL"},
        " s_min^L threshold for Local  (l+1)^2 for Spherical harmonics."
    };

    FHelpDescribeAndExit(argc, argv,
                         "Test the adaptive FMM.",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight,LocalOptionMinMultipoleThreshod,
                         LocalOptionMinLocalThreshod);


    typedef double FReal;
    const unsigned int P = 5 ;
    typedef FChebCell<FReal,P>                                        CellClass;
    //typedef FUnifCell<FReal,P>                                        CellClass;

    typedef FP2PParticleContainerIndexed<FReal>            ContainerClass;
    typedef FSimpleIndexedLeaf<FReal,ContainerClass>    LeafClass;
    typedef FInterpMatrixKernelR<FReal>                               MatrixKernelClass;
    typedef FAdaptiveChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,P> KernelClass;
    //typedef FAdaptiveUnifKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,P> KernelClass;
    typedef FAdaptiveCell< CellClass, ContainerClass >                                        CellWrapperClass;
    typedef FAdaptiveKernelWrapper< KernelClass, CellClass, ContainerClass >   KernelWrapperClass;
    typedef FOctree< FReal, CellWrapperClass, ContainerClass , LeafClass >                  OctreeClass;
    typedef FFmmAlgorithm<OctreeClass, CellWrapperClass, ContainerClass, KernelWrapperClass, LeafClass >     FmmClass;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test the FMM algorithm.\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 7);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    const int sminM    = FParameters::getValue(argc,argv,LocalOptionMinMultipoleThreshod.options, P*P*P);
    const int sminL     = FParameters::getValue(argc,argv,LocalOptionMinLocalThreshod.options, P*P*P);
    FTic counter;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    const FReal boxWidth = 1.0;
    const FPoint<FReal> boxCenter(0.0, 0.0, 0.0);
    OctreeClass tree(NbLevels, SizeSubLevels, boxWidth, boxCenter);

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    struct Particle{
        int idx;
        FPoint<FReal> pos;
        FReal physicalValue;
        FReal forces[3];
        FReal pot;
    };

    Particle particles[3];
    memset(particles, 0, sizeof(Particle) * 3);
    {
        particles[0].idx = 0;
        particles[0].pos = FPoint<FReal>(0.1-0.5, 0.1-0.5, 0.1-0.5);
        particles[0].physicalValue = 1.0;
        particles[1].idx = 1;
        particles[1].pos = FPoint<FReal>(0.1-0.5+0.0625, 0.1-0.5+0.0625, 0.1-0.5+0.0625);
        particles[1].physicalValue = 1.0;
        particles[2].idx = 2;
        particles[2].pos = FPoint<FReal>(0.5-0.1, 0.5-0.1, 0.5-0.1);
        particles[2].physicalValue = 1.0;

        tree.insert(particles[0].pos, 0, particles[0].physicalValue);
        tree.insert(particles[1].pos, 1, particles[1].physicalValue);
        tree.insert(particles[2].pos, 2, particles[2].physicalValue);
    }

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    const MatrixKernelClass MatrixKernel;
    KernelWrapperClass kernels(NbLevels, boxWidth, boxCenter,&MatrixKernel,sminM,sminL);
    FmmClass algo(&tree,&kernels);  //FFmmAlgorithm FFmmAlgorithmThread
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    for(int idxTarget = 0 ; idxTarget < 3 ; ++idxTarget){
        for(int idxSource = idxTarget+1 ; idxSource < 3 ; ++idxSource){
            FP2PR::MutualParticles(
                particles[idxTarget].pos.getX(),particles[idxTarget].pos.getY(), particles[idxTarget].pos.getZ(),
                particles[idxTarget].physicalValue, &particles[idxTarget].forces[0], &particles[idxTarget].forces[1],
                &particles[idxTarget].forces[2], &particles[idxTarget].pot,
                particles[idxSource].pos.getX(),particles[idxSource].pos.getY(), particles[idxSource].pos.getZ(),
                particles[idxSource].physicalValue, &particles[idxSource].forces[0], &particles[idxSource].forces[1],
                &particles[idxSource].forces[2], &particles[idxSource].pot);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    tree.forEachLeaf([&](LeafClass* leaf){
        const FSize nbParticles = leaf->getSrc()->getNbParticles();
        const FReal* posx = leaf->getSrc()->getPositions()[0];
        const FReal* posy = leaf->getSrc()->getPositions()[1];
        const FReal* posz = leaf->getSrc()->getPositions()[2];
        const FReal* fx = leaf->getSrc()->getForcesX();
        const FReal* fy = leaf->getSrc()->getForcesY();
        const FReal* fz = leaf->getSrc()->getForcesZ();
        const FReal* pot = leaf->getSrc()->getPotentials();
        const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();
        for(FSize idxPart = 0 ; idxPart < nbParticles ; ++idxPart){
            std::cout << "[FMM] Particle pos " << posx[idxPart] << " " << posy[idxPart] << " " << posz[idxPart] << "\n";
            std::cout << "\t>> res pot " << pot[idxPart] << " forces " << fx[idxPart] << " " << fy[idxPart] << " " << fz[idxPart] << "\n";
            std::cout << "[Direct] Particle pos " << particles[indexes[idxPart]].pos.getX() << " " << particles[indexes[idxPart]].pos.getY() << " " << particles[indexes[idxPart]].pos.getZ() << "\n";
            std::cout << "\t>> res pot " << particles[indexes[idxPart]].pot << " forces " << particles[indexes[idxPart]].forces[0] << " " << particles[indexes[idxPart]].forces[1] << " " << particles[indexes[idxPart]].forces[2] << "\n";
            std::cout << "\n";
        }
    });

    return 0;
}




