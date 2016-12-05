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

// ==== CMAKE =====
// @FUSE_MPI
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>


#include "../../Src/Kernels/Rotation/FRotationCell.hpp"
#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMemUtils.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Files/FRandomLoader.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"

#include "../../Src/Utils/FLeafBalance.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

/**
 * This program runs the FMM Algorithm Distributed with the Rotation kernel
 */

// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Test with MPI the chebyshev FMM and compare it to the direct computation for debugging purpose.",
                         FParameterDefinitions::NbParticles, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::NbThreads);

    typedef double FReal;


    FMpi app(argc,argv);

    const FSize nbParticles       = FParameters::getValue(argc,argv, FParameterDefinitions::NbParticles.options, 10000000ULL);
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    FTic time;

    std::cout << ">> This executable has to be used to test Proc Rotation Algorithm. \n";

    // init particles position and physical value
    struct TestParticle{
        FSize idxPart;
        FPoint<FReal> position;
        FReal physicalValue;
        const FPoint<FReal>& getPosition(){
            return position;
        }
    };

    // open particle file
    std::cout << "Creating : " << nbParticles << "\n" << std::endl;

    time.tic();
    const FSize totalNbParticles = nbParticles*app.global().processCount();
    TestParticle* particles = new TestParticle[totalNbParticles];
    memset(particles,0,(unsigned int) (sizeof(TestParticle)*totalNbParticles));
    for(int idxProc = 0 ; idxProc < app.global().processCount() ; ++idxProc){
        FRandomLoader<FReal> loader(nbParticles, 1.0, FPoint<FReal>(0,0,0), idxProc);
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particles[idxPart + idxProc*nbParticles].position);
            particles[idxPart + idxProc*nbParticles].physicalValue = 1.0;
            particles[idxPart + idxProc*nbParticles].idxPart = idxPart + idxProc*nbParticles;
        }
    }

    FVector<TestParticle> finalParticles;
    FLeafBalance balancer;
    FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(app.global(),&particles[app.global().processId()*nbParticles],
                                                                nbParticles,
                                                                FPoint<FReal>(0,0,0),
                                                                1.0,TreeHeight,
                                                                &finalParticles, &balancer);

    app.global().barrier();

    std::cout << "Testing : " << finalParticles.getSize()  << "\n" << std::endl;

    for(FSize idxRes = 0 ; idxRes < finalParticles.getSize() ; ++idxRes){
        FAssertLF(0 <= finalParticles[idxRes].idxPart, "idxRes ", idxRes, " finalParticles[idxRes].idxPart ", finalParticles[idxRes].idxPart);
        FAssertLF(finalParticles[idxRes].idxPart < totalNbParticles, "idxRes ", idxRes, " finalParticles[idxRes].idxPart ", finalParticles[idxRes].idxPart);

        const TestParticle correctPart = particles[finalParticles[idxRes].idxPart];
        const TestParticle testPart = finalParticles[idxRes];

        FAssertLF(testPart.idxPart == correctPart.idxPart);
        FAssertLF(testPart.position.getX() == correctPart.position.getX());
        FAssertLF(testPart.position.getY() == correctPart.position.getY());
        FAssertLF(testPart.position.getZ() == correctPart.position.getZ());
        FAssertLF(testPart.physicalValue == correctPart.physicalValue);
    }

    std::cout << "Done\n" << std::endl;

    app.global().barrier();

    std::unique_ptr<int[]> particlesExist(new int[totalNbParticles]);
    memset(particlesExist.get(), 0, sizeof(int)*totalNbParticles);

    for(FSize idxRes = 0 ; idxRes < finalParticles.getSize() ; ++idxRes){
        FAssertLF(particlesExist[finalParticles[idxRes].idxPart] == 0);
        particlesExist[finalParticles[idxRes].idxPart] = 1;
    }

    std::unique_ptr<int[]> particlesReduced(new int[totalNbParticles]);
    memset(particlesReduced.get(), 0, sizeof(int)*totalNbParticles);

    FAssert(totalNbParticles <= std::numeric_limits<int>::max());
    FMpi::Assert(MPI_Allreduce(particlesExist.get(), particlesReduced.get(), int(totalNbParticles),
                            MPI_INT, MPI_SUM,
                            app.global().getComm()), __LINE__);

    for(FSize idxPart = 0 ; idxPart < totalNbParticles ; ++idxPart){
        FAssertLF(particlesReduced[idxPart] == 1, idxPart, " " , particlesReduced[idxPart]);
    }

    return 0;
}
