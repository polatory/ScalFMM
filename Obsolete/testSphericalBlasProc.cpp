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
// @FUSE_BLAS
// ================
// Keep in private GIT
// @SCALFMM_PRIVATE


#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMpi.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Utils/FMath.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/Spherical/FSphericalBlasKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "../../Src/Files/FMpiFmaGenericLoader.hpp"
#include "../../Src/Files/FMpiTreeBuilder.hpp"

#include "../../Src/Utils/FLeafBalance.hpp"

#include "../../Src/Utils/FParameterNames.hpp"


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Test Spherical HArmonic kernel with blas and using MPI.",
                         FParameterDefinitions::InputFile, FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight, FParameterDefinitions::SHDevelopment,
                         FParameterDefinitions::NbThreads);

    typedef double FReal;
    typedef FSphericalCell<FReal>         CellClass;
    typedef FP2PParticleContainer<FReal>         ContainerClass;

    typedef FSimpleLeaf<FReal, ContainerClass >                     LeafClass;
    typedef FOctree<FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalBlasKernel<FReal, CellClass, ContainerClass >     KernelClass;

    typedef FFmmAlgorithmThreadProc<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;


    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app( argc, argv);

    const int DevP = FParameters::getValue(argc,argv,FParameterDefinitions::SHDevelopment.options, 8);
    const int NbLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeSubHeight.options, 3);
    FTic counter;
    const char* const defaultFilename = (sizeof(FReal) == sizeof(float))?
                "../Data/test20k.bin.fma.single":
                "../Data/test20k.bin.fma.double";
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFilename);
    const int nbThreads = FParameters::getValue(argc,argv,FParameterDefinitions::NbThreads.options,8);
    omp_set_num_threads(nbThreads);

    std::cout << "Opening : " << filename << "\n";

    FMpiFmaGenericLoader<FReal> loader(filename, app.global());
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    CellClass::Init(DevP);


    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    if( app.global().processCount() != 1){
        //////////////////////////////////////////////////////////////////////////////////
        // Build tree from mpi loader
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Build Tree ..." << std::endl;
        counter.tic();

        struct TestParticle{
            FPoint<FReal> position;
            FReal physicalValue;
            const FPoint<FReal>& getPosition(){
                return position;
            }
        };

        TestParticle* particles = new TestParticle[loader.getNumberOfParticles()];
        memset(particles, 0, sizeof(TestParticle) * loader.getNumberOfParticles());

        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&particles[idxPart].position,&particles[idxPart].physicalValue);
        }

        FVector<TestParticle> finalParticles;
        FLeafBalance balancer;
        // FMpiTreeBuilder< FReal,TestParticle >::ArrayToTree(app.global(), particles, loader.getNumberOfParticles(),
        //                                               tree.getBoxCenter(),
        //                                               tree.getBoxWidth(),
        //                                               tree.getHeight(), &finalParticles,&balancer);
        FMpiTreeBuilder< FReal,TestParticle >::DistributeArrayToContainer(app.global(),particles,
                                                                    loader.getMyNumberOfParticles(),
                                                                    tree.getBoxCenter(),
                                                                    tree.getBoxWidth(),tree.getHeight(),
                                                                    &finalParticles, &balancer);

        for(int idx = 0 ; idx < finalParticles.getSize(); ++idx){
            tree.insert(finalParticles[idx].position,finalParticles[idx].physicalValue);

        }

        delete[] particles;

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
    }
    else{
        FPoint<FReal> position;
        FReal physicalValue;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(&position,&physicalValue);
            tree.insert(position, physicalValue);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------
    std::cout << "Create kernel..." << std::endl;

    KernelClass kernels(DevP, NbLevels,loader.getBoxWidth(), loader.getCenterOfBox());

    std::cout << "Done  " << " in " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;

    FmmClass algo(app.global(),&tree,&kernels);

    counter.tic();
    algo.execute();
    counter.tac();

    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    { // get sum forces&potential

        FReal potential = 0;
        FReal fx = 0.0, fy = 0.0, fz = 0.0;

        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();

            for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                potential += potentials[idxPart];
                fx += forcesX[idxPart];
                fy += forcesY[idxPart];
                fz += forcesZ[idxPart];
            }
        });

        std::cout << "My potential is " << potential << std::endl;

        potential = app.global().reduceSum(potential);
        fx = app.global().reduceSum(fx);
        fy = app.global().reduceSum(fy);
        fz = app.global().reduceSum(fz);


        if(app.global().processId() == 0){
            std::cout << "Foces Sum  x = " << fx << " y = " << fy << " z = " << fz << std::endl;
            std::cout << "Potential Sum = " << potential << std::endl;
        }
    }
    return 0;
}
