
// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_FFT
// ================
// Keep in private GIT
// @SCALFMM_PRIVATE


#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/Core/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Uniform/FUnifKernel.hpp"

#include "../../Src/GroupTree/Uniform/FUnifCellPOD.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"
#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/GroupTree/Core/FGroupSeqAlgorithm.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskAlgorithm.hpp"
#ifdef SCALFMM_USE_OMP4
#include "../../Src/GroupTree/Core/FGroupTaskDepAlgorithm.hpp"
#endif
#ifdef SCALFMM_USE_STARPU
#include "../../Src/GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"
#include "../../Src/GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"
#endif
#include "../../Src/GroupTree/Core/FP2PGroupParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include <memory>

#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"
#include "Core/FFmmAlgorithmSectionTask.hpp"
#include "Core/FFmmAlgorithmTask.hpp"
#include "Core/FFmmAlgorithmThreadBalance.hpp"
#include "Components/FSimpleLeaf.hpp"

#include "Kernels/Uniform/FUnifCell.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Uniform/FUnifKernel.hpp"

#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FTemplate.hpp"

//#define RANDOM_PARTICLES

const FParameterNames LocalOrder { {"-order"}, "Order of the kernel"};
const FParameterNames LocalOptionOmpTask { {"-omp-task"}, "To use FFmmAlgorithmTask"};
const FParameterNames LocalOptionOmpSection { {"-omp-section"}, "To use FFmmAlgorithmSectionTask"};
const FParameterNames LocalOptionOmpBalance { {"-omp-balance"}, "To use FFmmAlgorithmThreadBalance"};

const FParameterNames LocalOptionClassic { {"-omp", "omp-classic"}, "In order to use classic parallelism"};
const FParameterNames LocalOptionBlocSize { {"-bs"}, "The size of the block of the blocked tree"};
const FParameterNames LocalOptionNoValidate { {"-no-validation"}, "To avoid comparing with direct computation"};

struct RunContainer{
    template <const int ORDER>
    static void Run(int argc, char* argv[]){
        std::cout << "Unif kernel ORDER " << ORDER << std::endl;
        // Initialize the types
        typedef double FReal;
        typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
        const MatrixKernelClass MatrixKernel;
        const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
        omp_set_num_threads(FParameters::getValue(argc,argv,FParameterDefinitions::NbThreads.options, omp_get_max_threads()));

        if(FParameters::existParameter(argc, argv, LocalOptionClassic.options)){
            std::cout << "\n>> Using " << omp_get_max_threads() << " omp threads.\n" << std::endl;
            const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);

            // init particles position and physical value
            struct TestParticle{
                FPoint<FReal> position;
                FReal forces[3];
                FReal physicalValue;
                FReal potential;
            };

            // open particle file
    #ifdef RANDOM_PARTICLES
            FRandomLoader<FReal> loader(FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 2000), 1.0, FPoint<FReal>(0,0,0), 0);
    #else
            const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
            FFmaGenericLoader<FReal> loader(filename);
    #endif
            FAssertLF(loader.isOpen());

            TestParticle* const particles = new TestParticle[loader.getNumberOfParticles()];
            for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                FPoint<FReal> position;
                FReal physicalValue = 0.0;
    #ifdef RANDOM_PARTICLES
                physicalValue = 0.10;
                loader.fillParticle(&position);
    #else
                loader.fillParticle(&position, &physicalValue);
    #endif
                // get copy
                particles[idxPart].position       = position;
                particles[idxPart].physicalValue  = physicalValue;
                particles[idxPart].potential      = 0.0;
                particles[idxPart].forces[0]      = 0.0;
                particles[idxPart].forces[1]      = 0.0;
                particles[idxPart].forces[2]      = 0.0;
            }

            ////////////////////////////////////////////////////////////////////

            FTic time;

            if(FParameters::existParameter(argc, argv, LocalOptionNoValidate.options) == false){
                // begin direct computation
                std::cout << "\nDirect computation ... " << std::endl;
                time.tic();
                {
                    for(FSize idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
                        for(FSize idxOther =  idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                            FP2P::MutualParticles(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                                  particles[idxTarget].position.getZ(), particles[idxTarget].physicalValue,
                                                  &particles[idxTarget].forces[0], &particles[idxTarget].forces[1],
                                    &particles[idxTarget].forces[2], &particles[idxTarget].potential,
                                    particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                    particles[idxOther].position.getZ(), particles[idxOther].physicalValue,
                                    &particles[idxOther].forces[0], &particles[idxOther].forces[1],
                                    &particles[idxOther].forces[2], &particles[idxOther].potential,
                                    &MatrixKernel);
                        }
                    }
                }
                time.tac();
                std::cout << "Done  " << "(@Direct computation = "
                          << time.elapsed() << "s)." << std::endl;

            } // end direct computation

            // typedefs
            typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
            typedef FSimpleLeaf<FReal, ContainerClass >  LeafClass;
            typedef FUnifCell<FReal,ORDER> CellClass;
            typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
            typedef FUnifKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;

            // init oct-tree
            OctreeClass tree(NbLevels, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());


            { // -----------------------------------------------------
                std::cout << "Creating & Inserting " << loader.getNumberOfParticles()
                          << " particles ..." << std::endl;
                std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SubTreeHeight << std::endl;
                time.tic();

                for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                    // put in tree
                    tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue);
                }

                time.tac();
                std::cout << "Done  " << "(@Creating and Inserting Particles = "
                          << time.elapsed() << "s)." << std::endl;
            } // -----------------------------------------------------

            { // -----------------------------------------------------
                std::cout << "\nLagrange/Uniform grid FMM (ORDER="<< ORDER << ") ... " << std::endl;
                KernelClass kernels(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(),&MatrixKernel);
                if(FParameters::existParameter(argc, argv, LocalOptionOmpBalance.options)){
                    typedef FFmmAlgorithmThreadBalance<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
                    std::cout << "Using FFmmAlgorithmThreadBalance " << std::endl;
                    FmmClass algorithm(&tree, &kernels);
                    time.tic();
                    algorithm.execute();
                    time.tac();
                    std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
                }
                else if(FParameters::existParameter(argc, argv, LocalOptionOmpTask.options)){
                    typedef FFmmAlgorithmTask<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
                    std::cout << "Using FFmmAlgorithmTask " << std::endl;
                    FmmClass algorithm(&tree, &kernels);
                    time.tic();
                    algorithm.execute();
                    time.tac();
                    std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
                }
                else if(FParameters::existParameter(argc, argv, LocalOptionOmpSection.options)){
                    typedef FFmmAlgorithmSectionTask<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
                    std::cout << "Using FFmmAlgorithmSectionTask " << std::endl;
                    FmmClass algorithm(&tree, &kernels);
                    time.tic();
                    algorithm.execute();
                    time.tac();
                    std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
                }
                else {
                    typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
                    std::cout << "Using FFmmAlgorithmThread " << std::endl;
                    FmmClass algorithm(&tree, &kernels);
                    time.tic();
                    algorithm.execute();
                    time.tac();
                    std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << "s)." << std::endl;
                }
            } // -----------------------------------------------------


            if(FParameters::existParameter(argc, argv, LocalOptionNoValidate.options) == false){
                // -----------------------------------------------------
                FMath::FAccurater<FReal> potentialDiff;
                FMath::FAccurater<FReal> fx, fy, fz;

                FReal checkPotential[20000];

                { // Check that each particle has been summed with all other

                    tree.forEachLeaf([&](LeafClass* leaf){
                        const FReal*const potentials = leaf->getTargets()->getPotentials();
                        const FReal*const forcesX = leaf->getTargets()->getForcesX();
                        const FReal*const forcesY = leaf->getTargets()->getForcesY();
                        const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
                        const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                        const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

                        for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                            const FSize indexPartOrig = indexes[idxPart];
                            //PB: store potential in nbParticles array
                            checkPotential[indexPartOrig]=potentials[idxPart];

                            potentialDiff.add(particles[indexPartOrig].potential,potentials[idxPart]);
                            fx.add(particles[indexPartOrig].forces[0],forcesX[idxPart]);
                            fy.add(particles[indexPartOrig].forces[1],forcesY[idxPart]);
                            fz.add(particles[indexPartOrig].forces[2],forcesZ[idxPart]);
                        }
                    });
                }

                // Print for information
                std::cout << "Potential " << potentialDiff << std::endl;
                std::cout << "Fx " << fx << std::endl;
                std::cout << "Fy " << fy << std::endl;
                std::cout << "Fz " << fz << std::endl;
            } // -----------------------------------------------------
        }
        else{
            typedef FUnifCellPODCore         GroupCellSymbClass;
            typedef FUnifCellPODPole<FReal,ORDER>  GroupCellUpClass;
            typedef FUnifCellPODLocal<FReal,ORDER> GroupCellDownClass;
            typedef FUnifCellPOD<FReal,ORDER>      GroupCellClass;

            std::cout << "\n>> Using " << omp_get_max_threads() << " omp threads.\n" << std::endl;

            typedef FP2PGroupParticleContainer<FReal>          GroupContainerClass;
            typedef FGroupTree< FReal, GroupCellClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass, 1, 4, FReal>  GroupOctreeClass;
    #ifdef SCALFMM_USE_STARPU
            typedef FStarPUAllCpuCapacities<FUnifKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER>> GroupKernelClass;
            typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
            typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper > GroupAlgorithm;
            std::cout << "Using FGroupTaskStarPUAlgorithm" << std::endl;
    #elif defined(SCALFMM_USE_OMP4)
            typedef FUnifKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER> GroupKernelClass;
            // Set the number of threads
            omp_set_num_threads(FParameters::getValue(argc,argv,FParameterDefinitions::NbThreads.options, omp_get_max_threads()));
            typedef FGroupTaskDepAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass,
                    GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
            std::cout << "Using FGroupTaskDepAlgorithm" << std::endl;
    #else
            typedef FUnifKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER> GroupKernelClass;
            //typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
            typedef FGroupTaskAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
            std::cout << "Using FGroupTaskAlgorithm" << std::endl;
    #endif
            // Get params
            const int groupSize     = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

            // Load the particles
    #ifdef RANDOM_PARTICLES
            FRandomLoader<FReal> loader(FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 2000), 1.0, FPoint<FReal>(0,0,0), 0);
    #else
            const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
            FFmaGenericLoader<FReal> loader(filename);
    #endif
            FAssertLF(loader.isOpen());
            FTic timer;

            FP2PParticleContainer<FReal> allParticles;
            for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                FPoint<FReal> particlePosition;
                FReal physicalValue;
    #ifdef RANDOM_PARTICLES
                physicalValue = 0.10;
                loader.fillParticle(&particlePosition);
    #else
                loader.fillParticle(&particlePosition, &physicalValue);
    #endif
                allParticles.push(particlePosition, physicalValue);
            }

            // Put the data into the tree
            timer.tic();
            GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
            groupedTree.printInfoBlocks();
            std::cout << "Tree created in " << timer.tacAndElapsed() << "s\n";

            // Run the algorithm
            GroupKernelClass groupkernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);
            GroupAlgorithm groupalgo(&groupedTree,&groupkernel);

            timer.tic();
            groupalgo.execute();
            std::cout << "Done  " << "(@Algorithm = " << timer.tacAndElapsed() << "s)." << std::endl;

            // Validate the result
            if(FParameters::existParameter(argc, argv, LocalOptionNoValidate.options) == false){
                FSize offsetParticles = 0;
                FReal*const allPhysicalValues = allParticles.getPhysicalValues();
                FReal*const allPosX = const_cast<FReal*>( allParticles.getPositions()[0]);
                FReal*const allPosY = const_cast<FReal*>( allParticles.getPositions()[1]);
                FReal*const allPosZ = const_cast<FReal*>( allParticles.getPositions()[2]);

                groupedTree.template forEachCellLeaf<FP2PGroupParticleContainer<FReal> >([&](GroupCellClass cellTarget, FP2PGroupParticleContainer<FReal> * leafTarget){
                    const FReal*const physicalValues = leafTarget->getPhysicalValues();
                    const FReal*const posX = leafTarget->getPositions()[0];
                    const FReal*const posY = leafTarget->getPositions()[1];
                    const FReal*const posZ = leafTarget->getPositions()[2];
                    const FSize nbPartsInLeafTarget = leafTarget->getNbParticles();

                    for(FSize idxPart = 0 ; idxPart < nbPartsInLeafTarget ; ++idxPart){
                        allPhysicalValues[offsetParticles + idxPart] = physicalValues[idxPart];
                        allPosX[offsetParticles + idxPart] = posX[idxPart];
                        allPosY[offsetParticles + idxPart] = posY[idxPart];
                        allPosZ[offsetParticles + idxPart] = posZ[idxPart];
                    }

                    offsetParticles += nbPartsInLeafTarget;
                });

                FAssertLF(offsetParticles == loader.getNumberOfParticles());

                FReal*const allDirectPotentials = allParticles.getPotentials();
                FReal*const allDirectforcesX = allParticles.getForcesX();
                FReal*const allDirectforcesY = allParticles.getForcesY();
                FReal*const allDirectforcesZ = allParticles.getForcesZ();

                for(int idxTgt = 0 ; idxTgt < offsetParticles ; ++idxTgt){
                    for(int idxMutual = idxTgt + 1 ; idxMutual < offsetParticles ; ++idxMutual){
                        FP2PR::MutualParticles(
                                    allPosX[idxTgt],allPosY[idxTgt],allPosZ[idxTgt], allPhysicalValues[idxTgt],
                                    &allDirectforcesX[idxTgt], &allDirectforcesY[idxTgt], &allDirectforcesZ[idxTgt], &allDirectPotentials[idxTgt],
                                    allPosX[idxMutual],allPosY[idxMutual],allPosZ[idxMutual], allPhysicalValues[idxMutual],
                                    &allDirectforcesX[idxMutual], &allDirectforcesY[idxMutual], &allDirectforcesZ[idxMutual], &allDirectPotentials[idxMutual]
                                    );
                    }
                }

                FMath::FAccurater<FReal> potentialDiff;
                FMath::FAccurater<FReal> fx, fy, fz;
                offsetParticles = 0;
                groupedTree.template forEachCellLeaf<FP2PGroupParticleContainer<FReal> >([&](GroupCellClass cellTarget, FP2PGroupParticleContainer<FReal> * leafTarget){
                    const FReal*const potentials = leafTarget->getPotentials();
                    const FReal*const forcesX = leafTarget->getForcesX();
                    const FReal*const forcesY = leafTarget->getForcesY();
                    const FReal*const forcesZ = leafTarget->getForcesZ();
                    const FSize nbPartsInLeafTarget = leafTarget->getNbParticles();

                    for(int idxTgt = 0 ; idxTgt < nbPartsInLeafTarget ; ++idxTgt){
                        potentialDiff.add(allDirectPotentials[idxTgt + offsetParticles], potentials[idxTgt]);
                        fx.add(allDirectforcesX[idxTgt + offsetParticles], forcesX[idxTgt]);
                        fy.add(allDirectforcesY[idxTgt + offsetParticles], forcesY[idxTgt]);
                        fz.add(allDirectforcesZ[idxTgt + offsetParticles], forcesZ[idxTgt]);
                    }

                    offsetParticles += nbPartsInLeafTarget;
                });

                std::cout << "Error : Potential " << potentialDiff << "\n";
                std::cout << "Error : fx " << fx << "\n";
                std::cout << "Error : fy " << fy << "\n";
                std::cout << "Error : fz " << fz << "\n";
            }
        }

    }
};


int main(int argc, char* argv[]){
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight, FParameterDefinitions::OctreeSubHeight,
                     #ifdef RANDOM_PARTICLES
                         FParameterDefinitions::NbParticles,
                     #else
                         FParameterDefinitions::InputFile,
                     #endif
                         FParameterDefinitions::NbThreads,
                         LocalOptionBlocSize, LocalOptionNoValidate, LocalOptionClassic,
                         LocalOptionOmpTask, LocalOptionOmpSection, LocalOptionOmpBalance,
                         LocalOrder);

    const int order = FParameters::getValue(argc,argv,LocalOrder.options, 5);
    FRunIf::Run<int, 3, 8, 1, RunContainer>(order, argc, argv);

    return 0;
}
