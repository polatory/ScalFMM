// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_FFT
// ================

// Keep in private GIT
//

#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <string>

#include "ScalFmmConfig.h"
//Generals include
#include "Files/FFmaGenericLoader.hpp"
#include "Containers/FVector.hpp"

//For interpolation kernel
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"

//For Cheb Sym kernel
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
//For Unif Kernel
#include "Kernels/Uniform/FUnifCell.hpp"
#include "Kernels/Uniform/FUnifKernel.hpp"
//For Rotation Kernel
#include "Kernels/Rotation/FRotationKernel.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"
// spherical kernel
#include "Kernels/Spherical/FSphericalCell.hpp"
#include "Kernels/Spherical/FSphericalBlasKernel.hpp"
#include "Kernels/Spherical/FSphericalBlockBlasKernel.hpp"
//Taylor Kernel
#include "../../Src/Kernels/Taylor/FTaylorCell.hpp"
#include "../../Src/Kernels/Taylor/FTaylorKernel.hpp"


#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Utils/FParameters.hpp"

#include "Containers/FOctree.hpp"

#ifdef _OPENMP
#include "Core/FFmmAlgorithmThread.hpp"
#else
#include "Core/FFmmAlgorithm.hpp"
#endif

#include "Utils/FTemplate.hpp"

typedef double FReal;

void usage() {
    std::cout << "Driver for all kernel (1/r kernel)" << std::endl;
    std::cout <<     "Options  "<< std::endl
              <<     "      -help         to see the parameters    " << std::endl
              <<        "      -depth       the depth of the octree   "<< std::endl
              <<        "      -subdepth  specifies the size of the sub octree   " << std::endl
              <<     "      -f   name    name specifies the name of the particle distribution" << std::endl
              <<     "      -t  n  specifies the number of threads used in the computations" << std::endl;
}

template <class FReal, class OctreeClass, class LeafClass>
void checkResAndPrint(OctreeClass * tree, FmaRWParticle<FReal, 8,8> * const particles,FReal energyD, std::string name,int order){
    FReal energy = 0.0;
    FMath::FAccurater<FReal> potentialDiff;
    FMath::FAccurater<FReal> fx, fy, fz;

    { // Check that each particle has been summed with all other

        tree->forEachLeaf([&](LeafClass* leaf){

                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                const FReal*const potentials        = leaf->getTargets()->getPotentials();
                const FReal*const forcesX            = leaf->getTargets()->getForcesX();
                const FReal*const forcesY            = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ            = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf           = leaf->getTargets()->getNbParticles();
                const FVector<FSize>& indexes       = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                    potentialDiff.add(particles[indexPartOrig].getPotential(),potentials[idxPart]);
                    fx.add(particles[indexPartOrig].getForces()[0],forcesX[idxPart]);
                    fy.add(particles[indexPartOrig].getForces()[1],forcesY[idxPart]);
                    fz.add(particles[indexPartOrig].getForces()[2],forcesZ[idxPart]);
                    energy += potentials[idxPart]*physicalValues[idxPart] ;
                }
            });
    }
    // Print for information

    std::cout << "FChebSymKernel Energy "  << FMath::Abs(energy-energyD) /energyD << std::endl;
    std::cout << "FChebSymKernel Potential " << potentialDiff << std::endl;
    std::cout << "FChebSymKernel Fx " << fx << std::endl;
    std::cout << "FChebSymKernel Fy " << fy << std::endl;
    std::cout << "FChebSymKernel Fz " << fz << std::endl;
    //Write res to file
    std::ofstream outfile;

    outfile.open(name, std::ios_base::app);
    if(outfile.tellp() == 0){
        outfile << name+"_Order \t" << name+"_Energy \t"
                << name+"_potential_L2 \t" << name+"_potential_RMS \t" << name+"_potential_Inf \t"
                << name+"_fx_L2 \t" << name+"_fx_RMS \t" << name+"_fx_Inf \t"
                << name+"_fy_L2 \t" << name+"_fy_RMS \t" << name+"_fy_Inf \t"
                << name+"_fz_L2 \t" << name+"_fz_RMS \t" << name+"_fz_Inf \n";
    }
    outfile << order << "\t" << FMath::Abs(energy-energyD) /energyD << "\t"
            << potentialDiff.getL2Norm() << "\t" << potentialDiff.getRMSError() << "\t" << potentialDiff.getInfNorm() << "\t"
            << fx.getL2Norm() << "\t" << fx.getRMSError() << "\t" << fx.getInfNorm() << "\t"
            << fy.getL2Norm() << "\t" << fy.getRMSError() << "\t" << fy.getInfNorm() << "\t"
            << fz.getL2Norm() << "\t" << fz.getRMSError() << "\t" << fz.getInfNorm() << "\n";
}


// Simply create particles and try the CHEBYSHEV kernels
struct ChebMainStruct{
    template <const unsigned int ORDER>
    static void For(int argc, char* argv[],FSize nbParticles,int TreeHeight, int SubTreeHeight,
                    FReal BoxWidth, FPoint<FReal>& CenterOfBox,FmaRWParticle<FReal, 8,8> * const particles,
                    FReal energyD,FReal totPhysicalValue)
    {
        //Global timer to be used
        FTic time;
        std::cout << "Kernels Chebyshev\n"
                  << "Current Order " << ORDER << "\n" << std::endl;

        //Start of kernels there
        {//Chebyshev
            typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
            typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
            typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
            typedef FChebCell<FReal,ORDER> CellClass;
            typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;

            typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
            typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
            // init oct-tree
            OctreeClass tree(TreeHeight, SubTreeHeight,  BoxWidth,CenterOfBox);
            // Create Matrix Kernel
            const MatrixKernelClass MatrixKernel;
            { // -----------------------------------------------------
                for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                    // put in tree
                    tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                }
            } // -----------------------------------------------------
            { // -----------------------------------------------------
                time.tic();
                KernelClass kernels(TreeHeight, BoxWidth, CenterOfBox,&MatrixKernel);
                FmmClass algorithm(&tree, &kernels);
                algorithm.execute();
                time.tac();
                std::cout <<"(FChebSymKernel @Algorithm = " << time.elapsed() << " s)." << std::endl;
            } // -----------------------------------------------------
            { // Check that each particle has been summed with all other
                checkResAndPrint<FReal,OctreeClass,LeafClass>(&tree,particles,energyD,std::string("Chebyshev.res"),ORDER);
            }

        }
    }

};

// Simply create particles and try the CHEBYSHEV kernels
struct UnifMainStruct{
    template <const unsigned int ORDER>
    static void For(int argc, char* argv[],FSize nbParticles,int TreeHeight, int SubTreeHeight,
                    FReal BoxWidth, FPoint<FReal>& CenterOfBox,FmaRWParticle<FReal, 8,8> * const particles,
                    FReal energyD,FReal totPhysicalValue)
    {
        //Global timer to be used
        FTic time;
        std::cout << "Kernels Lagrange\n"
                  << "Current Order " << ORDER << "\n" << std::endl;

        //Start of kernels there
        {//Lagrange

            typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
            typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
            typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
            typedef FUnifCell<FReal,ORDER> CellClass;
            typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;

            typedef FUnifKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
            typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
            // init oct-tree
            OctreeClass tree(TreeHeight, SubTreeHeight,  BoxWidth,CenterOfBox);
            // Create Matrix Kernel
            const MatrixKernelClass MatrixKernel;
            { // -----------------------------------------------------
                for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                    // put in tree
                    tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                }
            } // -----------------------------------------------------
            { // -----------------------------------------------------
                time.tic();
                KernelClass kernels(TreeHeight, BoxWidth, CenterOfBox,&MatrixKernel);
                FmmClass algorithm(&tree, &kernels);
                algorithm.execute();
                time.tac();
                std::cout <<"(FUnifKernel @Algorithm = " << time.elapsed() << " s)." << std::endl;
            } // -----------------------------------------------------
            { // Check that each particle has been summed with all other
                checkResAndPrint<FReal,OctreeClass,LeafClass>(&tree,particles,energyD,std::string("Lagrange.res"),ORDER);
            }

        }
    }

};

struct RotMainStruct{
    template <const unsigned int ORDER>
    static void For(int argc, char* argv[],FSize nbParticles,int TreeHeight, int SubTreeHeight,
                    FReal BoxWidth, FPoint<FReal>& CenterOfBox,FmaRWParticle<FReal, 8,8> * const particles,
                    FReal energyD,FReal totPhysicalValue)
    {
        //Global timer to be used
        FTic time;
        std::cout << "Kernels Rotation\n"
                  << "Current Order " << ORDER << "\n" << std::endl;

        //Start of kernels there
        {//Rotation

            typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
            typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
            typedef FRotationCell<FReal,ORDER> CellClass;
            typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;

            typedef FRotationKernel<FReal, CellClass,ContainerClass,ORDER> KernelClass;
            typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
            // init oct-tree
            OctreeClass tree(TreeHeight, SubTreeHeight,  BoxWidth,CenterOfBox);
            { // -----------------------------------------------------
                for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                    // put in tree
                    tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                }
            } // -----------------------------------------------------
            { // -----------------------------------------------------
                time.tic();
                KernelClass kernels(TreeHeight, BoxWidth, CenterOfBox);
                FmmClass algorithm(&tree, &kernels);
                algorithm.execute();
                time.tac();
                std::cout <<"(FRotationKernel @Algorithm = " << time.elapsed() << " s)." << std::endl;
            } // -----------------------------------------------------
            { // Check that each particle has been summed with all other
                checkResAndPrint<FReal,OctreeClass,LeafClass>(&tree,particles,energyD,std::string("Rotation.res"),ORDER);
            }

        }
    }

};

struct TaylorMainStruct{
    template <const unsigned int ORDER>
    static void For(int argc, char* argv[],FSize nbParticles,int TreeHeight, int SubTreeHeight,
                    FReal BoxWidth, FPoint<FReal>& CenterOfBox,FmaRWParticle<FReal, 8,8> * const particles,
                    FReal energyD,FReal totPhysicalValue)
    {
        //Global timer to be used
        FTic time;
        std::cout << "Kernels Taylor\n"
                  << "Current Order " << ORDER << "\n" << std::endl;

        //Start of kernels there
        {//Taylor

            typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
            typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
            typedef FTaylorCell<FReal,ORDER,1> CellClass;
            typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;

            typedef FTaylorKernel<FReal,CellClass,ContainerClass,ORDER,1> KernelClass;
            typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
            // init oct-tree
            OctreeClass tree(TreeHeight, SubTreeHeight,  BoxWidth,CenterOfBox);
            { // -----------------------------------------------------
                for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                    // put in tree
                    tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                }
            } // -----------------------------------------------------
            { // -----------------------------------------------------
                time.tic();
                KernelClass kernels(TreeHeight, BoxWidth, CenterOfBox);
                FmmClass algorithm(&tree, &kernels);
                algorithm.execute();
                time.tac();
                std::cout <<"(FTaylorKernel @Algorithm = " << time.elapsed() << " s)." << std::endl;
            } // -----------------------------------------------------
            { // Check that each particle has been summed with all other
                checkResAndPrint<FReal,OctreeClass,LeafClass>(&tree,particles,energyD,std::string("Taylor.res"),ORDER);
            }

        }
    }

};


// Simply create particles and try the SPHERICAL kernels
struct SphericalBlasMainStruct{
    template <const unsigned int ORDER>
    static void For(int argc, char* argv[],FSize nbParticles,int TreeHeight, int SubTreeHeight,
                    FReal BoxWidth, FPoint<FReal>& CenterOfBox,FmaRWParticle<FReal, 8,8> * const particles,
                    FReal energyD,FReal totPhysicalValue)
    {
        //Global timer to be used
        FTic time;
        std::cout << "Kernels Spherical\n"
                  << "Current Order " << ORDER << "\n" << std::endl;

        //Start of kernels there
        {//Spherical

            typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
            typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
            typedef FSphericalCell<FReal> CellClass;
            typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
            CellClass::Init(ORDER);
            typedef FSphericalBlasKernel<FReal,CellClass,ContainerClass> KernelClass;
            typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
            // init oct-tree
            OctreeClass tree(TreeHeight, SubTreeHeight,  BoxWidth,CenterOfBox);
            { // -----------------------------------------------------
                for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                    // put in tree
                    tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                }
            } // -----------------------------------------------------
            { // -----------------------------------------------------
                time.tic();
                KernelClass kernels(ORDER, TreeHeight, BoxWidth, CenterOfBox);
                FmmClass algorithm(&tree, &kernels);
                algorithm.execute();
                time.tac();
                std::cout <<"(FSphericalBlasKernel @Algorithm = " << time.elapsed() << " s)." << std::endl;
            } // -----------------------------------------------------
            { // Check that each particle has been summed with all other
                checkResAndPrint<FReal,OctreeClass,LeafClass>(&tree,particles,energyD,std::string("Blas.res"),ORDER);
            }

        }
    }

};

// Simply create particles and try the SPHERICAL kernels
struct SphericalBlockBlasMainStruct{
    template <const unsigned int ORDER>
    static void For(int argc, char* argv[],FSize nbParticles,int TreeHeight, int SubTreeHeight,
                    FReal BoxWidth, FPoint<FReal>& CenterOfBox,FmaRWParticle<FReal, 8,8> * const particles,
                    FReal energyD,FReal totPhysicalValue)
    {
        //Global timer to be used
        FTic time;
        std::cout << "Kernels Spherical Block Blas\n"
                  << "Current Order " << ORDER << "\n" << std::endl;

        //Start of kernels there
        {//Spherical

            typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
            typedef FSimpleLeaf<FReal, ContainerClass> LeafClass;
            typedef FSphericalCell<FReal> CellClass;
            typedef FOctree<FReal, CellClass,ContainerClass,LeafClass> OctreeClass;
            CellClass::Init(ORDER);
            typedef FSphericalBlockBlasKernel<FReal,CellClass,ContainerClass> KernelClass;
            typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
            // init oct-tree
            OctreeClass tree(TreeHeight, SubTreeHeight,  BoxWidth,CenterOfBox);
            { // -----------------------------------------------------
                for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                    // put in tree
                    tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                }
            } // -----------------------------------------------------
            { // -----------------------------------------------------
                time.tic();
                KernelClass kernels(ORDER, TreeHeight, BoxWidth, CenterOfBox);
                FmmClass algorithm(&tree, &kernels);
                algorithm.execute();
                time.tac();
                std::cout <<"(FSphericalBlockBlasKernel @Algorithm = " << time.elapsed() << " s)." << std::endl;
            } // -----------------------------------------------------
            { // Check that each particle has been summed with all other
                checkResAndPrint<FReal,OctreeClass,LeafClass>(&tree,particles,energyD,std::string("BlockBlas.res"),ORDER);
            }

        }
    }

};

int main(int argc, char** argv){
    if(FParameters::existParameter(argc, argv, "-h")||FParameters::existParameter(argc, argv, "-help")){
        usage() ;
        std::cout << "Driver for testing different approximations  for the  1/r kernel" << std::endl;
        exit(EXIT_SUCCESS);
    }

    // get info from commande line
    const std::string  filename(FParameters::getStr(argc,argv,"-f", "../Data/UTest/unitCubeRef20kDouble.bfma"));
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-depth", 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-subdepth", 2);
    const unsigned int NbThreads      = FParameters::getValue(argc, argv, "-t", omp_get_max_threads());



    //Open files
    FFmaGenericLoader<FReal> loader(filename);
    if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");
    FSize nbParticles = loader.getNumberOfParticles() ;
    FmaRWParticle<FReal, 8,8>* const particles = new FmaRWParticle<FReal, 8,8>[nbParticles];
    //
    loader.fillParticle(particles,nbParticles);

    //Write a header
    std::cout << "Number of particule "<< nbParticles << "\n"
              << "TreeHeight " << TreeHeight << "\n"
              << "NbThreads " << NbThreads << "\n"
              <<std::endl;

    //Direct computation
    //  Compute direct energy
    FReal energyD =0.0, totPhysicalValue =0.0;

#pragma omp parallel for reduction(+:energyD,totPhysicalValue)
    for(FSize idx = 0 ; idx <  loader.getNumberOfParticles()  ; ++idx){
        energyD             +=  particles[idx].getPotential()*particles[idx].getPhysicalValue() ;
        totPhysicalValue += particles[idx].getPhysicalValue() ;
    }
    const int startInter = 2;
    const int endInter = 13;
    const int start = 2;
    const int end = 30;
    const int step = 1;
    const int step2 = 2;

    FForAll::For<unsigned int, startInter, endInter, step, ChebMainStruct>( argc, argv,
                                                         nbParticles,TreeHeight, SubTreeHeight,
                                                         loader.getBoxWidth(), loader.getCenterOfBox(),particles,
                                                         energyD,totPhysicalValue);
    FForAll::For<unsigned int, startInter, endInter, step, UnifMainStruct>( argc, argv,
                                                         nbParticles,TreeHeight, SubTreeHeight,
                                                         loader.getBoxWidth(), loader.getCenterOfBox(),particles,
                                                         energyD,totPhysicalValue);
    FForAll::For<unsigned int, start, end, step2, RotMainStruct>( argc, argv,
                                                         nbParticles,TreeHeight, SubTreeHeight,
                                                         loader.getBoxWidth(), loader.getCenterOfBox(),particles,
                                                         energyD,totPhysicalValue);
    FForAll::For<unsigned int, startInter, endInter, step, TaylorMainStruct>( argc, argv,
                                                           nbParticles,TreeHeight, SubTreeHeight,
                                                           loader.getBoxWidth(), loader.getCenterOfBox(),particles,
                                                           energyD,totPhysicalValue);

    FForAll::For<unsigned int, start, end, step2, SphericalBlasMainStruct>( argc, argv,
                                                                  nbParticles,TreeHeight, SubTreeHeight,
                                                                  loader.getBoxWidth(), loader.getCenterOfBox(),particles,
                                                                  energyD,totPhysicalValue);
    FForAll::For<unsigned int, start, end, step2, SphericalBlockBlasMainStruct>( argc, argv,
                                                                  nbParticles,TreeHeight, SubTreeHeight,
                                                                  loader.getBoxWidth(), loader.getCenterOfBox(),particles,
                                                                  energyD,totPhysicalValue);

    return 0;
}
