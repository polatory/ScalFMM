// See LICENCE file at project root

// ==== CMAKE =====
// @FUSE_BLAS
// ================

// Keep in private GIT
//

/**
 * @brief This executable execute the chebyshev interpolation for
 * differents orders and differentscompression epsilon (see Matthias
 * Messner's papers for more details.)
 */

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



template <class OctreeClass, class LeafClass>
void checkResAndPrint(OctreeClass * tree, FmaRWParticle<FReal, 8,8> * const particles,
                      FReal energyD, std::string name,int order,double timeSpent, FReal epsilon){
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
    std::ostringstream ostr;
    ostr << order;
    ostr << name;
    outfile.open(ostr.str(), std::ios_base::app);
    if(outfile.tellp() == 0){
        outfile << name+"_Order \t" << name+"_Epsilon" <<"\t" << name+"_Energy \t"
                << name+"_potential_L2 \t" << name+"_potential_RMS \t" << name+"_potential_Inf \t"
                << name+"_fx_L2 \t" << name+"_fx_RMS \t" << name+"_fx_Inf \t"
                << name+"_fy_L2 \t" << name+"_fy_RMS \t" << name+"_fy_Inf \t"
                << name+"_fz_L2 \t" << name+"_fz_RMS \t" << name+"_fz_Inf \t"
                << name+"Time \n";
    }
    outfile << order << "\t" << epsilon <<"\t"<< FMath::Abs(energy-energyD) /energyD << "\t"
            << potentialDiff.getL2Norm() << "\t" << potentialDiff.getRMSError() << "\t" << potentialDiff.getInfNorm() << "\t"
            << fx.getL2Norm() << "\t" << fx.getRMSError() << "\t" << fx.getInfNorm() << "\t"
            << fy.getL2Norm() << "\t" << fy.getRMSError() << "\t" << fy.getInfNorm() << "\t"
            << fz.getL2Norm() << "\t" << fz.getRMSError() << "\t" << fz.getInfNorm() << "\t" << timeSpent << "\n";
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
            typedef FSimpleLeaf<FReal,ContainerClass> LeafClass;
            typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
            typedef FChebCell<FReal,ORDER> CellClass;
            typedef FOctree<FReal,CellClass,ContainerClass,LeafClass> OctreeClass;

            typedef FChebSymKernel<FReal,CellClass,ContainerClass,MatrixKernelClass,ORDER> KernelClass;
            typedef FFmmAlgorithmThread<OctreeClass,CellClass,ContainerClass,KernelClass,LeafClass> FmmClass;
            // Create Matrix Kernel
            //Variation on compression espilon

            FReal powerOf10 = ORDER-2;
            while(powerOf10 <= ORDER+2)
                { // -----------------------------------------------------
                    // init oct-tree
                    OctreeClass tree(TreeHeight, SubTreeHeight,  BoxWidth,CenterOfBox);

                    FReal epsilon = FMath::pow(10.0,static_cast<FReal>(-powerOf10));
                    const MatrixKernelClass MatrixKernel;
                    { // -----------------------------------------------------
                        for(FSize idxPart = 0 ; idxPart < nbParticles; ++idxPart){
                            // put in tree
                            tree.insert(particles[idxPart].getPosition(), idxPart, particles[idxPart].getPhysicalValue());
                        }
                    } // -----------------------------------------------------

                    time.tic();
                    KernelClass* kernels = new KernelClass(TreeHeight, BoxWidth, CenterOfBox,&MatrixKernel,epsilon);
                    FmmClass algorithm(&tree, kernels);
                    algorithm.execute();
                    time.tac();
                    std::cout <<"(FChebSymKernel @Algorithm = " << time.elapsed() << " s). ORDER "<< ORDER << " Epsilon "<< epsilon
                              << std::endl;
                    delete kernels;

                    { // Check that each particle has been summed with all other
                        checkResAndPrint<OctreeClass,LeafClass>(&tree,particles,energyD,std::string("Chebyshev.res"),
                                                                ORDER,time.elapsed(),epsilon);
                    }
                    powerOf10 += 0.5;
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
    FmaRWParticle<FReal,8,8>* const particles = new FmaRWParticle<FReal,8,8>[nbParticles];
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
    const int startInter = 4;
    const int endInter = 14;
    const int step = 1;

    FForAll::For<unsigned int, startInter, endInter, step, ChebMainStruct>( argc, argv,
                                                                            nbParticles,TreeHeight, SubTreeHeight,
                                                                            loader.getBoxWidth(), loader.getCenterOfBox(),particles,
                                                                            energyD,totPhysicalValue);
    return 0;
}
