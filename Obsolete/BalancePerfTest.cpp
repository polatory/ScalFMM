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

// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

// ==== CMAKE =====
// @FUSE_BLAS
// ================


/**
 * \file
 * \authors B. Bramas, O. Coulaud
 * \brief This program runs the balanced FMM Algorithm with the interpolation
 * kernel based on Chebyshev interpolation (1/r kernel)
 *
 * This program runs the FMM Algorithm with the Chebyshev kernel and compares
 * the results with a direct computation.
 *
 *
 *  This code is a short example to use the Chebyshev Interpolation approach for the 1/r kernel
 */


#include <iostream>
#include <string>

// Utilities
#include "ScalFmmConfig.h"
#include "Files/FFmaGenericLoader.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

// Data structures
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Containers/FOctree.hpp"
#include "Components/FSimpleLeaf.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

// Kernels
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#include "BalanceTree/FChebSymCostKernel.hpp"

// Algorithms
#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"
#include "Core/FFmmAlgorithmTask.hpp"
#include "BalanceTree/FFmmAlgorithmThreadBalanced.hpp"
#include "BalanceTree/FCostZones.hpp"


// typedefs
using FReal = double;

// Chebyshev accuracy
const unsigned int ORDER = 7;


class AbstractPerfTest {
protected:
    FTic time;

    template <typename... Args>
    struct false_type {
        bool value = false;
    };  

    template <typename...Args>
    void loadTree(Args...) {
        static_assert(false_type<Args...>::value,
                      "I don't know how to load this tree with this loader...");
    }

    template <class OctreeClass>
    void loadTree(FFmaGenericLoader<FReal>& loader, OctreeClass& tree) {
        std::cout << "Creating & inserting particles" << std::flush;

        time.tic();

        FPoint<FReal> position;
        FReal physicalValue = 0.0;
        for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart) {
            // Read particle per particle from file
            loader.fillParticle(&position,&physicalValue);
            // put particle in octree
            tree.insert(position, idxPart, physicalValue);
        }

        time.tac();
        std::cout << "Done  (" << time.elapsed() << "s)." << std::endl;
    }

    virtual void setup() = 0;
    virtual void runAlgo() = 0;
    virtual void finalize() = 0;
    
    template <class LeafClass, class OctreeClass, class FmmClass, class LoaderClass>
    void finalize(OctreeClass& tree, FmmClass& algo, LoaderClass& loader) {
        std::cout << "Timers Far Field \n"
                  << "P2M " << algo.getTime(FAlgorithmTimers::P2MTimer) << " seconds\n"
                  << "M2M " << algo.getTime(FAlgorithmTimers::M2MTimer) << " seconds\n"
                  << "M2L " << algo.getTime(FAlgorithmTimers::M2LTimer) << " seconds\n"
                  << "L2L " << algo.getTime(FAlgorithmTimers::L2LTimer) << " seconds\n"
                  << "P2P and L2P " << algo.getTime(FAlgorithmTimers::NearTimer) << " seconds\n"
                  << std::endl;

        std::cout << "Done  " << "(@Algorithm = " << time.elapsed() << " s) ." << std::endl;

        FSize N1 = 0, N2 = loader.getNumberOfParticles()/2, N3 = loader.getNumberOfParticles() - 1;
        FReal energy = 0.0;
        //
        //   Loop over all leaves
        //
        std::cout << std::endl;
        std::cout << std::scientific;
        std::cout.precision(10) ;

        tree.forEachLeaf([&](LeafClass* leaf){
                const FReal*const posX = leaf->getTargets()->getPositions()[0];
                const FReal*const posY = leaf->getTargets()->getPositions()[1];
                const FReal*const posZ = leaf->getTargets()->getPositions()[2];

                const FReal*const potentials = leaf->getTargets()->getPotentials();
                const FReal*const forcesX = leaf->getTargets()->getForcesX();
                const FReal*const forcesY = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
                const FSize nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();

                const FVector<FSize>& indexes = leaf->getTargets()->getIndexes();

                for(FSize idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const FSize indexPartOrig = indexes[idxPart];
                    if ((indexPartOrig == N1) || (indexPartOrig == N2) || (indexPartOrig == N3)  ) {
                        std::cout << "Index "<< indexPartOrig <<"  potential  " << potentials[idxPart]
                                  << " Pos "<<posX[idxPart]<<" "<<posY[idxPart]<<" "<<posZ[idxPart]
                                  << "   Forces: " << forcesX[idxPart] << " " << forcesY[idxPart] << " "<< forcesZ[idxPart] <<std::endl;
                    }
                    energy += potentials[idxPart]*physicalValues[idxPart] ;
                }
            });
        std::cout <<std::endl<<"Energy: "<< energy<<std::endl;        
    }   
    
public:
    virtual ~AbstractPerfTest(){};

    void run() {
        this->setup();
        this->runAlgo();
        this->finalize();
    }

};

template < template<typename...> class Algo > class PerfTest;

template <>
class PerfTest<FFmmAlgorithmThread> : public AbstractPerfTest {
public: // typedefs
    using CellClass          = FChebCell<FReal, ORDER>;
    using ContainerClass     = FP2PParticleContainerIndexed<FReal>;
    using LeafClass          = FSimpleLeaf<FReal, ContainerClass >;
    using OctreeClass        = FOctree<FReal, CellClass, ContainerClass, LeafClass>;
    using MatrixKernelClass  = FInterpMatrixKernelR<FReal>;
    using KernelClass        = FChebSymKernel      <FReal, CellClass, ContainerClass,
                                                    MatrixKernelClass, ORDER>;

    using FmmClass = FFmmAlgorithmThread<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>;

protected:    
    int _nbThreads;
    FFmaGenericLoader<FReal> _loader;
    OctreeClass _tree;
    FmmClass* _algo;

    bool _ompStaticScheduling;

public:
    PerfTest(const std::string& fileName, const int nbThreads, const int treeHeight, const int subTreeHeight, bool ompStaticScheduling) :
        _nbThreads(nbThreads) ,
        _loader(fileName),
        _tree(treeHeight, subTreeHeight,_loader.getBoxWidth(),_loader.getCenterOfBox()),
        _ompStaticScheduling(ompStaticScheduling) {
    }

    ~PerfTest() {
        if(_algo != nullptr)
            delete _algo;
    }

protected:
    virtual void setup() {
        omp_set_num_threads(_nbThreads);
        std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;

        loadTree(_loader,_tree);
    }

    virtual void runAlgo() {
        time.tic();
        const MatrixKernelClass MatrixKernel;
        KernelClass kernels(_tree.getHeight(), _loader.getBoxWidth(), _loader.getCenterOfBox(),&MatrixKernel);
        _algo = new FmmClass(&_tree, &kernels,_ompStaticScheduling);

        _algo->execute();
        time.tac();
    }

    void finalize() {
        AbstractPerfTest::finalize<LeafClass>(_tree, *_algo, _loader);
    }
};

template <>
class PerfTest<FFmmAlgorithmTask> : public AbstractPerfTest {
public: // typedefs
    using CellClass          = FChebCell<FReal, ORDER>;
    using ContainerClass     = FP2PParticleContainerIndexed<FReal>;
    using LeafClass          = FSimpleLeaf<FReal, ContainerClass >;
    using OctreeClass        = FOctree<FReal, CellClass, ContainerClass, LeafClass>;
    using MatrixKernelClass  = FInterpMatrixKernelR<FReal>;
    using KernelClass        = FChebSymKernel      <FReal, CellClass, ContainerClass,
                                                    MatrixKernelClass, ORDER>;

    using FmmClass = FFmmAlgorithmTask<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>;

protected:    
    int _nbThreads;
    FFmaGenericLoader<FReal> _loader;
    OctreeClass _tree;
    FmmClass* _algo;

public:
    PerfTest(const std::string& fileName, const int nbThreads, const int treeHeight, const int subTreeHeight) :
        _nbThreads(nbThreads) ,
        _loader(fileName),
        _tree(treeHeight, subTreeHeight, _loader.getBoxWidth(), _loader.getCenterOfBox()) {
    }

    ~PerfTest() {
        if(_algo != nullptr)
            delete _algo;
    }

protected:
    virtual void setup() {
        omp_set_num_threads(_nbThreads);
        std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;

        loadTree(_loader,_tree);
    }

    virtual void runAlgo() {
        time.tic();
        const MatrixKernelClass MatrixKernel;
        KernelClass kernels(_tree.getHeight(), _loader.getBoxWidth(), _loader.getCenterOfBox(),&MatrixKernel);
        _algo = new FmmClass(&_tree, &kernels);

        _algo->execute();
        time.tac();
    }

    void finalize() {
        AbstractPerfTest::finalize<LeafClass>(_tree, *_algo, _loader);
    }    
};



template <>
class PerfTest<FFmmAlgorithmThreadBalanced> : public AbstractPerfTest {
public: // typedefs
    using ContainerClass     = FP2PParticleContainerIndexed<FReal>;
    using CellClass          = FCostCell  <FChebCell<FReal, ORDER>>;
    using LeafClass          = FSimpleLeaf<FReal, ContainerClass >;
    using OctreeClass        = FOctree    <FReal, CellClass, ContainerClass, LeafClass>;
    using MatrixKernelClass  = FInterpMatrixKernelR<FReal>;
    using KernelClass        = FChebSymKernel      <FReal, CellClass, ContainerClass,
                                                    MatrixKernelClass, ORDER>;
    using CostKernelClass = FChebSymCostKernel     <FReal, CellClass, ContainerClass, 
                                                    MatrixKernelClass, ORDER, OctreeClass>;

    template <template <typename...> class T, typename Kernel>
    using FmmAlgoClass = T<OctreeClass, CellClass, ContainerClass, Kernel, LeafClass>;

    using FmmClass     = FmmAlgoClass<FFmmAlgorithmThreadBalanced, KernelClass>;
    using CostFmmClass = FmmAlgoClass<FFmmAlgorithm, CostKernelClass>;

    const FReal epsilon = 1e-4;

protected:    

    int _nbThreads;
    FFmaGenericLoader<FReal> _loader;
    OctreeClass _tree;
    FmmClass* _algo;

public:
    PerfTest<FFmmAlgorithmThreadBalanced>(
        const std::string& fileName, const int nbThreads,
        const int treeHeight, const int subTreeHeight) :
        _nbThreads(nbThreads) ,
        _loader(fileName),
        _tree(treeHeight, subTreeHeight, _loader.getBoxWidth(), _loader.getCenterOfBox()) {
    }

    ~PerfTest<FFmmAlgorithmThreadBalanced>() {
        if(_algo != nullptr)
            delete _algo;
    }

protected:
    virtual void setup() {
        omp_set_num_threads(_nbThreads);
        std::cout << "\n>> Using " << omp_get_max_threads() << " threads.\n" << std::endl;

        loadTree(_loader,_tree);
    }

    virtual void runAlgo() {
        // Compute tree cells costs
        CostKernelClass balanceKernel(&_tree, epsilon);
        CostFmmClass    costAlgo(&_tree, &balanceKernel);

        time.tic();
        costAlgo.execute();
        time.tac();
        std::cout << "Generating tree cost: " << time.elapsed() << "s.\n";

        FCostZones<OctreeClass, CellClass> costzones(&_tree, omp_get_max_threads());

        time.tic();
        costzones.run();
        time.tac();
        std::cout << "Generating cost zones: " << time.elapsed() << "s.\n";

        time.tic();
        const MatrixKernelClass MatrixKernel;
        KernelClass kernels(_tree.getHeight(), _loader.getBoxWidth(), _loader.getCenterOfBox(),&MatrixKernel);
        _algo = new FmmClass(&_tree, &kernels, costzones.getZoneBounds(), costzones.getLeafZoneBounds());

        _algo->execute();
        time.tac();
    }

    void finalize() {
        AbstractPerfTest::finalize<LeafClass>(_tree, *_algo, _loader);
    }
};



// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
    FHelpDescribeAndExit(argc, argv,
                         "Driver for Chebyshev interpolation kernel  (1/r kernel).",
                         FParameterDefinitions::InputFile,
                         FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::NbThreads,
                         {{"--algo"},"Algorithm to run (costzones, basic-static, basic-dynamic, task)"});

    const std::string  defaultFile("../Data/unitCubeXYZQ100.bfma" );
    const std::string  filename =
        FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, defaultFile.c_str());
    const unsigned int TreeHeight =
        FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 5);
    const unsigned int SubTreeHeight =
        FParameters::getValue(argc, argv, FParameterDefinitions::OctreeSubHeight.options, 2);
    const unsigned int NbThreads =
        FParameters::getValue(argc, argv, FParameterDefinitions::NbThreads.options, 1);
    const std::string  algoChoice =
        FParameters::getStr(argc,argv,{"--algo"},"costzones");

    std::cout << "file: " << filename << "  tree height: " << TreeHeight
              << "(" << SubTreeHeight << ")  algo: " << algoChoice << std::endl;

    if(algoChoice == "costzones") {
        PerfTest<FFmmAlgorithmThreadBalanced>
            balancePerfTest(filename, NbThreads, TreeHeight, SubTreeHeight);
        balancePerfTest.run();
    } else if (algoChoice == "basic-static") {
        PerfTest<FFmmAlgorithmThread>
            threadPerfTestStatic(filename, NbThreads, TreeHeight, SubTreeHeight, true);
        threadPerfTestStatic.run();
    } else if (algoChoice == "basic-dynamic") {
        PerfTest<FFmmAlgorithmThread>
            threadPerfTestDynamic(filename, NbThreads, TreeHeight, SubTreeHeight,false);
        threadPerfTestDynamic.run();
    } else if (algoChoice == "task") {
        PerfTest<FFmmAlgorithmTask>
            taskPerfTest(filename, NbThreads, TreeHeight, SubTreeHeight);
        taskPerfTest.run();
    } else {
        std::cerr << "Wrong algorithm choice. Try 'basic' or 'costzones'." << std::endl;
    }

}
