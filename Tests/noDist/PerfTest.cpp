
/**
 * \file
 * \author Quentin Khan
 *
 * This program is used to run different performance tests for the various
 * algorithms that have been implemented for ScalFMM.
 *
 * See the PerfUtils.hpp file classes for some more in depth information. Run
 * with argument --help for usage information.
 */


#include <iostream>
#include <string>

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#include "PerfTestUtils.hpp"

#include "TreeLoaderFCheb.hpp"

#include "KernelLoaderFChebSym.hpp"

#include "AlgoLoaderThread.hpp"
#include "AlgoLoaderTask.hpp"
#include "AlgoLoaderCostZones.hpp"

/**
 * \brief Runs a generic sequence of actions to use an algorithm.
 *
 * This function runs the basic steps that are needed to run an FMM algorithm
 * over a set of particles. It does the following steps :
 *
 *    - Load a tree using the class defined as a TreeLoader
 *    - Prepares the needed kernels using the KernelLoader
 *    - Prepares and runs the algorithm using the AlgorithmLoader
 *
 * See documentation of FTreeLoader, FKernelLoader, FAlgoLoader.
 */
template <class TreeLoader,
          template <typename TL> class KernelLoader,
          template <typename TL, template <typename TL> class KL> class AlgoLoader>
void runperf(FPerfTestParams& params)
{
    TreeLoader treeLoader(params);
    KernelLoader<TreeLoader> kernelLoader(params, treeLoader);
    AlgoLoader<TreeLoader, KernelLoader> algoLoader(params, treeLoader, kernelLoader);
    algoLoader.run();
}

int main (int argc, char** argv)
{
    FHelpDescribeAndExit(argc, argv,
                         "Driver for Chebyshev interpolation kernel  (1/r kernel).",
                         FParameterDefinitions::InputFile,
                         FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::NbThreads,
                         {{"--algo"},"Algorithm to run (costzones, basic, task)"},
                         {{"--schedule"},"OpenMP scheduling policy (static, dynamic)."});
    FPerfTestParams params;
    {
        using namespace FParameterDefinitions;
        using namespace FParameters;
        params.filename = getStr(argc,argv,InputFile.options,
                                 "../Data/unitCubeXYZQ100.bfma");
        params.treeHeight = getValue(argc, argv, OctreeHeight.options, 5);
        params.subTreeHeight = getValue(argc, argv, OctreeSubHeight.options, 2);
        params.nbThreads = getValue(argc, argv, NbThreads.options, 1);
        params.algo = getStr(argc,argv,{"--algo"},"task");
        params.omp_static_schedule =
            getStr(argc,argv,{"--schedule"},"dynamic") == std::string("static");
    }


    if( "basic-dynamic" == params.algo ) {
        runperf<TreeLoaderFCheb<>, KernelLoaderFChebSym, AlgoLoaderThread>(params);
    } else if( "task" == params.algo ) {
        runperf<TreeLoaderFCheb<>, KernelLoaderFChebSym, AlgoLoaderTask>(params);
    } else if ( "costzones" == params.algo ) {
        runperf<TreeLoaderFCheb<>, KernelLoaderFChebSym, AlgoLoaderCostZones>(params);
    } else {
        std::cout << "Unknown algorithm: " << params.algo << std::endl;
    }
    

}