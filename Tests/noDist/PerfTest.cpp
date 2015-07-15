// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE


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

#include "PerfTest/PerfTestUtils.hpp"

#include "PerfTest/TreeLoaderBasic.hpp"
#include "PerfTest/TreeLoaderFCheb.hpp"
#include "PerfTest/TreeLoaderMpiSplitFCheb.hpp"
#include "PerfTest/TreeLoaderMpiGenericFCheb.hpp"

#include "PerfTest/KernelLoaderFChebSym.hpp"

#include "PerfTest/AlgoLoaderThread.hpp"
#include "PerfTest/AlgoLoaderTask.hpp"
#include "PerfTest/AlgoLoaderSectionTask.hpp"
#include "PerfTest/AlgoLoaderCostZones.hpp"
#include "PerfTest/AlgoLoaderThreadBalance.hpp"
#include "PerfTest/AlgoLoaderThreadProc.hpp"

#define HOST_NAME_MAX 64

/**
 * \brief Runs a generic sequence of actions to use an algorithm.
 *
 * This function runs the basic steps that are needed to run an FMM algorithm
 * over a set of particles. It does the following steps :
 *
 *    - Load a tree using the class defined as a TreeLoader
 *    - Prepare the needed kernels using the KernelLoader
 *    - Prepare and run the algorithm using the AlgorithmLoader
 *
 * See documentation of FTreeLoader, FKernelLoader, FAlgoLoader.
 */
template <class TreeLoader,
          template <typename TL_1> class KernelLoader,
          template <typename TL_2, template <typename TL_3> class KL> class AlgoLoader>
void runperf(FPerfTestParams& params)
{
    TreeLoader treeLoader(params);
    KernelLoader<TreeLoader> kernelLoader(params, treeLoader);
    AlgoLoader<TreeLoader, KernelLoader> algoLoader(params, treeLoader, kernelLoader);
    algoLoader.run();

    char hostname[HOST_NAME_MAX];
    memset(hostname,'\0',HOST_NAME_MAX);
    if ( -1 == gethostname(hostname, HOST_NAME_MAX-1) ) {
        perror("Could not get hostname");
        strncpy(hostname, "unknown", HOST_NAME_MAX);
    }

    std::cout << "@@ "
              << "host:" << hostname << " "
              << "algo:" << params.algo << " "
              << "file:" << params.filename.substr(
                  params.filename.find_last_of('/')+1 ) << " "
              << "particles:" << treeLoader._loader.getNumberOfParticles() << " "
              << "procs:"     << params.nbProcs                            << " "
              << "threads:"   << params.nbThreads                          << " "
              << "height:"    << params.treeHeight                         << " "
              << "subheight:" << params.subTreeHeight                      << " "
              << algoLoader.getRunInfoString()
              << "P2M:" << algoLoader.getCumulatedTime(FAlgorithmTimers::P2MTimer)     << " "
              << "M2M:" << algoLoader.getCumulatedTime(FAlgorithmTimers::M2MTimer)     << " "
              << "M2L:" << algoLoader.getCumulatedTime(FAlgorithmTimers::M2LTimer)     << " "
              << "L2L:" << algoLoader.getCumulatedTime(FAlgorithmTimers::L2LTimer)     << " "
              << "P2PL2P:" << algoLoader.getCumulatedTime(FAlgorithmTimers::NearTimer) << " "
              << std::endl;
}

namespace ParName {
    const FParameterNames Algo = {{"--algo"},"Algorithm to run (basic, task, costzones, sectiontask, autobalance, mpi)"};
    const FParameterNames Schedule = {{"--schedule"},"OpenMP scheduling policy (static, dynamic)."};
    const FParameterNames ChunkSize = {{"--chunk-size"},"OpenMP chunk size for basic dynamic algorithm."};
}

int main (int argc, char** argv)
{
    // Parameter handling //////////////
    FHelpDescribeAndExit(argc, argv,
                         "Driver for Chebyshev interpolation kernel  (1/r kernel).",
                         FParameterDefinitions::InputFile,
                         FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::OctreeSubHeight,
                         FParameterDefinitions::NbThreads,
                         ParName::Algo,
                         ParName::Schedule,
                         ParName::ChunkSize);
    FPerfTestParams params;
    {
        using namespace FParameterDefinitions;
        using namespace FParameters;
        params.filename      = getStr(argc,argv,InputFile.options,
                                 "../Data/unitCubeXYZQ100.bfma");
        params.treeHeight    = getValue(argc, argv, OctreeHeight.options, 5);
        params.subTreeHeight = getValue(argc, argv, OctreeSubHeight.options, 2);
        params.nbThreads     = getValue(argc, argv, NbThreads.options, 1);
        params.algo = getStr(argc,argv,ParName::Algo.options,"task");
        params.omp_chunk_size = getValue(argc, argv, ParName::ChunkSize.options, 0);

#ifdef SCALFMM_USE_MPI
        std::string prefix("mpi-");
        if( params.algo.substr(0, prefix.size()) == prefix ) {
            params.mpiContext = new FMpi(argc,argv);
            params.nbProcs = params.mpiContext->global().processCount();
        }
#endif
    }
    // End of Parameter handling ///////

    char hostname[HOST_NAME_MAX];
    memset(hostname,'\0',HOST_NAME_MAX);
    if ( -1 == gethostname(hostname, HOST_NAME_MAX-1) ) {
        perror("Could not get hostname");
        strncpy(hostname, "unknown", HOST_NAME_MAX);
    }
    std::cout << "Hostname: " << hostname << std::endl;

    omp_set_num_threads(params.nbThreads);

    using FReal = double;
    constexpr const int ORDER = 7;

    if( "basic" == params.algo ) {
        runperf<TreeLoaderFCheb<FReal,ORDER>,
                KernelLoaderFChebSym,
                AlgoLoaderThread>
            (params);
    } else if( "task" == params.algo ) {
        runperf<TreeLoaderFCheb<FReal,ORDER>,
                KernelLoaderFChebSym,
                AlgoLoaderTask>
            (params);
    } else if ( "costzones" == params.algo ) {
        runperf<TreeLoaderFCheb<FReal,ORDER>,
                KernelLoaderFChebSym,
                AlgoLoaderCostZones>
            (params);
    } else if ( "sectiontask" == params.algo ) {
        runperf<TreeLoaderFCheb<FReal,ORDER>,
                KernelLoaderFChebSym,
                AlgoLoaderSectionTask>
            (params);
    } else if ( "autobalance" == params.algo ) {
        runperf<TreeLoaderFCheb<FReal,ORDER>,
                KernelLoaderFChebSym,
                AlgoLoaderThreadBalance>
            (params);
    } else if ( "mpi-split" == params.algo ) {
        runperf<TreeLoaderMpiSplitFCheb<FReal,ORDER>,
                KernelLoaderFChebSym,
                AlgoLoaderThreadProc>
            (params);        
    } else if ( "mpi-generic" == params.algo ) {
        runperf<TreeLoaderMpiGenericFCheb<FReal,ORDER>,
                KernelLoaderFChebSym,
                AlgoLoaderThreadProc>
            (params);        
    } else {
        std::cout << "Unknown algorithm: " << params.algo << std::endl;
    }
    
    if( nullptr != params.mpiContext ) {
        delete params.mpiContext;
    }

}
