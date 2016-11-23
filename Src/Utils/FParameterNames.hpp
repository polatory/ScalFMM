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
#ifdef FPARAMETERNAMES_HPP
#error FParameterNames must be included only once by each main file!
#else
#define FPARAMETERNAMES_HPP

#include "Utils/FGlobal.hpp"

#include "Utils/FParameters.hpp"

#include <iostream>
#include <vector>
#include <chrono>
#include <unistd.h>

/**
 * This file contains some useful classes/functions to manage the parameters,
 * but it also contains all the parameters definitions!
 * A Scalfmm executable must define a static object to propose a parameter to the user.
 */

/** A parameter descriptor */
struct FParameterNames {
    std::vector<const char*> options;
    const char* description;
};

/** The namespace where we put all the definitions */
namespace FParameterDefinitions {

static const FParameterNames Help = {
    {"-help", "--help"} ,
     "To have print the options used by the application."
};

static const FParameterNames Compile = {
    {"-show-compile", "--show-compile", "--flags"} ,
     "To have the list of flags and lib linked to scalfmm."
};

static const FParameterNames DateHost = {
    {"-show-info", "--show-host", "--datehost"} ,
     "To have to print the current host and the execution date."
};

static const FParameterNames UserParams = {
    {"-show-param", "--show-param", "-show-params", "--show-params"} ,
     "To print out the paremeters passed to the command line."
};

static const FParameterNames NbParticles = {
    {"-nb", "--number-of-particles", "-N"} ,
     "The number of particles if they are generated by the executable."
};

static const FParameterNames OctreeHeight = {
    {"-h", "--height", "-depth"} ,
     "The number of levels in the octree (at least 2 for the root and the leaves)."
};

static const FParameterNames OctreeSubHeight = {
    {"-sh", "--sub-height", "-subdepth"} ,
     "The number of allocated levels in the sub octree."
};

static const FParameterNames InputFile = {
    {"-f", "-fin", "--input-filename", "-filename"} ,
     "To give an input file."
};

static const FParameterNames InputFileOne = {
    {"-f1", "-fin1", "--file-one"} ,
     "To give the first input file."
};

static const FParameterNames InputFileTwow = {
    {"-f2", "-fin2", "--file-two"} ,
     "To give the second input file."
};

static const FParameterNames InputBinFormat = {
    {"-binin", "-bininput", "--binary-input" } ,
     "To input is in binary format."
};

static const FParameterNames OutputFile = {
    {"-fout", "--output-filename"} ,
     "To give the output filename with extension .fma or  .bfma (binary fma format)."
};

static const FParameterNames OutputVisuFile = {
    {"-fvisuout"} ,
     "Specify the name (with extension) and the format of the file for visualization purpose.\n                     Available format are  vtk, vtp, cvs or cosmo. vtp is the default."
};
//static const FParameterNames FormatVisuFile{
//	        {"-visufmt","-visu-fmt"},
//	        "To specify format for the visu file (vtk, vtp, cvs or cosmo). vtp is the default"
//	    };


static const FParameterNames OutputBinFormat = {
    {"-binout", "-binoutput"} ,
     "To output in binary format."
};

static const FParameterNames NbThreads = {
    {"-t", "-nbthreads"} ,
     "To choose the number of threads."
};

static const FParameterNames SequentialFmm = {
    {"-sequential", "--sequential-fmm"} ,
     "No parallelization in the FMM algorithm."
};

static const FParameterNames TaskFmm = {
    {"-task", "--task-fmm"} ,
     "Task parallelization in the FMM algorithm."
};

static const FParameterNames SHDevelopment = {
    {"-devp", "-sh-p"} ,
     "The degree of development for the spherical harmonic kernel (P)."
};

static const FParameterNames EnabledVerbose = {
    {"-verbose", "--verbose"} ,
     "To have a high degree of verbosity."
};

static const FParameterNames PeriodicityNbLevels = {
    {"-per", "--periodic-degree"} ,
     "The number of level upper to the root to proceed."
};

static const FParameterNames PeriodicityDisabled = {
    {"-noper", "--no-periodicity"} ,
     "To disable the periodicity."
};

static const FParameterNames DeltaT = {
    {"-dt", "--delta-time"} ,
     "The time step between iterations."
};

static const FParameterNames RotationKernel = {
    {"-rotation", "--rotation-kernel"} ,
     "To use the rotation kernel (based on spherical harmonics)."
};

static const FParameterNames SphericalKernel = {
    {"-spherical", "--spherical-kernel"} ,
     "To use the spherical harmonics old kernel."
};

static const FParameterNames ChebyshevKernel = {
    {"-chebyshev", "--chebyshev-kernel"} ,
     "To use the Chebyshev kernel."
};

static const FParameterNames Epsilon = {
    {"-epsilon", "--epsilon"} ,
     "The epsilon needed for the application."
};

static const FParameterNames PhysicalValue = {
    {"-pv", "--physical-value"} ,
     "The physical value of the particles."
};

static const FParameterNames SeparationCriterion = {
    {"-sep", "--separation-criterion"} ,
     "Specify number of clusters separing 2 well-separated clusters."
};

/** To print a list of parameters */
inline void PrintUsedOptions(const std::vector<FParameterNames>& options){
    std::cout << ">> Here is the list of the parameters you can pass to this application :\n";
    for(const FParameterNames& option : options ){
        std::cout << "\t";
        for(const char* name : option.options ){
            std::cout << name << ", ";
        }
        std::cout << "\n\t\t" << option.description << "\n";
        std::cout << "\n";
    }
}

inline bool CheckValidParameters(const int argc, char* argv[], const std::vector<FParameterNames> options){
    bool isValide = true;
    bool previousIsCorrectParameter = false;
    for(int idxParameter = 1 ; idxParameter < argc ; ++idxParameter){
        bool paramExist = false;
        for(unsigned idxTest = 0 ; idxTest < options.size() ; ++idxTest){
            paramExist = FParameters::existParameter(1, &argv[idxParameter], options[idxTest].options);
            if(paramExist){
                break;
            }
        }
        if(paramExist == true){
            previousIsCorrectParameter = true;
        }
        else if(previousIsCorrectParameter){
            previousIsCorrectParameter = false;
        }
        else{
            previousIsCorrectParameter = false;
            isValide = false;
            std::cout << "[PARAMETER-ERROR] Parameter " << (idxParameter-1) << " = \"" << argv[idxParameter] << "\" seems incorrect." << "\n";
        }
    }
    if(isValide == false){
        std::cout << "[ScalFMM] To know more about correct parameters ask for help by passing:\n\t";
        for(const char* param : FParameterDefinitions::Help.options){
            std::cout << "\t" << param << ",";
        }
        std::cout << std::endl;
    }
    return isValide;
}

inline void PrintFlags(){
    std::cout << "[ScalFMM] This executable has been compiled with:\n";
    std::cout << "    Flags: " << SCALFMMCompileFlags << "\n";
    std::cout << "    Libs;   " << SCALFMMCompileLibs   << "\n";
    std::cout << std::endl;     std::cout.flush();
}

inline void PrintDateHost(){
    std::cout << "[ScalFMM] This execution is on:\n";
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "    Date: " << std::ctime(&now) ;
    char hostname[1024];
    gethostname(hostname, 1024);
    std::cout << "    Host: " << hostname << "\n";
    std::cout << std::endl;     std::cout.flush();
}

inline void PrintGivenParams(int argc, const char* const * const argv){
    std::cout << "[ScalFMM] Given Parameters are:\n";
    for(int idx = 0 ; idx < argc ; ++idx){
        std::cout << " | " << argv[idx];
    }
    std::cout << " | \n\n";
    std::cout.flush();
}

}// End of namespace

/** This should be include at the beginin of all main file
 *  @code FHelpAndExit(argc, argv, FParameterDefinitions::NbParticles, FParameterNames OctreeSubHeight );
 */
#define FHelpAndExit(argc, argv, ...) \
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::Compile.options)) {\
        FParameterDefinitions::PrintFlags();\
    } \
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::DateHost.options)) {\
        FParameterDefinitions::PrintDateHost();\
    } \
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::UserParams.options)) {\
        FParameterDefinitions::PrintGivenParams(argc, argv);\
    } \
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::Help.options)) {\
        const std::vector<FParameterNames> optionsvec = {FParameterDefinitions::Compile, FParameterDefinitions::DateHost, \
                                                         FParameterDefinitions::UserParams, __VA_ARGS__};\
        FParameterDefinitions::PrintUsedOptions(optionsvec);\
        exit(0);                                              \
    } \
    if(FParameterDefinitions::CheckValidParameters(argc, argv, {FParameterDefinitions::Compile, FParameterDefinitions::DateHost, \
                    FParameterDefinitions::UserParams, __VA_ARGS__}) == false){ \
        return 121;                                                     \
    }                                                                   \

/** This should be include at the beginin of all main file
 *  @code FHelpDescribeAndExit(argc, argv,
 *  @code       "This executable is doing this and this.",
 *  @code       FParameterDefinitions::NbParticles, FParameterNames OctreeSubHeight );
 */
#define FHelpDescribeAndExit(argc, argv, description, ...) \
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::Compile.options)) {\
        FParameterDefinitions::PrintFlags();\
    } \
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::DateHost.options)) {\
        FParameterDefinitions::PrintDateHost();\
    } \
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::UserParams.options)) {\
        FParameterDefinitions::PrintGivenParams(argc, argv);\
    } \
    if(FParameters::existParameter(argc, argv, FParameterDefinitions::Help.options)) { \
        std::cout << argv[0] << " : " << description << "\n";           \
        const std::vector<FParameterNames> optionsvec = {FParameterDefinitions::Compile, FParameterDefinitions::DateHost, \
                                                         FParameterDefinitions::UserParams, __VA_ARGS__}; \
        FParameterDefinitions::PrintUsedOptions(optionsvec);            \
        exit(0);                                                          \
    }                                                                   \
    if(FParameterDefinitions::CheckValidParameters(argc, argv, {FParameterDefinitions::Compile, FParameterDefinitions::DateHost, \
                                         FParameterDefinitions::UserParams, __VA_ARGS__}) == false){ \
        return 121;\
    } \
    {\
        std::cout << "[ScalFMM] To have the help for this executable pass: "; \
        for(auto pr: FParameterDefinitions::Help.options) std::cout << pr << ", "; \
        std::cout << "\n";\
    }


#endif // FPARAMETERNAMES_HPP
