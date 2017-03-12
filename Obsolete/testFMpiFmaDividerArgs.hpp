// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _TESTFMPIFMADIVIDERARGS_HPP_
#define _TESTFMPIFMADIVIDERARGS_HPP_

#include <string>
#include <omp.h>

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"




class testFMpiFmaDividerArgs {
    const int _treeHeightInit = 5;
    const int _zoneCountInit = omp_get_max_threads();
    const int _verboseInit = 0;
    const char* _inFileNameInit = "";
    const char* _outFileNameInit = "";
    const char* _outFileNameExtInit = "";
    const char* _dispatchPolicyInit = "particles";

    int _argc;
    char** _argv;

    const FParameterNames OutputFileBasename =
    {{"--output-file-basename", "-fout-base"},
     "Output files' basename. One file is created for each level in "
     "the tree. Each file has a level-in-tree based extension."};

    const FParameterNames OutputFileExtension =
    {{"--output-file-extention", "-fout-ext"},
     "Output files extension. One file is created for each level in the"
     " tree. Each file has a 'basename_$nbZones$z.$i$.extension' "
     "extension where $i is the level."};

    const FParameterNames ZoneCount = 
    {{"--zone-count","-z"},"Number of zones to create."};

    const FParameterNames DispatchPolicy = 
    {{"--dispatch","-d"},"What should be counted to dispatch the work"
     " ('particles' or 'boxes')."};

public:
    int treeHeight() const {
        using namespace FParameterDefinitions;
        using namespace FParameters;

        return getValue(_argc, _argv, OctreeHeight.options, _treeHeightInit);
    }

    int zoneCount()     const {
        using namespace FParameterDefinitions;
        using namespace FParameters;

        return getValue(_argc, _argv, ZoneCount.options, _zoneCountInit);
    }

    int verboseLevel()  const {
        using namespace FParameterDefinitions;
        using namespace FParameters;

        return getValue(_argc, _argv, EnabledVerbose.options, _verboseInit);
    }

    std::string inFileName()  const {
        using namespace FParameterDefinitions;
        using namespace FParameters;

        return getStr(_argc, _argv, InputFile.options, _inFileNameInit);
    } 

    std::string outFileName() const {
        using namespace FParameterDefinitions;
        using namespace FParameters;

        if( ! existParameter(_argc, _argv, OutputFileBasename.options) ) {
            std::string str(getStr(_argc, _argv, InputFile.options, _inFileNameInit));
            return str.substr(0, str.find_last_of("."));
        }

        return getStr(_argc, _argv, OutputFileBasename.options, _outFileNameInit);
    } 

    std::string outFileExt() const {
        using namespace FParameterDefinitions;
        using namespace FParameters;

        if( ! existParameter(_argc, _argv, OutputFileExtension.options) ) {
            std::string str(getStr(_argc, _argv, InputFile.options, _inFileNameInit));
            return str.substr(str.find_last_of("."));
        }


        
        std::string ext = getStr(_argc, _argv, OutputFileExtension.options,
                                 _outFileNameExtInit);
        if ( ext.at(0) != '.' )
            return '.' + ext;
        return ext;
    } 

    std::string dispatchPolicyString() const {
        using namespace FParameterDefinitions;
        using namespace FParameters;

        std::string str = getStr(_argc, _argv, DispatchPolicy.options, _dispatchPolicyInit);
        return str;
    }

    int dispatchPolicy() const {
        std::string str = dispatchPolicyString();

        if ( "particles" == str ) {
            return 0;
        } else {
            return 1;
        }
    }

    testFMpiFmaDividerArgs(int argc, char** argv) : _argc(argc), _argv(argv) {
        parse();
    }

    int parse() {
    using namespace FParameterDefinitions;
    using namespace FParameters;

        FHelpDescribeAndExit
            (_argc, _argv,
             "Loads an FMA file into a tree and runs a pseudo FMM algorithm "
             "through it to compute load balancing.",
             OctreeHeight,
             InputFile,
             OutputFileBasename,
             OutputFileExtension,
             ZoneCount,
             DispatchPolicy

             );
        return 0;
    }

};


#endif
