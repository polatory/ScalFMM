// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _LOADFMAANDRUNFMMARGS_HPP_
#define _LOADFMAANDRUNFMMARGS_HPP_

#include <string>
#include <omp.h>

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"


class loadFMAAndRunFMMArgs {
    const int _treeHeightInit = 5;
    const int _subTreeHeightInit = 1;
    const int _zoneCountInit = omp_get_max_threads();
    const int _verboseInit = 0;
    const char* _inFileNameInit = "";
    const char* _outFileNameInit = "balancetest";
    const char* _outFileNameExtInit = "csv";

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

public:
    int treeHeight() const {
        using namespace FParameterDefinitions;
        using namespace FParameters;

        return getValue(_argc, _argv, OctreeHeight.options, _treeHeightInit);
    }

    int subTreeHeight() const {
        using namespace FParameterDefinitions;
        using namespace FParameters;

        return getValue(_argc, _argv, OctreeSubHeight.options,
                        _subTreeHeightInit);
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

        return getStr(_argc, _argv, OutputFileBasename.options, _outFileNameInit);
    } 

    std::string outFileExt() const {
        using namespace FParameterDefinitions;
        using namespace FParameters;
        
        std::string ext = getStr(_argc, _argv, OutputFileExtension.options,
                                 _outFileNameExtInit);
        if ( ext.at(0) != '.' )
            return '.' + ext;
        return ext;
    } 

    loadFMAAndRunFMMArgs(int argc, char** argv) : _argc(argc), _argv(argv) {
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
             OctreeSubHeight,
             InputFile,
             OutputFileBasename,
             OutputFileExtension,
             ZoneCount,
             EnabledVerbose
             );
        return 0;
    }

};


#endif
