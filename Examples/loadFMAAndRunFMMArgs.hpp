#ifndef _LOADFMAANDRUNFMMARGS_HPP_
#define _LOADFMAANDRUNFMMARGS_HPP_

#include "tclap++11/CmdLine.h"
#include "tclap++11/LinuxOutput.hpp"
#include "tclap++11/CompletionVisitor.h"

class loadFMAAndRunFMMArgs {
    using vs = std::vector<std::string>;
    
    const int _treeHeightInit = 5;
    const int _subTreeHeightInit = 1;
    const int _zoneCountInit = 4;
    const std::string _inFileNameInit = "";
    const std::string _outFileNameInit = "balancetest";
    const std::string _outFileNameExtInit = "csv";

    TCLAP::CmdLine
    _cmd {  "Loads an FMA file into a tree and runs a pseudo FMM algorithm "
            "through it to compute load balancing.", ' ', "0.0"};
    
    TCLAP::LinuxOutput _cmdFormat;
    
    TCLAP::BashCompletion::Visitor _compVistor;
    TCLAP::SwitchArg
    _compArg { "",
            "completion",
            "Show completion arguments",
            _cmd, false, &_compVistor};


    TCLAP::ValueArg <int>
    _subTreeHeight { "", vs{"subtree-height", "sth"},
            "Subtree height.", false, _subTreeHeightInit, "sub-tree height", _cmd};
    
    TCLAP::ValueArg <int>
    _treeHeight { "H", vs{"tree-height", "th"},
            "Tree height.", true, _treeHeightInit, "tree height", _cmd};

    TCLAP::ValueArg <int>
    _zoneCount { "z", vs{"zone-count"},
            "Number of zones to create", true, _zoneCountInit, "zone count", _cmd};

    TCLAP::ValueArg <std::string>
    _outFileExt { "x", "output-file-extension",
            "Output files extension. One file is created for each level in the"
            " tree. Each file has a 'basename_$nbZones$z.$i$.extension' "
            "extension where $i is the level. Default value is "
            + _outFileNameExtInit + ".",
            false, _outFileNameExtInit, "suffix", _cmd};

    TCLAP::ValueArg <std::string>
    _outFileName { "o", "output-file-basename",
            "Output files' basename. One file is created for each level in "
            "the tree. Each file has a level-in-tree based extension.",
            false, _outFileNameInit, "basename", _cmd};
    
    TCLAP::ValueArg <std::string>
    _inFileName { "f", "input-file",
            "Input file name.", true, _inFileNameInit, "filename", _cmd};
    
public:
    int treeHeight()    {return _treeHeight.getValue();}
    int subTreeHeight() {return _subTreeHeight.getValue();}
    int zoneCount() {return _zoneCount.getValue();}
    std::string inFileName() {return _inFileName.getValue();} 
    std::string outFileName() {return _outFileName.getValue();} 
    std::string outFileExt() {
        std::string ext = _outFileExt.getValue();
        if ( ext.at(0) != '.' )
            return '.' + ext;
        return ext;
    } 

    loadFMAAndRunFMMArgs(int argc, char** argv) 
        : _compVistor(argc, argv, &_cmd, {
                { &_inFileName,{TCLAP::BashCompletion::FILE, "!*.fma" }},
                {&_outFileExt, {"csv"}}}
            ){
        int columns = 80;
        struct winsize w;
        if (ioctl(0, TIOCGWINSZ, &w) == 0)
            columns = w.ws_col;

        _compArg.setHideDesc(true);
        _cmdFormat.setTextWidth(columns);                        
        _cmd.setOutput(&_cmdFormat);
        _cmd.parse(argc,argv);
    }
};


#endif
