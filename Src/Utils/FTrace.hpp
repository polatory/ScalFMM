#ifndef FTRACE_HPP
#define FTRACE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FGlobal.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTrace
* Please read the license
*
* This class is used to print trace data during processing.
* You have to use the FTRACE macro as shown in the example.
*
* <code>
* FTRACE( FTrace::FRegion regionTrace("Sub part of function", __FUNCTION__ , __FILE__ , __LINE__) );
* ...
* FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
* </code>
*/

#ifndef SCALFMM_USE_TRACE

#define FTRACE( X )

class FTrace {
public:
    class FRegion {
    public:
        FRegion(const char*, const char*, const char* = 0, const int = 0){}
        void end(){}
    };
    class FFunction {
    public:
        FFunction(const char*, const char*, const char* = 0, const int = 0){}
        void end(){}
    };
};


#else

#define FTRACE( X ) X

#ifdef SCALFMM_USE_ITAC

#include <VT.h>

class FTrace {
public:
    class FRegion {
        VT_Region region;
    public:
        FRegion(const char*const regionName, const char*const className)
            : region( regionName, className) {}
        FRegion(const char*const regionName, const char*const className, const char* const file, const int line)
            : region( regionName, className, file, line ) {}
        void end(){
            region.end();
        }
    };

    class FFunction {
        VT_Function function;
    public:
        FFunction(const char*const functionName, const char*const className)
            : function( functionName, className) {}
        FFunction(const char*const functionName, const char*const className, const char* const file, const int line)
            : function( functionName, className, file, line ) {}
        void end(){
            function.end();
        }
    };
};

#else

#include <iostream>
#include <iomanip>

#include "FTic.hpp"

class FTrace{
    static int Deep;
    static FTic TimeSinceBegining;

    static void PrintTab(){
        std::cout << "{" << std::setw( 6 ) << TimeSinceBegining.tacAndElapsed() << "s} ";
        for(int idxDeep = 0 ; idxDeep < Deep ; ++idxDeep){
            std::cout << '\t';
        }
    }

public:
    class FRegion {
        bool closed;
        void close(){
            if(!closed){
                closed = true;
                --FTrace::Deep;
            }
        }
    public:
        FRegion(const char*const regionName, const char*const className)
                : closed(false) {
            FTrace::PrintTab();
            std::cout << "@Region: " << regionName << " (" << className << ")\n";
            ++FTrace::Deep;
        }
        FRegion(const char*const regionName, const char*const className, const char* const file, const int line)
                : closed(false) {
            FTrace::PrintTab();
            std::cout << "@Region: " << regionName << " (" << className << ")" << " -- line " << line << " file " << file << "\n";
            ++FTrace::Deep;
        }
        ~FRegion(){
            close();
        }
        void end(){
            close();
        }
    };

    class FFunction {
        bool closed;
        void close(){
            if(!closed){
                closed = true;
                --FTrace::Deep;
            }
        }
    public:
        FFunction(const char*const functionName, const char*const className)
                : closed(false){
            FTrace::PrintTab();
            std::cout << "@Function: " << functionName << " (" << className << ")\n";
            ++FTrace::Deep;
        }
        FFunction(const char*const functionName, const char*const className, const char* const file, const int line)
                : closed(false) {
            FTrace::PrintTab();
            std::cout << "@Function: " << functionName << " (" << className << ")" << " -- line " << line << " file " << file << "\n";
            ++FTrace::Deep;
        }
        ~FFunction(){
            close();
        }
        void end(){
            close();
        }
    };

    friend class FRegion;
    friend class FFunction;
};

#endif //SCALFMM_USE_ITAC

#endif //SCALFMM_USE_TRACE

#endif //FTRACE_HPP

// [--LICENSE--]
