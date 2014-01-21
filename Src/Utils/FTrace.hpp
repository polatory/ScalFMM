// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FTRACE_HPP
#define FTRACE_HPP


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

#ifndef ScalFMM_USE_TRACE

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
        };
    };


#else

    #define FTRACE( X ) X

    #ifdef ScalFMM_USE_ITAC

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
            };
        };

    #else
        #ifdef ScalFMM_USE_EZTRACE

            #include <eztrace.h>

            class FTrace {
                static const unsigned IdModule = 0xCC00; // must be between 0x1000 and 0xff00
                static unsigned BuildMask(const char* phrase){
                    unsigned mask = 0;
                    if(phrase){
                        while( *phrase ){
                            mask = (mask<<1) ^ (*phrase++) ^ ((mask>>15)&1);
                        }
                    }
                    return ((mask & 0xFF) | IdModule) & ~0x2; // 0x[IdModule][Mask] last bits has to 01
                }
            public:
                class FRegion {
                    const unsigned mask;
                    bool hadFinished;
                public:
                    FRegion(const char*const regionName, const char*const className)
                        : mask(BuildMask(regionName)), hadFinished(false) {
                        EZTRACE_EVENT2(mask, regionName, className);
                    }
                    FRegion(const char*const regionName, const char*const className, const char* const file, const int line)
                        : mask(BuildMask(regionName)), hadFinished(false) {
                        EZTRACE_EVENT4(mask, regionName, className, file, line);
                    }
                    ~FRegion(){
                        end();
                    }
                    void end(){
                        if( !hadFinished ){
                            hadFinished = true;
                            EZTRACE_EVENT0(mask^0x3);
                        }
                    }
                };

                class FFunction {
                    const unsigned mask;
                public:
                    FFunction(const char*const functionName, const char*const className)
                        : mask(BuildMask(functionName)) {
                        EZTRACE_EVENT2(mask, functionName, className);
                    }
                    FFunction(const char*const functionName, const char*const className, const char* const file, const int line)
                        : mask(BuildMask(functionName)) {
                        EZTRACE_EVENT4(mask, functionName, className, file, line);
                    }
                    ~FFunction(){
                        EZTRACE_EVENT0(mask^0x3);
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
                };

                friend class FRegion;
                friend class FFunction;
            };

        #endif //ScalFMM_USE_EZTRACE

    #endif //ScalFMM_USE_ITAC

#endif //ScalFMM_USE_TRACE

#endif //FTRACE_HPP


