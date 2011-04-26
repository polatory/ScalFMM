#ifndef FTRACE_HPP
#define FTRACE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FGlobal.hpp"

#ifndef SCALFMM_USE_TRACE

#define FTRACE( X )

#else

#define FTRACE( X ) X

#include <iostream>
#include <fstream>
#include <sstream>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTrace
* Please read the license
*
* This class is used to print trace data durint processing.
* You have to use the FTRACE macro as shown in the example.
*
* <code>
* FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
* ...
* FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
* </code>
*/
class FTrace{
private:    
    ///////////////////////////////////////////////////
    // Constant & config
    ///////////////////////////////////////////////////

    static const int MaxDeep = 4;               //< Max deep enabled
    static const int PrintParameters = true;    //< Enable parameters printing
    static const int PrintLineNumber = true;    //< Enable parameters printing
    static const int PrintFileName = true;      //< Enable parameters printing

    /** This tags are used to enable/disable
      * trace depending on the functions that
      * calls it.
      */
    enum Tags {
        FMM = 1,
        FILE = 2,
        KERNELS = 4,
        UTILS = 8,
        UNDEFINED = 16
    };

    static const int EnabledTags = FMM | FILE ;  //< Which tag can print

    ///////////////////////////////////////////////////
    // Real class
    ///////////////////////////////////////////////////

    int currentDeep;            //< Current deep in the code
    std::ostream* stream;	//< Standart c++ ostream

    /** Default constructor forbiden */
    FTrace() : stream(&std::cout), currentDeep(0) {
    }

    /** Default destructor forbiden */
    virtual ~FTrace(){
        close();
    }

    /**
	* Close the current debug stream
	* dealloc the stream if differents from cout.
	* after this call stream is useless
	*/
    void close(){
        flush();
        if(this->stream != &std::cout) delete(this->stream);
    }

    /**
	* Copy constructor forbiden
	* @param other the source class to copy
	*/
    FTrace(const FTrace& ){}

    /**
	* Copy constructor forbiden
	* @param other the source class to copy
	* @return this a reference to the current class
	*/
    FTrace& operator=(const FTrace& ){
        return *this;
    }

    /** Print tab depending on the deep
      *
      */
    void printTab(){
        for(int idxDeep = 0 ; idxDeep < currentDeep ; ++idxDeep ){
            (*this->stream) << "\t";
        }
    }

public:
    static FTrace Controller; 	//< Singleton

    /** Call this function when entering a function
      * @param inTag the tag at this position
      * @param inName the function's name
      * @param inFile the file's name
      * @param inLine the line number
      * <code>
      * FTRACE( FTrace::Controller.enterFunction(FTrace::FMM, __FUNCTION__ , __FILE__ , __LINE__) );
      * </code>
      */
    template <class T, class TF, class TL>
    FTrace& enterFunction(const Tags inTag, const T& inName, const TF& inFile , const TL& inLine){
        if(inTag & EnabledTags){
            if(currentDeep < MaxDeep){
                printTab();
                (*this->stream) << "@Method " << inName << "\n";
                ++currentDeep;
                if(PrintLineNumber){
                    printTab();
                    (*this->stream) << "@line " << inLine << "\n";
                }
                if(PrintFileName){
                    printTab();
                    (*this->stream) << "@file " << inFile << "\n";
                }
            }
            else{
                ++currentDeep;
            }
        }
        return *this;
    }

    /** Call this function when leaving a function
      * @param inTag the tag at this position
      * @param inName the function's name
      * <code>
      * FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM, __FUNCTION__ ) );
      * </code>
      */
    template <class T>
    FTrace& leaveFunction(const Tags inTag, const T& inName){
        if(inTag & EnabledTags){
            --currentDeep;
            if(currentDeep < MaxDeep){
                printTab();
                (*this->stream) << "@end " << inName << "\n";
            }
        }
        return *this;
    }

    /** Call this function when leaving a function
      * @param inTag the tag at this position
      * <code>
      * FTRACE( FTrace::Controller.leaveFunction(FTrace::FMM) );
      * </code>
      */
    FTrace& leaveFunction(const Tags inTag){
        if(inTag & EnabledTags){
            --currentDeep;
            if(currentDeep < MaxDeep){
                printTab();
                (*this->stream) << "@end\n";
            }
        }
        return *this;
    }

    /** Call this function when leaving a function
      * @param inTag the tag at this position
      * @param inName the parameter name
      * @param inValue the parameter value
      * <code>
      * FTRACE( FTrace::Controller.printParameter(FTrace::FMM, "toto" , toto) );
      * </code>
      */
    template <class T, class TV>
    FTrace& printParameter(const Tags inTag, const T& inName, const TV& inValue){
        if(inTag & EnabledTags){
            if(PrintParameters && currentDeep <= MaxDeep){
                printTab();
                (*this->stream) << "@param " << inName << "\t = " << inValue << ")\n";
            }
        }
        return *this;
    }


    /**
	* To set the debug stream to write into a file
	* @param filename the file to write
	*/
    void writeToFile(const char* const filename){
        close();

        std::ofstream* const file = new std::ofstream();
        file->open(filename);

        this->stream = file;
    }

    /**
	* To set the debug stream to write to std::cout
	*/
    void writeToCout(){
        close();
        this->stream = &std::cout;
    }

    /** Flush data into stream */
    void flush(){
        this->stream->flush();
    }

};

#endif //SCALFMM_USE_TRACE

#endif //FTRACE_HPP

// [--LICENSE--]
