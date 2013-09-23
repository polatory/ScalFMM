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
#ifndef FDEBUG_HPP
#define FDEBUG_HPP
#include <iostream>


#include "FGlobal.hpp"
#include "FNoCopyable.hpp"

#ifndef ScalFMM_USE_LOG

#define FDEBUG( X )

#else

#define FDEBUG( X ) X

#include <iostream>
#include <fstream>
#include <sstream>

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FDebug
* Please read the license
*
* This class is used to print debug data durint processing.
* You have to use the DEBUG macro as shown in the example.
*
* <code>
* FDEBUG( FDebug::Controller.writeFromLine("hello World", __LINE__, __FILE__); ) <br>
* FDEBUG( FDebug::Controller << "I want to debug " << toto << "\n"; )
* </code>
*/
class FDebug : public FNoCopyable {
private:
	std::ostream* stream;	//< Standart c++ ostream

	/** Default constructor forbiden */
        FDebug() : stream(&std::cout) {
	}

	/** Default destructor forbiden */
	virtual ~FDebug(){
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

public:
	static FDebug Controller; 	//< Singleton

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

	/**
	* stream operator to print debug data
	* @param inMessage a message - from any type - to print
	* @return current FDebug
	*/
	template <class T>
	FDebug& operator<<(const T& inMessage){
                return write(inMessage);
	}

	/**
	* to write debug data
	* @param inMessage a message - from any type - to print
	* @return current FDebug
	*/
	template <class T>
	FDebug& write(const T& inMessage){
                (*this->stream) << inMessage;
		return *this;
	}

        /** Flush data into stream */
        void flush(){
            this->stream->flush();
        }

        enum FlushType{
            Flush,
            FlushWithLine
        };

        /**
        * stream operator to flush debug data
        * @param inType flush type
        * @return current FDebug
        */
        FDebug& write(const FlushType inType){
            if(inType == FlushWithLine) (*this->stream) << '\n';
            flush();
            return *this;
        }

	/**
	* to write debug data with line & file
	* @param inMessage a message - from any type - to print
	* @param inLinePosition line number
	* @param inFilePosition file name
	* @return current FDebug
	*
        * <code> FDebug::Controller.writeFromLine("hello World", __LINE__, __FILE__); </code>
	*
	* To prevent use from multiple thread we use a ostringstream before printing
	*/
	template <class T, class Tline, class Tfile>
	FDebug& writeFromLine(const T& inMessage, const Tline& inLinePosition, const Tfile& inFilePosition){
		std::ostringstream oss;
		oss << "Message from " << inFilePosition << " (at line " << inLinePosition <<")\n";
		oss << ">> " << inMessage << "\n";

                (*this->stream) << oss.str();
		return *this;
	}

	/**
	* to write debug data with line & file
	* @param inVariable variable name
	* @param inValue variable value
	* @param inLinePosition line number
	* @param inFilePosition file name
	* @return current FDebug
	*
        * <code> FDebug::Controller.writeVariableFromLine( "toto", toto, __LINE__, __FILE__); </code>
	*
	* To prevent use from multiple thread we use a ostringstream before printing
	*/
	template <class T, class Tline, class Tfile>
	FDebug& writeVariableFromLine(const char* const inVariable, const T& inValue, const Tline& inLinePosition, const Tfile& inFilePosition){
		std::ostringstream oss;
		oss << "[Value] " << inVariable << " = " << inValue << " at line " << inLinePosition <<" (file " << inFilePosition << ")\n";

                (*this->stream) << oss.str();
		return *this;
	}

};

#endif //SCALFMM_USE_DEBUG

#endif //FDEBUG_HPP


