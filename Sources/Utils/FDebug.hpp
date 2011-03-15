#ifndef FDEBUG_HPP
#define FDEBUG_HPP
// /!\ Please, you must read the license at the bottom of this page

// Uncomment next line to use debug mode
#define FUSE_DEBUG

#ifndef FUSE_DEBUG
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
class FDebug{
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
                this->stream->flush();
                if(this->stream != &std::cout) delete(this->stream);
	}

	/**
	* Copy constructor forbiden
	* @param other the source class to copy
	*/
        FDebug(const FDebug& ){}

	/**
	* Copy constructor forbiden
	* @param other the source class to copy
	* @return this a reference to the current class
	*/
        FDebug& operator=(const FDebug& ){return *this;}


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
                (*this->stream) << inMessage;
		return *this;
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

#endif //FUSE_DEBUG

#endif //FDEBUG_HPP

// [--LICENSE--]
