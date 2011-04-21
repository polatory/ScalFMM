#ifndef FABSTRACTAPPLICATION_HPP
#define FABSTRACTAPPLICATION_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <sstream>
#include <iostream>
#include <string.h>

#include "FAssertable.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractApplication
* Please read the license
*
* This class is an interface for main application.
* It represents the core of the system
* Only one app can be instancied by program.
*
* @warning you have to implement run() and call execute to start the app
*
* Please refere to testApplication.cpp to see an example.
* <code>
* </code>
*/
class FAbstractApplication {
private:	
	const int argc;			//< argc from command line

	char ** const argv;		//< argv from command line

protected:
	/**
	* This will be called as the main method
        * @warning Must be impleted in the derived class
	*/
	virtual void run() = 0;

	/**
	* This function is called before the run if the process is the master one
	*/
	virtual void initMaster(){}

	/**
	* This function is called before the run if the process is a slave
	*/
	virtual void initSlave(){}

        /**
          * Send data to another process
          */
        virtual void sendData(const int inReceiver, const int inSize, void* const inData, const int inTag) = 0;

        /**
          * Receive from any process
          */
        virtual void receiveData(const int inSize, void* const inData, int* const inSource, int* const inTag, int* const inFilledSize) = 0;

private:

	/** Forbiden (private) default constructor */
	FAbstractApplication():argc(0), argv(0){}

	/** Forbiden (private) copy constructor */
	FAbstractApplication(const FAbstractApplication&):argc(0), argv(0){}
	
	/** Forbiden (private) copy */
        FAbstractApplication& operator=(const FAbstractApplication&){return *this;}

public:
	/**
	* Constructor
	* @param inArgc argc from command line
	* @param inArgv argv from command line
        * This will also set the current app in the assert system
	*/
	FAbstractApplication(const int inArgc, char ** const inArgv )
		: argc(inArgc), argv(inArgv) {
		FAssertable::SetCurrentApp(this);
	}

	/** Destructor */
	virtual ~FAbstractApplication(){}

	/**
	* This function has to be called to execute the process run	
	* @return 0 if success
	*/
	int execute(){
		if( isMaster() ) initMaster();
		else initSlave();
		run();
		return 0;
	}

	/**
	* To get the current process id	
	* @return the process numeber [0 ; processCount [
	*/
	virtual int processId() const = 0;

	/**
	* To get the number of process
	* @return process count
	*/
	virtual int processCount() const = 0;

	/**
	* To make a barrier between process
	*/
	virtual void processBarrier() const = 0;

	/**
	* To kill all process
	* @param inErrorCode the error to return to OS (default is 1)
	*/
	virtual void abort(const int inErrorCode = 1) const = 0;	

	/**
	* This function has to be used to know if the current app is the master
	* @return true if id() == 0
	*/	
	bool isMaster() const {
		return !processId();
	}

	/**
	* This function has to be used to know if the current app is alone
	* @return true if processCount() == 1
	*/	
	bool isAlone() const {
		return processCount() == 1;
	}

	/**
	* This function has to be used to know if the current app is a slave
	* @return true if id() != 0
	*/
	bool isSlave() const {
		return processId();
	}

	/**
	* This function gives the number of parameters (argc)
	* @return argc
	*/
	int userParemetersCount() const{
		return this->argc;
	}

	/**
	* This function gives a parameter
	* @parameter inArg parameter position has to be strictly less than argc/userParemetersCount
	* @return argv[inArg]
	*/
	const char* userParemeterAt(const int inArg) const{
		return this->argv[inArg];
	}

	/**
	* This function gives a parameter in a standart type
	* @parameter inArg parameter position has to be strictly less than argc/userParemetersCount
	* @return argv[inArg] in the template VariableType form
        * @warning VariableType need to work with istream >> operator
        * <code> const int argInt = userParemetersAt<int>(1,-1); </code>
	*/
	template <class VariableType>
	const VariableType userParemeterAt(const int inArg, const VariableType& defaultValue = VariableType()) const{
		std::istringstream iss(this->argv[inArg],std::istringstream::in);
		VariableType value;
		iss >> value;
		if( /*iss.tellg()*/ iss.eof() ) return value;
		return defaultValue;
	}

	/**
	* This function gives the parameter in a standart type after a key parameter
	* Do not use pointer in template type!
	* @parameter inArg parameter key
	* @return argv[inArg.position + 1] in the template VariableType form
        * @warning VariableType need to work with istream >> operator
        * <code> const int argInt = userParemetersAt<int>(1,-1); </code>
	*/
	template <class VariableType>
	const VariableType userParemeterFromKey(const char* const inKey, const VariableType& defaultValue = VariableType(), bool* const inState = 0) const{
		const int keysArgc= this->argc - 1;
		// loop from 1 to argc  1
		for(int indexArg = 1 ; indexArg < keysArgc ; ++indexArg){
			// if argv == inArg
			if(strcmp(this->argv[indexArg] , inKey) == 0){
				// the argv + 1 => variable to use
				std::istringstream iss(this->argv[indexArg + 1],std::istringstream::in);
				VariableType value;
				iss >> value;
				// if we can cast to the template type
				if( iss.eof() ){
					if( inState ) *inState = true;
					return value;
				}
				break;
			}
		}
		// cannot cast to template or key not found
		if( inState ) *inState = false;
		return defaultValue;
	}

};

#endif //FABSTRACTAPPLICATION_HPP

// [--LICENSE--]
