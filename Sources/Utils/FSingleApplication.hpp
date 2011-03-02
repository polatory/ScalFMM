#ifndef FSINGLEAPPLICATION_HPP
#define FSINGLEAPPLICATION_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <stdlib.h>


#include "FAbstractApplication.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FSingleApplication
* Please read the license
*
* This class is an application of abstract application with as a standalone process
*
* @warning you have to implement run() and call execute to start the app
*
* Please refere to testApplication.cpp to see an example
* <code>
* </code>
*/
class FSingleApplication : public FAbstractApplication {
protected:
	/**
	* This will be called as the main method
        * @warning Must be impleted in the derived class
	*/
	virtual void run() = 0;


public:
	/**
	* Constructor
	* @param inArgc argc from command line
	* @param inArgv argv from command line
	*/
	FSingleApplication(const int inArgc, char ** const inArgv )
		: FAbstractApplication(inArgc,inArgv) {
	}

	/** Destructor */
	virtual ~FSingleApplication(){}

	/**
	* To get the current process id	
	* @return the process numeber [0 ; processCount [
	*/
	int processId() const {
		return 0;
	}

	/**
	* To get the number of process
	* @return process count
	*/
	int processCount() const {
		return 1;
	}

	/**
	* To make a barrier between process
	*/
	void processBarrier() const{}

	/**
	* To kill all process
	* @param inErrorCode the error to return to OS (default is 1)
	*/
	void abort(const int inErrorCode = 1) const {
		exit(inErrorCode);
	}	

};

#endif //FSINGLEAPPLICATION_HPP

// [--LICENSE--]
