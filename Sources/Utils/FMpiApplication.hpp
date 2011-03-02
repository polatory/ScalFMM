#ifndef FMPIAPPLICATION_HPP
#define FMPIAPPLICATION_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <mpi.h>

#include "FAbstractApplication.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FMpiApplication
* Please read the license
*
* This class is an implementation of the abstract application with mpi
*
* @warning you have to implement run() and call execute to start the app
*
* Please refere to testApplication.cpp to see an example
* <code>
* </code>
*/
class FMpiApplication : public FAbstractApplication {
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
	FMpiApplication(int inArgc, char **  inArgv )
		: FAbstractApplication(inArgc,inArgv) {
		MPI_Init(&inArgc,&inArgv);
	}

	/** Destructor */
	virtual ~FMpiApplication(){
		MPI_Finalize();
	}

	/**
	* To get the current process id	
	* @return the process numeber [0 ; processCount [
	*/
	int processId() const {
		int id;
    		MPI_Comm_rank(MPI_COMM_WORLD,&id); 
		return id;
	}

	/**
	* To get the number of process
	* @return process count
	*/
	int processCount() const {
		int count;
		MPI_Comm_size(MPI_COMM_WORLD,&count);
		return count;
	}

	/**
	* To make a barrier between process
	*/
	void processBarrier() const{
		MPI_Barrier(MPI_COMM_WORLD);
	}

	/**
	* To kill all process
	* @param inErrorCode the error to return to OS (default is 1)
	*/
	void abort(const int inErrorCode = 1) const {
		MPI_Abort(MPI_COMM_WORLD, inErrorCode);
	}	

	

};

#endif //FMPIAPPLICATION_HPP

// [--LICENSE--]
