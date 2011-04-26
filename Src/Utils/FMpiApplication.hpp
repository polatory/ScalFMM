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


        void sendData(const int inReceiver, const int inSize, void* const inData, const int inTag){
            MPI_Request request;
            MPI_Isend(inData, inSize, MPI_CHAR , inReceiver, inTag, MPI_COMM_WORLD, &request);
        }

        void receiveData(const int inSize, void* const inData, int* const inSource, int* const inTag, int* const inFilledSize){
            MPI_Status status;
            MPI_Recv(inData, inSize, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, &status);
            *inSource = status.MPI_SOURCE;
            *inTag = status.MPI_TAG;
            MPI_Get_count(&status,MPI_CHAR,inFilledSize);
        }

        /**
          * To know if some data arrived
          */
        virtual bool receivedData(){
            int flag;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
            return flag;
        }

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
