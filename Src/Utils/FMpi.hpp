#ifndef FMPI_HPP
#define FMPI_HPP
// /!\ Please, you must read the license at the bottom of this page


#include "FGlobal.hpp"

#ifdef SCALFMM_USE_MPI
#include <mpi.h>
#endif

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FMpi
* Please read the license
*
* This namespace is to compile with or without mpi
*/

class FMpi {
public:
#ifdef SCALFMM_USE_MPI

////////////////////////////////////////////////////////
// Use MPI
////////////////////////////////////////////////////////

    FMpi(int inArgc, char **  inArgv ) {
        MPI_Init(&inArgc,&inArgv);
    }

    ~FMpi(){
        MPI_Finalize();
    }

    void sendData(const int inReceiver, const int inSize, void* const inData, const int inTag){
        //MPI_Request request;
        //MPI_Isend(inData, inSize, MPI_CHAR , inReceiver, inTag, MPI_COMM_WORLD, &request);
        MPI_Send(inData, inSize, MPI_CHAR , inReceiver, inTag, MPI_COMM_WORLD);
    }

    void receiveData(const int inSize, void* const inData, int* const inSource, int* const inTag, int* const inFilledSize){
        MPI_Status status;
        MPI_Recv(inData, inSize, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, &status);
        *inSource = status.MPI_SOURCE;
        *inTag = status.MPI_TAG;
        MPI_Get_count(&status,MPI_CHAR,inFilledSize);
    }

    void receiveData(const int inSize, const int inTag, void* const inData, int* const inSource, int* const inFilledSize){
        MPI_Status status;
        MPI_Recv(inData, inSize, MPI_CHAR, MPI_ANY_SOURCE, inTag, MPI_COMM_WORLD, &status);
        *inSource = status.MPI_SOURCE;
        MPI_Get_count(&status,MPI_CHAR,inFilledSize);
    }

    bool receivedData(){
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
        return flag;
    }

    bool receivedData(int* const tag){
        MPI_Status status;
        int flag;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        *tag = status.MPI_TAG;
        return flag;
    }

    int processId() {
        int id;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        return id;
    }

    int processCount() {
        int count;
        MPI_Comm_size(MPI_COMM_WORLD,&count);
        return count;
    }

    void processBarrier() {
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void abort(const int inErrorCode = 1) {
        MPI_Abort(MPI_COMM_WORLD, inErrorCode);
    }


    double reduceSum(double data){
        double result;
        MPI_Reduce( &data, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
        return result;
    }

#else
////////////////////////////////////////////////////////
// Without MPI
////////////////////////////////////////////////////////

    FMpi(int inArgc, char **  inArgv ) {}

    void sendData(const int, const int, void* const, const int ){}

    void receiveData(const int, void* const, int* const inSource, int* const inTag, int* const inFilledSize){
        *inSource = 0;
        *inTag = 0;
        *inFilledSize = 0;
    }

    void receiveData(const int , const int , void* const , int* const inSource, int* const inFilledSize){
        *inSource = 0;
        *inFilledSize = 0;
    }

    bool receivedData(){
        return false;
    }

    bool receivedData(int* const){
        return false;
    }

    int processId() {
        return 0;
    }

    int processCount() {
        return 1;
    }

    void processBarrier() {}

    void abort(const int inErrorCode = 1) {
        exit(inErrorCode);
    }

    double reduceSum(double data){
        return data;
    }

#endif

////////////////////////////////////////////////////////
// To use in any case
////////////////////////////////////////////////////////
    bool isMaster() {
        return !processId();
    }

    bool isAlone() {
        return processCount() == 1;
    }

    bool isSlave() {
        return processId();
    }
};

#endif //FMPI_HPP

// [--LICENSE--]
