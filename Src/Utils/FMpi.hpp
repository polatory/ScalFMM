#ifndef FMPI_HPP
#define FMPI_HPP
// /!\ Please, you must read the license at the bottom of this page


#include "FGlobal.hpp"
#include "FMath.hpp"

#ifdef SCALFMM_USE_MPI
#include <mpi.h>
#endif

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FMpi
* Please read the license
*
* This namespace is to compile with or without mpi
* It defines a class to access MPI data, if the lib is compiled
* without mpi support then simulate data.
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

    void receiveData(const int inSize, void* const inData, int* const inSource = 0, int* const inTag = 0, int* const inFilledSize = 0){
        MPI_Status status;
        MPI_Recv(inData, inSize, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, &status);
        if(inSource) *inSource = status.MPI_SOURCE;
        if(inTag) *inTag = status.MPI_TAG;
        if(inFilledSize) MPI_Get_count(&status,MPI_CHAR,inFilledSize);
    }

    void receiveDataFromTag(const int inSize, const int inTag, void* const inData, int* const inSource = 0, int* const inFilledSize = 0){
        MPI_Status status;
        MPI_Recv(inData, inSize, MPI_CHAR, MPI_ANY_SOURCE, inTag, MPI_COMM_WORLD, &status);
        if(inSource) *inSource = status.MPI_SOURCE;
        if(inFilledSize) MPI_Get_count(&status,MPI_CHAR,inFilledSize);
    }

    void receiveDataFromTagAndSource(const int inSize, const int inTag, const int inSource, void* const inData, int* const inFilledSize = 0){
        MPI_Status status;
        MPI_Recv(inData, inSize, MPI_CHAR, inSource, inTag, MPI_COMM_WORLD, &status);
        if(inFilledSize) MPI_Get_count(&status,MPI_CHAR,inFilledSize);
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

    template < class T >
    MPI_Datatype getType(){
        return MPI_INT;
    }

    template< class T >
    T reduceSum(T data){
        T result;
        MPI_Reduce( &data, &result, 1, getType<T>(), MPI_SUM, 0, MPI_COMM_WORLD );
        return result;
    }

    template< class T >
    T reduceMin(T myMin){
        T result;
        MPI_Reduce( &myMin, &result, 1, getType<T>(), MPI_MIN, 0, MPI_COMM_WORLD );
        return result;
    }

    template< class T >
    T reduceMax(T myMax){
        T result;
        MPI_Reduce( &myMax, &result, 1, getType<T>(), MPI_MAX, 0, MPI_COMM_WORLD );
        return result;
    }

    template< class T >
    T broadcast(T value, const int from = 0){
        MPI_Bcast ( &value, 1, getType<T>(), from, MPI_COMM_WORLD );
        return value;
    }

#else
////////////////////////////////////////////////////////
// Without MPI
////////////////////////////////////////////////////////

    FMpi(int inArgc, char **  inArgv ) {}

    void sendData(const int, const int, void* const, const int ){}

    void receiveData(const int, void* const, int* const inSource = 0, int* const inTag = 0, int* const inFilledSize = 0){
        if(inSource) *inSource = 0;
        if(inTag) *inTag = 0;
        if(inFilledSize) *inFilledSize = 0;
    }

    void receiveDataFromTag(const int , const int , void* const , int* const inSource = 0, int* const inFilledSize = 0){
        if(inSource) *inSource = 0;
        if(inFilledSize) *inFilledSize = 0;
    }

    void receiveDataFromTagAndSource(const int , const int , const int , void* const , int* const inFilledSize = 0){
        if(inFilledSize) *inFilledSize = 0;
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

    template< class T >
    T reduceSum(T data){
        return data;
    }

    template< class T >
    T reduceMin(T myMin){
        return myMin;
    }

    template< class T >
    T reduceMax(T myMax){
        return myMax;
    }

    template< class T >
    T broadcast(T value, const int = 0){
        return value;
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

    template< class T >
    T getLeft(const T inSize) {
        const float step = (float(inSize) / processCount());
        return T(FMath::Ceil(step * processId()));
    }

    template< class T >
    T getRight(const T inSize) {
        const float step = (float(inSize) / processCount());
        const T res = T(FMath::Ceil(step * (processId()+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    template< class T >
    T getOtherRight(const T inSize, const int other) {
        const float step = (float(inSize) / processCount());
        const T res = T(FMath::Ceil(step * (other+1)));
        if(res > inSize) return inSize;
        else return res;
    }

    template< class T >
    int getProc(const int position, const T inSize) {
        const float step = (float(inSize) / processCount());
        return int(position/step);
    }
};

#ifdef SCALFMM_USE_MPI
template <>
MPI_Datatype FMpi::getType<long long>(){
    return MPI_LONG_LONG;
}

template <>
MPI_Datatype FMpi::getType<double>(){
    return MPI_DOUBLE;
}

template <>
MPI_Datatype FMpi::getType<int>(){
    return MPI_INT;
}
#endif

#endif //FMPI_HPP

// [--LICENSE--]
