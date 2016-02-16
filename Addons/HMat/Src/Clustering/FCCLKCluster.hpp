// @SCALFMM_PRIVATE

#ifndef FCCLKCLUSTER_HPP
#define FCCLKCLUSTER_HPP

#include "./Utils/FGlobal.hpp"
#include "./Utils/FAssert.hpp"
#include "./Utils/FMath.hpp"
#include "../Utils/FHUtils.hpp"

#include <stack>
#include <vector>
#include <functional>
#include <queue>
#include <limits>
#include <algorithm>

extern "C" {
#include <cluster.h>
}


namespace CCL {
    enum ClusterCenterMethod {
        CCL_CCM_ARITHMETIC_MEAN,
        CCL_CCM_MEDIAN,
        CCL_CCM_DUMMY
    };

    inline char ClusterCenterMethodToChar(const ClusterCenterMethod method){
        switch (method) {
        case CCL_CCM_ARITHMETIC_MEAN:
            return 'a';
            break;
        case CCL_CCM_MEDIAN:
            return 'm';
            break;
        default:
            break;
        }
        return '?';
    }

    enum Distance {
        CCL_DIST_MEAN,
        CCL_DIST_MEDIAN,
        CCL_DIST_SHORTEST,
        CCL_DIST_LONGEST,
        CCL_DIST_AVG,
        CCL_DIST_DUMMY
    };

    inline char DistanceToChar(const Distance method){
        switch (method) {
        case CCL_DIST_MEAN:
            return 'a';
            break;
        case CCL_DIST_MEDIAN:
            return 'm';
            break;
        case CCL_DIST_SHORTEST:
            return 's';
            break;
        case CCL_DIST_LONGEST:
            return 'x';
            break;
        case CCL_DIST_AVG:
            return 'v';
            break;
        default:
            break;
            }
        return '?';
    }
}


template <class FReal>
class FCCLKCluster {
protected:
    const int nbPartitions;
    const int nbElements;
    const int nbDim;
    const int nbPass; //< Number of call to EM algorithm
    CCL::ClusterCenterMethod method;
    CCL::Distance distance;
    int* partitions;

public:
    FCCLKCluster(const int inNbPartitions, const int inNbElements, const FReal inDistMat[], const int inNbPass = 0)
        : nbPartitions(inNbPartitions), nbElements(inNbElements), nbDim(0), nbPass(inNbPass), method(CCL::CCL_CCM_DUMMY), distance(CCL::CCL_DIST_DUMMY), partitions(nullptr) {

        double** distMatPtrs = new double*[nbElements];

        // Build mask, everyone is here
        for(int idxRow = 0 ; idxRow < nbElements ; ++idxRow){
            distMatPtrs[idxRow] = new double[idxRow+1];
            for(int idxCol = 0 ; idxCol <= idxRow ; ++idxCol){
                distMatPtrs[idxRow][idxCol] = double(inDistMat[idxCol*nbElements + idxRow]);
            }
        }

        // allocate partitions
        partitions = new int[nbElements];

        // Errors
        double* error = new double[nbElements];
        // Nb of times the optimal clustering was found
        int* ifound = new int[nbElements];

        kmedoids (nbPartitions, nbElements, distMatPtrs, nbPass, partitions, error, ifound);

        for(int idxRow = 0 ; idxRow < nbElements ; ++idxRow){
            delete[] distMatPtrs[idxRow];
        }
        delete[] distMatPtrs;

    }

    FCCLKCluster(const int inNbPartitions, const int inNbElements, const int inNbDim, const FReal inDataMat[], const CCL::ClusterCenterMethod inMethod, const CCL::Distance inDistance, const int inNbPass = 0)
        : nbPartitions(inNbPartitions), nbElements(inNbElements), nbDim(inNbDim), nbPass(inNbPass), method(inMethod), distance(inDistance), partitions(nullptr) {

        double** dataMatPtrs = new double*[nbElements];
        int** mask = new int*[nbElements];

        // Build mask, everyone is here
        for(int idxRow = 0 ; idxRow < nbElements ; ++idxRow){
            mask[idxRow]      = new int[idxRow+1];
            dataMatPtrs[idxRow] = new double[nbDim];
            for(int idxCol = 0 ; idxCol < nbDim ; ++idxCol){
                mask[idxRow][idxCol] = 1;
                dataMatPtrs[idxRow][idxCol] = double(inDataMat[idxCol*nbElements + idxRow]);
            }
        }

        // allocate partitions
        partitions = new int[nbElements];

        // Errors
        double* error = new double[nbElements];
        // Nb of times the optimal clustering was found
        int* ifound = new int[nbElements];

        // Weights
        double* weights = new double[nbElements];
        for(int idxRow = 0 ; idxRow < nbElements ; ++idxRow)
            weights[idxRow]=double(1.0);

        kcluster(nbPartitions, nbElements, nbDim, dataMatPtrs, mask, weights, 0, nbPass, ClusterCenterMethodToChar(method), DistanceToChar(distance), partitions, error, ifound);

        for(int idxRow = 0 ; idxRow < nbElements ; ++idxRow){
            delete[] mask[idxRow];
            delete[] dataMatPtrs[idxRow];
        }
        delete[] mask;
        delete[] dataMatPtrs;

    }

    ~FCCLKCluster(){
        delete[] partitions;
    }

    int getPartitions(const int inNbPartitions, int inNbIdxInPartitions[]) const{

        /// Map partitions to 0.. nbPartitions
        FAssertLF(inNbPartitions == nbPartitions);

        // Copy partitions
        int* sortedPartitions = new int[nbElements];
        for(int idx = 0 ; idx < nbElements ; ++idx){
            sortedPartitions[idx]=partitions[idx];
        }
        // sort partitions
        std::sort(sortedPartitions,sortedPartitions+nbElements);

        // Map partitions to 0..nbPartitions
        int counterPartition=0;
        int* mapPartitions = new int[inNbPartitions];
        int currentPartition=sortedPartitions[0];
        for(int idx = 0 ; idx < nbElements ; ++idx){
            mapPartitions[counterPartition]=currentPartition;
            if(sortedPartitions[idx+1]!=currentPartition){
                currentPartition=sortedPartitions[idx+1];
                ++counterPartition;
            }
        }
        FAssertLF(counterPartition == inNbPartitions);

        /// Count particles in each partition
        int totalGiven = 0;
        for(int idxPartition = 0 ; idxPartition < inNbPartitions ; ++idxPartition){

            inNbIdxInPartitions[idxPartition]=0;

            for(int idx = 0 ; idx < nbElements ; ++idx){

                if(partitions[idx]==mapPartitions[idxPartition])
                    inNbIdxInPartitions[idxPartition]+=1;
            
            }
            totalGiven +=inNbIdxInPartitions[idxPartition];
        }
        FAssertLF(totalGiven == nbElements);

        return 0; // no empty partition in kclusters/kmedoids algorithms
    }


};

#endif // FCCLKCLUSTER_HPP

