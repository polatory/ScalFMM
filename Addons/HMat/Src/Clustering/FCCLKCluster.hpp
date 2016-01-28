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
        CCL_CCM_MEDIAN
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
        CCL_DIST_AVG
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
    const int dim;
    CCL::ClusterCenterMethod method;
    int* partitions;

public:
    FCCLKCluster(const int inNbPartitions, const int inDim, const FReal inDistMat[], const CCL::ClusterCenterMethod inMethod)
        : nbPartitions(inNbPartitions), dim(inDim), method(inMethod), partitions(nullptr) {

        double** distMatPtrs = new double*[dim];

        // Build mask, everyone is here
        for(int idxRow = 0 ; idxRow < dim ; ++idxRow){
            distMatPtrs[idxRow] = new double[idxRow+1];
            for(int idxCol = 0 ; idxCol <= idxRow ; ++idxCol){
                distMatPtrs[idxRow][idxCol] = double(inDistMat[idxCol*dim + idxRow]);
            }
        }

        // allocate partitions
        partitions = new int[dim];

        // Number of call to EM algorithm
        const int nbPass = 1;
        // Errors
        double* error = new double[dim];
        // Nb of times the optimal clustering was found
        int* ifound = new int[dim];

        kmedoids (nbPartitions, dim, distMatPtrs, nbPass, partitions, error, ifound);

        for(int idxRow = 0 ; idxRow < dim ; ++idxRow){
            delete[] distMatPtrs[idxRow];
        }
        delete[] distMatPtrs;
    }

    ~FCCLKCluster(){
        delete[] partitions;
    }

    int getPartitions(const int inNbPartitions, int inNbIdxInPartitions[]) const{

        /// Map partitions to 0.. nbPartitions
        FAssertLF(inNbPartitions == nbPartitions);

        // Copy partitions
        int* sortedPartitions = new int[dim];
        for(int idx = 0 ; idx < dim ; ++idx){
            sortedPartitions[idx]=partitions[idx];
        }
        // sort partitions
        std::sort(sortedPartitions,sortedPartitions+dim);

        // Map partitions to 0..nbPartitions
        int counterPartition=0;
        int* mapPartitions = new int[inNbPartitions];
        int currentPartition=sortedPartitions[0];
        for(int idx = 0 ; idx < dim ; ++idx){
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

            for(int idx = 0 ; idx < dim ; ++idx){

                if(partitions[idx]==mapPartitions[idxPartition])
                    inNbIdxInPartitions[idxPartition]+=1;
            
            }
            totalGiven +=inNbIdxInPartitions[idxPartition];
        }
        FAssertLF(totalGiven == dim);

        return 0; // no empty partition in kclusters/kmedoids algorithms
    }


};

#endif // FCCLKCLUSTER_HPP

