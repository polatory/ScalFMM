#ifndef FCONNEXCLUSTERING_HPP
#define FCONNEXCLUSTERING_HPP

// @SCALFMM_PRIVATE

#include "./Utils/FGlobal.hpp"
#include "./Utils/FAssert.hpp"
#include "./Utils/FMath.hpp"
#include "./Containers/FBoolArray.hpp"

#include "FClusterTree.hpp"

#include <stack>
#include <vector>
#include <functional>
#include <queue>
#include <limits>
#include <memory>



template <class FReal>
class FConnexClustering {
protected:
    const int dim;

    int* permsNewToOrig;
    int* permsOrigToNew;
    int* partitions;
    int* partitionsOffset;
    int nbPartitions;

public:
    FConnexClustering(const int inDim, const FReal inDistMat[], const FReal thresh)
        : dim(inDim),
          permsNewToOrig(new int[dim]),
          permsOrigToNew(new int[dim]),
          partitions(new int[dim]),
          partitionsOffset(new int[dim+1]),
          nbPartitions (0){

        std::unique_ptr<int[]> partitionsMappin(new int[dim]);

        for(int idx = 0 ; idx < dim ; ++idx){
            partitionsMappin[idx] = -1;
        }

        partitionsOffset[0] = 0;
        partitions[0] = 0;

        for(int idx = 0 ; idx < dim ; ++idx){
            if(partitionsMappin[idx] == -1){
                FAssertLF(nbPartitions < dim);

                partitionsOffset[nbPartitions+1] = partitionsOffset[nbPartitions]+1;
                partitionsMappin[idx] = nbPartitions;

                int idxPartitionElement = partitionsOffset[nbPartitions];
                permsNewToOrig[idxPartitionElement] = idx;

                while(idxPartitionElement < partitionsOffset[nbPartitions+1]){
                    FAssertLF(idxPartitionElement < dim);

                    for(int idxOther = 0 ; idxOther < dim ; ++idxOther){
                        if(partitionsMappin[idxOther] == -1
                               && inDistMat[permsNewToOrig[idxPartitionElement]*dim + idxOther] < thresh){
                            partitionsMappin[idxOther] = nbPartitions;
                            permsNewToOrig[partitionsOffset[nbPartitions+1]] = idxOther;
                            permsOrigToNew[idxOther] = partitionsOffset[nbPartitions+1];
                            partitionsOffset[nbPartitions+1] += 1;
                            FAssertLF(partitionsOffset[nbPartitions+1] <= dim);
                        }
                    }

                    idxPartitionElement += 1;
                }

                partitions[nbPartitions] = partitionsOffset[nbPartitions+1]-partitionsOffset[nbPartitions];
                nbPartitions += 1;
            }
        }

        FAssertLF(partitionsOffset[nbPartitions] == dim);
    }

    ~FConnexClustering(){
        delete[] permsNewToOrig;
        delete[] permsOrigToNew;
        delete[] partitionsOffset;
        delete[] partitions;
    }


    void fillPermutations(int* inPermuts, int* invPermuts = nullptr) const {
        memcpy(inPermuts, permsOrigToNew, sizeof(int)*dim);
        if(invPermuts){
            memcpy(invPermuts, permsNewToOrig, sizeof(int)*dim);
        }
    }

    int getNbPartitions() const{
        return nbPartitions;
    }

    void getPartitions(const int inNbPartitions, int inNbIdxInPartitions[]) const{
        FAssertLF(nbPartitions  == inNbPartitions);
        memcpy(inNbIdxInPartitions, partitions, sizeof(int)*nbPartitions);
    }
};

#endif // FCONNEXCLUSTERING_HPP

