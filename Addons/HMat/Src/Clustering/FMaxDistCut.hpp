#ifndef FMAXDISTCUT_HPP
#define FMAXDISTCUT_HPP

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
class FMaxDistCut {
protected:
    const int dim;

    int* permutations;
    std::pair<int, int>* gclusters;

public:
    FMaxDistCut(const int inDim, const FReal inDistMat[])
        : dim(inDim),
          permutations(new int[dim]),
          gclusters (new std::pair<int, int>[dim-1]){

        for(int idx = 0 ; idx < dim ; ++idx){
            permutations[idx] = idx;
        }

        std::queue<std::pair<int, int>> intervals;
        intervals.push(std::pair<int,int>(0, dim));

        int idxCluster = dim-2;
        int idxChildCluster = dim-3;
        int countToUnknowns = 0;

        while(intervals.size()){
            FAssertLF(idxCluster >= 0);

            const std::pair<int,int> currentCluster = intervals.front();
            intervals.pop();

            const int sizeInterval = currentCluster.second-currentCluster.first;
            FAssertLF(sizeInterval != 1);

            int firstIdxMax  = permutations[0 + currentCluster.first];
            int secondIdxMax = permutations[1 + currentCluster.first];
            for(int idxCol = 0 ; idxCol < sizeInterval ; ++idxCol){
                const int idxColReal = permutations[idxCol+currentCluster.first];
                for(int idxRow = idxCol+1 ; idxRow < sizeInterval ; ++idxRow){
                    const int idxRowReal = permutations[idxRow+currentCluster.first];
                    if(inDistMat[idxColReal*dim+ idxRowReal]
                            > inDistMat[firstIdxMax*dim + secondIdxMax]){
                        firstIdxMax  = idxColReal;
                        secondIdxMax = idxRowReal;
                    }
                }
            }

            {
                int idxLeft  = 0;
                int idxRight = sizeInterval-1;
                while(idxLeft <= idxRight){
                    const int idxLeftReal  = permutations[idxLeft + currentCluster.first];
                    const int idxRightReal = permutations[idxRight + currentCluster.first];
                    if(inDistMat[idxLeftReal*dim+ firstIdxMax]
                            < inDistMat[idxLeftReal*dim + secondIdxMax]){
                        idxLeft += 1;
                    }
                    else if(inDistMat[idxRightReal*dim+ firstIdxMax]
                            > inDistMat[idxRightReal*dim + secondIdxMax]){
                        idxRight -= 1;
                    }
                    else{
                        const int unkPerm = permutations[idxLeft+currentCluster.first];
                        permutations[idxLeft+currentCluster.first]  = permutations[idxRight+currentCluster.first];
                        permutations[idxRight+currentCluster.first] = unkPerm;

                        idxLeft  += 1;
                        idxRight -= 1;
                    }
                }
                // idxLeft is on the first 1
                if(idxLeft == 1){
                    gclusters[idxCluster].first = permutations[currentCluster.first];
                    countToUnknowns += 1;
                    FAssertLF(countToUnknowns <= dim);
                }
                else if(idxLeft > 1){
                    FAssertLF(idxChildCluster >= 0);
                    gclusters[idxCluster].first = (-idxChildCluster)-1;
                    idxChildCluster -= 1;
                    intervals.push(std::pair<int,int>(currentCluster.first, currentCluster.first+idxLeft));
                }
                if(idxRight == sizeInterval-2){
                    gclusters[idxCluster].second = permutations[sizeInterval-1+currentCluster.first];
                    countToUnknowns += 1;
                    FAssertLF(countToUnknowns <= dim);
                }
                else if(idxRight < sizeInterval-2){
                    FAssertLF(idxChildCluster >= 0);
                    gclusters[idxCluster].second = (-idxChildCluster)-1;
                    idxChildCluster -= 1;
                    intervals.push(std::pair<int,int>(currentCluster.first+idxLeft, currentCluster.first+sizeInterval));
                }
            }

            idxCluster -= 1;
        }
        FAssertLF(idxCluster == -1);
        FAssertLF(idxChildCluster == -1);
        FAssertLF(countToUnknowns == dim);
    }

    ~FMaxDistCut(){
        delete[] permutations;
        delete[] gclusters;
    }

    void fillClusterTree(FClusterTree<FReal>* ctree){
        int* permsOrigToNew = new int[dim];
        int* permsNewToOrig = new int[dim];

        {
            std::stack<int> depthFirst;
            depthFirst.push(gclusters[dim-2].second);
            depthFirst.push(gclusters[dim-2].first);

            int idxPerm = 0;
            while(depthFirst.size()){
                const int current = depthFirst.top();
                depthFirst.pop();
                if(0 <= current){
                    permsOrigToNew[current] = idxPerm;
                    permsNewToOrig[idxPerm] = current;
                    FAssertLF(permsNewToOrig[idxPerm] == permutations[idxPerm]);
                    idxPerm += 1;
                }
                else{
                    depthFirst.push(gclusters[-1-current].second);
                    depthFirst.push(gclusters[-1-current].first);
                }
            }
        }

        typename FClusterTree<FReal>::Leaf* leaves = new typename FClusterTree<FReal>::Leaf[dim];
         {
             for(int idxUnk = 0 ; idxUnk < dim ; ++idxUnk){
                leaves[idxUnk].id = idxUnk;
                leaves[idxUnk].offset= idxUnk;
                leaves[idxUnk].size  = 1;
             }
         }

        typename FClusterTree<FReal>::Node* clusters = new typename FClusterTree<FReal>::Node[dim-1];

         for(int idxCluster = 0 ; idxCluster < dim-1 ; ++idxCluster){
            typename FClusterTree<FReal>::Node& currentNd = clusters[idxCluster];
            const std::pair<int,int>& srcNd = gclusters[idxCluster];
            currentNd.id    = (-idxCluster)-1;
            currentNd.size  = 0;
            currentNd.left  = srcNd.first;
            if(0 <= srcNd.first){
                currentNd.size += 1;
                leaves[srcNd.first].parent = currentNd.id;
            }
            else{
                currentNd.size += clusters[(-srcNd.first)-1].size;
                clusters[(-srcNd.first)-1].parent = currentNd.id;
            }
            currentNd.right  = srcNd.second;
            if(0 <= srcNd.second){
                currentNd.size  += 1;
                leaves[srcNd.second].parent = currentNd.id;
            }
            else{
                currentNd.size  += clusters[(-srcNd.second)-1].size;
                clusters[(-srcNd.second)-1].parent = currentNd.id;
            }
            currentNd.score  = 0;
         }
         clusters[dim-2].parent = 0;

         ctree->setData(dim, dim-1, clusters, dim, leaves, permsOrigToNew, permsNewToOrig);
         // We do not deallocate, ctree is in charge of this
    }


};


#endif // FMAXDISTCUT_HPP

