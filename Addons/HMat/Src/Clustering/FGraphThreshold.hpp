#ifndef FGRAPHTHRESHOLD_HPP
#define FGRAPHTHRESHOLD_HPP

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

#include <scotch.h>

template <class FReal>
class FGraphThreshold {
protected:
    const int dim;
    const FReal treshold;

    int* permutations;
    std::pair<int, int>* gclusters;

public:
    FGraphThreshold(const int inDim, const FReal inDistMat[], const FReal inTreshold)
        : dim(inDim), treshold(inTreshold),
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

            int nbValuesUnderThreshold = 0;
            for(int idxCol = 0 ; idxCol < sizeInterval ; ++idxCol){
                const int idxColReal = permutations[idxCol+currentCluster.first];
                for(int idxRow = 0 ; idxRow < sizeInterval ; ++idxRow){
                    const int idxRowReal = permutations[idxRow+currentCluster.first];
                    nbValuesUnderThreshold += (inDistMat[idxColReal*dim+ idxRowReal] < treshold && idxRowReal != idxColReal? 1 : 0);
                }
            }
            FAssertLF(nbValuesUnderThreshold != 0);


            // Base value for all array indexings
            const SCOTCH_Num baseval = 0;
            // Number of vertices in graph
            const SCOTCH_Num vertnbr = sizeInterval;
            // Number of arcs in graph. Since edges are represented by both of their ends,
            // the number of edge data in the graph is twice the number of graph edges.
            const SCOTCH_Num edgenbr = nbValuesUnderThreshold;// it is already symmetric
            // Array of start indices in edgetab of vertex adjacency sub-arrays
            SCOTCH_Num* verttab = new SCOTCH_Num[vertnbr+1];
            SCOTCH_Num* vendtab = &verttab[1];
            // edgetab[verttab[i]] to edgetab[vendtab[i] − 1]
            SCOTCH_Num* edgetab = new SCOTCH_Num[edgenbr];
            // Optional array, of size vertnbr, holding the integer load associated with every vertex.
            SCOTCH_Num* velotab = nullptr;
            // Optional array, of a size equal at least to (max i (vendtab[i]) − baseval), holding the integer load associated with every arc.
            SCOTCH_Num* edlotab = nullptr;
            SCOTCH_Num* vlbltab = nullptr;

            verttab[0] = 0;
            for(int idxCol = 0 ; idxCol < sizeInterval ; ++idxCol){
                const int idxColReal = permutations[idxCol+currentCluster.first];
                verttab[idxCol+1] = verttab[idxCol];
                for(int idxRow = 0 ; idxRow < sizeInterval ; ++idxRow){
                    const int idxRowReal = permutations[idxRow+currentCluster.first];
                    if(inDistMat[idxColReal*dim+ idxRowReal] < treshold
                                && idxColReal != idxRowReal){
                        edgetab[verttab[idxCol+1]] = idxRow;
                        verttab[idxCol+1] += 1;
                    }
                }
            }
            FAssertLF(verttab[vertnbr] == edgenbr);

            SCOTCH_Graph grafdat;
            FAssertLF(SCOTCH_graphInit(&grafdat) == 0);

            FAssertLF(SCOTCH_graphBuild(&grafdat, baseval, vertnbr, verttab, vendtab,
                              velotab, vlbltab, edgenbr, edgetab, edlotab) == 0);
            FAssertLF(SCOTCH_graphCheck(&grafdat) == 0);

            SCOTCH_Strat straptr;
            FAssertLF(SCOTCH_stratInit(&straptr) == 0);


            //const SCOTCH_Num pwgtmax = std::numeric_limits<SCOTCH_Num>::max(); //maximum cluster vertex weight
            //const double densmin = 0; // the minimum edge density
            //const double bbalval = 0.2;// bipartition imbalance ratio
            //SCOTCH_stratGraphClusterBuild (&straptr, SCOTCH_STRATDEFAULT, pwgtmax, densmin, bbalval);

            const int partnbr   = 2;
            SCOTCH_Num* parttab = new SCOTCH_Num[vertnbr];
            FAssertLF(SCOTCH_graphPart(&grafdat, partnbr, &straptr, parttab) == 0);

            {
                int idxLeft  = 0;
                int idxRight = sizeInterval-1;
                while(idxLeft <= idxRight){
                    if(parttab[idxLeft] == 0){
                        idxLeft += 1;
                    }
                    else if(parttab[idxRight] == 1){
                        idxRight -= 1;
                    }
                    else{
                        const int unk     = parttab[idxLeft];
                        parttab[idxLeft]  = parttab[idxRight];
                        parttab[idxRight] = unk;

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
                fflush(stdout);
            }

            SCOTCH_stratExit(&straptr);
            SCOTCH_graphExit(&grafdat);

            delete[] parttab;
            delete[] verttab;
            delete[] edgetab;

            idxCluster -= 1;
        }
        FAssertLF(idxCluster == -1);
        FAssertLF(idxChildCluster == -1);
        FAssertLF(countToUnknowns == dim);
    }

    ~FGraphThreshold(){
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



#endif // FGRAPHTHRESHOLD_HPP

