// @SCALFMM_PRIVATE

#ifndef FCCLTREECLUSTER_HPP
#define FCCLTREECLUSTER_HPP

#include "./Utils/FGlobal.hpp"
#include "./Utils/FAssert.hpp"
#include "./Utils/FMath.hpp"
#include "../Utils/FHUtils.hpp"

#include "FClusterTree.hpp"

#include <stack>
#include <vector>
#include <functional>
#include <queue>
#include <limits>

extern "C" {
#include <cluster.h>
}


namespace CCL {
    enum TreeMethod {
        CCL_TM_SINGLE,
        CCL_TM_MAXIMUM,
        CCL_TM_AVG_LINKAGE
    };

    inline char TreeMethodToChar(const TreeMethod method){
        switch (method) {
        case CCL_TM_SINGLE:
            return 's';
            break;
        case CCL_TM_MAXIMUM:
            return 'm';
            break;
        case CCL_TM_AVG_LINKAGE:
            return 'a';
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
class FCCLTreeCluster {
protected:
    const int dim;
    CCL::TreeMethod method;
    Node* croot;

public:
    FCCLTreeCluster(const int inDim, const FReal inDistMat[], const CCL::TreeMethod inMethod)
        : dim(inDim), method(inMethod), croot(nullptr){

        double** distMatPtrs = new double*[dim];
        int** mask = new int*[dim];

        // Build mask, everyone is here
        for(int idxRow = 0 ; idxRow < dim ; ++idxRow){
            mask[idxRow]      = new int[idxRow+1];
            distMatPtrs[idxRow] = new double[idxRow+1];
            for(int idxCol = 0 ; idxCol <= idxRow ; ++idxCol){
                mask[idxRow][idxCol] = 1;
                distMatPtrs[idxRow][idxCol] = double(inDistMat[idxCol*dim + idxRow]);
            }
        }

        croot = treecluster (dim, dim, nullptr, mask, nullptr, 0, '?', TreeMethodToChar(method), distMatPtrs);
        FAssertLF(croot);

        for(int idxRow = 0 ; idxRow < dim ; ++idxRow){
            delete[] mask[idxRow];
            delete[] distMatPtrs[idxRow];
        }
        delete[] mask;
        delete[] distMatPtrs;
    }

    ~FCCLTreeCluster(){
        free(croot);
    }

    void fillClusterTree(FClusterTree<FReal>* ctree) const {
        int* permsOrigToNew = new int[dim];
        int* permsNewToOrig = new int[dim];

        {
            std::stack<int> depthFirst;
            depthFirst.push(croot[dim-2].right);
            depthFirst.push(croot[dim-2].left);

            int idxPerm = 0;
            while(depthFirst.size()){
                const int current = depthFirst.top();
                depthFirst.pop();
                if(0 <= current){
                    permsOrigToNew[current] = idxPerm;
                    permsNewToOrig[idxPerm] = current;
                    idxPerm += 1;
                }
                else{
                    depthFirst.push(croot[-1-current].right);
                    depthFirst.push(croot[-1-current].left);
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
            const Node& srcNd = croot[idxCluster];
            currentNd.id    = (-idxCluster)-1;
            currentNd.size  = 0;
            currentNd.left  = srcNd.left;
            if(0 <= srcNd.left){
                currentNd.size += 1;
                leaves[srcNd.left].parent = currentNd.id;
            }
            else{
                currentNd.size += clusters[(-srcNd.left)-1].size;
                clusters[(-srcNd.left)-1].parent = currentNd.id;
            }
            currentNd.right  = srcNd.right;
            if(0 <= srcNd.right){
                currentNd.size  += 1;
                leaves[srcNd.right].parent = currentNd.id;
            }
            else{
                currentNd.size  += clusters[(-srcNd.right)-1].size;
                clusters[(-srcNd.right)-1].parent = currentNd.id;
            }
            currentNd.score  = srcNd.distance;
         }
         clusters[dim-2].parent = 0;

         ctree->setData(dim, dim-1, clusters, dim, leaves, permsOrigToNew, permsNewToOrig);
         // We do not deallocate, ctree is in charge of this
    }


};

#endif // FCCLTREECLUSTER_HPP

