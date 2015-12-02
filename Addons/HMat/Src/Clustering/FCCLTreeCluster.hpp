// @SCALFMM_PRIVATE

#ifndef FCCLTREECLUSTER_HPP
#define FCCLTREECLUSTER_HPP

#include "./Utils/FGlobal.hpp"
#include "./Utils/FAssert.hpp"
#include "./Utils/FMath.hpp"
#include "../Utils/FHUtils.hpp"

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

    void fillPermutations(int* inPermuts, int* invPermuts = nullptr) const {
        std::stack<int> depthFirst;
        depthFirst.push(croot[dim-2].right);
        depthFirst.push(croot[dim-2].left);

        int idxPerm = 0;
        while(depthFirst.size()){
            const int current = depthFirst.top();
            depthFirst.pop();
            if(0 <= current){
                inPermuts[current] = idxPerm;
                if(invPermuts) invPermuts[idxPerm] = current;
                idxPerm += 1;
            }
            else{
                depthFirst.push(croot[-1-current].right);
                depthFirst.push(croot[-1-current].left);
            }
        }
    }

    void originalIdxToCluster(const int inNbClusters, int idxToClusters[]) const{
        cuttree(dim, croot, inNbClusters, idxToClusters);
    }

    void getPartitions(const int inHeight, const int inNbPartitions, int inNbIdxInPartitions[]) const{
        // Here we select the inNbPartitions last partitions
        // But we have to take them in the right order (left to right)
        // To ensure coherency with the permutations
        // if inNbPartitions, we should have inNbIdxInPartitions filled with 1
        FAssertLF(FMath::pow2(inHeight-1)  == inNbPartitions);

        std::vector<int> sizeOfClusters(dim-1, 0);
        for(int idxCluster = 0 ; idxCluster < dim-1 ; ++idxCluster){
            if(croot[idxCluster].left < 0){
                sizeOfClusters[idxCluster] += sizeOfClusters[-1-croot[idxCluster].left] ;
            }
            else{
                sizeOfClusters[idxCluster] += 1;
            }
            if(croot[idxCluster].right < 0){
                sizeOfClusters[idxCluster] += sizeOfClusters[-1-croot[idxCluster].right] ;
            }
            else{
                sizeOfClusters[idxCluster] += 1;
            }
        }
        FAssertLF(sizeOfClusters[dim-2] == dim);

        const int noNodeFlag = std::numeric_limits<int>::max();
        std::queue<int> breadthFirst;
        breadthFirst.push(croot[dim-2].left);
        breadthFirst.push(croot[dim-2].right);

        for(int idxLevel = 1 ; idxLevel < inHeight-1 ; ++idxLevel){
            std::queue<int> breadthFirstLower;

            while(breadthFirst.size()){
                const int current = breadthFirst.front();
                breadthFirst.pop();

                if(current == noNodeFlag){
                    breadthFirstLower.push(noNodeFlag);
                    breadthFirstLower.push(noNodeFlag);
                }
                else if(0 <= current){
                    breadthFirstLower.push(current);
                    breadthFirstLower.push(noNodeFlag);
                }
                else{
                    breadthFirstLower.push(croot[-1-current].left);
                    breadthFirstLower.push(croot[-1-current].right);
                }
            }

            breadthFirst = std::move(breadthFirstLower);
        }
        FAssertLF(int(breadthFirst.size()) == inNbPartitions);

        int counterPartition = 0;
        while(breadthFirst.size()){
            const int current = breadthFirst.front();
            breadthFirst.pop();

            if(current == noNodeFlag){
                inNbIdxInPartitions[counterPartition] = 0;
            }
            else if(0 <= current){
                inNbIdxInPartitions[counterPartition] = 1;
            }
            else{
                inNbIdxInPartitions[counterPartition] = sizeOfClusters[-1-current];
            }
            counterPartition += 1;
        }
    }


    void saveToXml(const char inDirname[], const char inFilename[]) const{
        char buffer[1024];
        sprintf(buffer, "%s/%s", inDirname, inFilename);
        saveToXml(buffer);
    }

    void saveToXml(const char inFilename[]) const{
        std::vector<int> sizeOfClusters(dim-1, 0);
        for(int idxCluster = 0 ; idxCluster < dim-1 ; ++idxCluster){
            if(croot[idxCluster].left < 0){
                sizeOfClusters[idxCluster] += sizeOfClusters[-1-croot[idxCluster].left] ;
            }
            else{
                sizeOfClusters[idxCluster] += 1;
            }
            if(croot[idxCluster].right < 0){
                sizeOfClusters[idxCluster] += sizeOfClusters[-1-croot[idxCluster].right] ;
            }
            else{
                sizeOfClusters[idxCluster] += 1;
            }
        }
        FAssertLF(sizeOfClusters[dim-2] == dim);


        FILE* fxml = fopen(inFilename, "w");
        FAssertLF(fxml);

        fprintf(fxml,"<?xml version=\"1.0\"?>\n");
        fprintf(fxml,"<data>\n");

        for(int idxCluster = dim-2 ; idxCluster >= 0 ; --idxCluster){
            fprintf(fxml,"\t<cluster id=\"%d\">\n", -idxCluster-1);
            fprintf(fxml,"\t\t<size>%d</size>\n", sizeOfClusters[idxCluster]);
            fprintf(fxml,"\t\t<distance>%e</distance>\n", croot[idxCluster].distance);
            fprintf(fxml,"\t\t<child id=\"%d\" direction=\"left\"/>\n", croot[idxCluster].left);
            fprintf(fxml,"\t\t<child id=\"%d\" direction=\"right\"/>\n", croot[idxCluster].right);
            fprintf(fxml,"\t</cluster>;\n");
        }

        for(int idxUnk = 0 ; idxUnk < dim ; ++idxUnk){
            fprintf(fxml,"\t<cluster id=\"%d\">\n", idxUnk);
            fprintf(fxml,"\t\t<size>%d</size>\n", 1);
            fprintf(fxml,"\t\t<indexes>%d</indexes>\n", idxUnk);
            fprintf(fxml,"\t</cluster>;\n");
        }

        fprintf(fxml,"</data>\n");
        fclose(fxml);
    }

    void saveToDot(const char inDirname[], const char inFilename[]) const{
        char buffer[1024];
        sprintf(buffer, "%s/%s", inDirname, inFilename);
        saveToDot(buffer);
    }

    void saveToDot(const char inFilename[]) const{
        std::vector<int> sizeOfClusters(dim-1, 0);
        for(int idxCluster = 0 ; idxCluster < dim-1 ; ++idxCluster){
            if(croot[idxCluster].left < 0){
                sizeOfClusters[idxCluster] += sizeOfClusters[-1-croot[idxCluster].left] ;
            }
            else{
                sizeOfClusters[idxCluster] += 1;
            }
            if(croot[idxCluster].right < 0){
                sizeOfClusters[idxCluster] += sizeOfClusters[-1-croot[idxCluster].right] ;
            }
            else{
                sizeOfClusters[idxCluster] += 1;
            }
        }
        FAssertLF(sizeOfClusters[dim-2] == dim);


        FILE* fdot = fopen(inFilename, "w");
        FAssertLF(fdot);

        fprintf(fdot,"# dot -Tsvg %s -o %s.svg\n",
                inFilename, inFilename);
        fprintf(fdot,"digraph BST {\n");

        for(int idxCluster = dim-2 ; idxCluster >= 0 ; --idxCluster){
            fprintf(fdot,"\t %d [label=\"%d Size=%d Dist=%e\"];\n",
                    (-idxCluster-1)+dim-2, idxCluster, sizeOfClusters[idxCluster], croot[idxCluster].distance);
            fprintf(fdot,"\t %d -> %d;\n", (-idxCluster-1)+dim-2, croot[idxCluster].left+dim-2);
            fprintf(fdot,"\t %d -> %d;\n", (-idxCluster-1)+dim-2, croot[idxCluster].right+dim-2);
        }

        for(int idxUnk = 0 ; idxUnk < dim ; ++idxUnk){
            fprintf(fdot,"\t%d [label=\"%d\"];\n",
                    idxUnk+dim-2, idxUnk);
        }

        fprintf(fdot,"}\n");
        fclose(fdot);
    }
};

#endif // FCCLTREECLUSTER_HPP

