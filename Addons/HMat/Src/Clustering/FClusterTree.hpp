#ifndef FCLUSTERTREE_HPP
#define FCLUSTERTREE_HPP

// @SCALFMM_PRIVATE
#include "./Utils/FGlobal.hpp"
#include "./Containers/FBoolArray.hpp"
#include "./Utils/FAssert.hpp"

#include <queue>
#include <vector>

template <class FReal>
class FClusterTree {
public:
    struct Node {
        int id;
        int left;
        int right;
        int parent;
        int size;
        FReal score;
    };

    struct Leaf {
        int id;
        int parent;
        int size;
        int offset;
    };

protected:

    FClusterTree(const FClusterTree&) = delete;
    FClusterTree& operator=(const FClusterTree&) = delete;

    int dim;
    int nbClusters;
    const Node* clusters;
    int nbLeaves;
    const Leaf* leaves;
    const int* permsOrigToNew;
    const int* permsNewToOrig;

    void release(){
        dim = 0;
        nbClusters = 0;
        nbLeaves   = 0;
        delete[] clusters;
        delete[] leaves;
        delete[] permsOrigToNew;
        delete[] permsNewToOrig;
    }

public:
    FClusterTree()
        : dim(0), nbClusters(0), clusters(nullptr), nbLeaves(0), leaves(nullptr),
           permsOrigToNew(nullptr), permsNewToOrig(nullptr) {

    }

    ~FClusterTree(){
        release();
    }

    bool isInit() const {
        return clusters != nullptr && leaves != nullptr;
    }

    void setData(const int inDim,
                 const int inNbClusters, const Node* inClusters,
                 const int inNbLeaves, const Leaf* inLeaves,
                 const int* inPermsOrigToNew, const int* inPermsNewToOrig){
        release();
        dim = inDim;
        nbClusters = inNbClusters;
        clusters = inClusters;
        nbLeaves = inNbLeaves;
        leaves = inLeaves;
        permsOrigToNew = inPermsOrigToNew;
        permsNewToOrig = inPermsNewToOrig;
    }

    void checkData() const{
        { // test permutations
            for(int idx = 0 ; idx < dim ; ++idx){
                FAssertLF(permsNewToOrig[permsOrigToNew[idx]] == idx);
            }
        }
        FAssertLF(clusters[nbClusters-1].size == dim);
        {
            FBoolArray unkTests(dim);
            FAssertLF(nbClusters < dim);
            for(int idx = 0 ; idx < nbClusters ; ++idx){
                const Node& nd = clusters[idx];
                FAssertLF(nd.id == (-idx)-1);
                FAssertLF(0 <= nd.size);
                FAssertLF(-nbClusters-1 < nd.left && nd.left < nbLeaves);
                FAssertLF(-nbClusters-1 < nd.right && nd.right < nbLeaves);

                int size = 0;
                if(0 <= nd.left){
                    const Leaf& lf = leaves[nd.left];
                    size += lf.size;
                    FAssertLF(nd.id == lf.parent);
                    FAssertLF(lf.size + lf.offset <= dim);
                    for(int idxUnk = 0 ; idxUnk < lf.size ; ++idxUnk){
                        FAssertLF(unkTests.get(idxUnk + lf.offset) == false);
                        unkTests.set(idxUnk + lf.offset, true);
                    }
                }
                else{
                    FAssertLF(nd.left < idx);
                    size += clusters[(-nd.left)-1].size;
                    FAssertLF(nd.id == clusters[(-nd.left)-1].parent);
                }
                if(0 <= nd.right){
                    const Leaf& lf = leaves[nd.right];
                    size += lf.size;
                    FAssertLF(nd.id == lf.parent);
                    FAssertLF(lf.size + lf.offset <= dim);
                    for(int idxUnk = 0 ; idxUnk < lf.size ; ++idxUnk){
                        FAssertLF(unkTests.get(idxUnk + lf.offset) == false);
                        unkTests.set(idxUnk + lf.offset, true);
                    }
                }
                else{
                    FAssertLF(nd.right < idx);
                    size += clusters[(-nd.right)-1].size;
                    FAssertLF(nd.id == clusters[(-nd.right)-1].parent,
                            nd.id ,"==", clusters[(-nd.right)-1].parent);
                }
                FAssertLF(size == nd.size);
            }

            for(int idxUnk = 0 ; idxUnk < dim ; ++idxUnk){
                FAssertLF(unkTests.get(idxUnk) == true);
            }
        }
        { // test permutations
            int offset = 0;
            for(int idx = 0 ; idx < nbLeaves ; ++idx){
                FAssertLF(leaves[idx].id == idx);
                FAssertLF((clusters[-1-leaves[idx].parent].left == idx)
                            ^ (clusters[-1-leaves[idx].parent].right == idx));
                FAssertLF(leaves[idx].offset == offset);
                offset += leaves[idx].size;
            }
            FAssertLF(dim == offset);
        }
    }

    void fillPermutations(int* inPermuts, int* invPermuts = nullptr) const {
        memcpy(inPermuts, permsOrigToNew, sizeof(int)*dim);
        if(invPermuts){
            memcpy(invPermuts, permsNewToOrig, sizeof(int)*dim);
        }
    }

    int getPartitions(const int inHeight, const int inNbPartitions, int inNbIdxInPartitions[]) const{
        // Here we select the inNbPartitions last partitions
        // But we have to take them in the right order (left to right)
        // To ensure coherency with the permutations
        // if inNbPartitions, we should have inNbIdxInPartitions filled with 1
        FAssertLF(FMath::pow2(inHeight-1)  == inNbPartitions);

        const int noNodeFlag = std::numeric_limits<int>::max();
        std::queue<int> breadthFirst;
        breadthFirst.push(clusters[dim-2].left);
        breadthFirst.push(clusters[dim-2].right);

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
                    breadthFirstLower.push(clusters[(-current)-1].left);
                    breadthFirstLower.push(clusters[(-current)-1].right);
                }
            }

            breadthFirst = std::move(breadthFirstLower);
        }
        FAssertLF(int(breadthFirst.size()) == inNbPartitions);

        int counterPartition = 0;
        int totalGiven = 0;
        int emptyPartitions = 0;
        while(breadthFirst.size()){
            const int current = breadthFirst.front();
            breadthFirst.pop();

            if(current == noNodeFlag){
                inNbIdxInPartitions[counterPartition] = 0;
                emptyPartitions += 1;
            }
            else if(0 <= current){
                inNbIdxInPartitions[counterPartition] = leaves[current].size;
            }
            else{
                inNbIdxInPartitions[counterPartition] = clusters[(-current)-1].size;
            }
            totalGiven +=inNbIdxInPartitions[counterPartition];
            counterPartition += 1;
        }
        FAssertLF(totalGiven == dim);

        return emptyPartitions;
    }


    void saveToXml(const char inDirname[], const char inFilename[]) const{
        char buffer[1024];
        sprintf(buffer, "%s/%s", inDirname, inFilename);
        saveToXml(buffer);
    }

    void saveToXml(const char inFilename[]) const{
        FILE* fxml = fopen(inFilename, "w");
        FAssertLF(fxml);

        fprintf(fxml,"<?xml version=\"1.0\"?>\n");
        fprintf(fxml,"<data>\n");

        for(int idxCluster = 0 ; idxCluster < nbClusters ; ++idxCluster){
            fprintf(fxml,"\t<cluster id=\"%d\">\n", clusters[idxCluster].id);
            fprintf(fxml,"\t\t<size>%d</size>\n", clusters[idxCluster].size);
            fprintf(fxml,"\t\t<score>%e</score>\n", clusters[idxCluster].score);
            fprintf(fxml,"\t\t<child id=\"%d\" direction=\"left\"/>\n", clusters[idxCluster].left);
            fprintf(fxml,"\t\t<child id=\"%d\" direction=\"right\"/>\n", clusters[idxCluster].right);
            fprintf(fxml,"\t</cluster>;\n");
        }

        for(int idxUnk = 0 ; idxUnk < nbLeaves ; ++idxUnk){
            fprintf(fxml,"\t<cluster id=\"%d\">\n", idxUnk);
            fprintf(fxml,"\t\t<size>%d</size>\n", leaves[idxUnk].size);
            fprintf(fxml,"\t\t<offset>%d</offset>\n", leaves[idxUnk].offset);
            fprintf(fxml,"\t\t<indexes>\n");
            for(int idxIndex = 0 ; idxIndex < leaves[idxUnk].size ; ++idxIndex){
                fprintf(fxml,"<unknown>%d</unknown>\n", permsNewToOrig[leaves[idxUnk].offset + idxIndex]);
            }
            fprintf(fxml,"\t\t</indexes>\n");
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
        FILE* fdot = fopen(inFilename, "w");
        FAssertLF(fdot);

        fprintf(fdot,"# dot -Tsvg %s -o %s.svg\n",
                inFilename, inFilename);
        fprintf(fdot,"digraph BST {\n");

        for(int idxCluster = 0 ; idxCluster < nbClusters ; ++idxCluster){
            FAssertLF(idxCluster == -(clusters[idxCluster].id+1));
            fprintf(fdot,"\t %d [label=\"%d Size=%d\", tooltip=\"score %e\"];\n",
                    idxCluster, idxCluster, clusters[idxCluster].size, clusters[idxCluster].score);
            fprintf(fdot,"\t %d -> %d;\n", idxCluster, (-clusters[idxCluster].left)-2+dim);
            fprintf(fdot,"\t %d -> %d;\n", idxCluster, (-clusters[idxCluster].right)-2+dim);
        }

        for(int idxUnk = 0 ; idxUnk < nbLeaves ; ++idxUnk){
            fprintf(fdot,"\t%d [label=\"%d size = %d offset = %d\"];\n",
                    idxUnk+dim-2, idxUnk, leaves[idxUnk].size, leaves[idxUnk].offset);
        }

        fprintf(fdot,"}\n");
        fclose(fdot);
    }
};


#endif // FCLUSTERTREE_HPP

