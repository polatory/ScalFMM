#ifndef FDIV2BISSECTION_HPP
#define FDIV2BISSECTION_HPP

// @SCALFMM_PRIVATE

#include "Utils/FGlobal.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FAssert.hpp"

#include "../Utils/FHUtils.hpp"

#include <functional>



template <class FReal, class LeafClass, class CellClass >
class FDiv2Bissection {
protected:
    struct LeafNode {
        FBlockDescriptor infos;
        LeafClass leaf;
    };

    struct CellNode {
        FBlockDescriptor infos;
        CellClass cell;
    };

    const int dim;
    const int height;

    int* nbCells;
    CellNode** cells;

    int nbLeaves;
    LeafNode* leaves;

    int totalNbBlocks;

    FDiv2Bissection(const FDiv2Bissection&) = delete;
    FDiv2Bissection& operator=(const FDiv2Bissection&) = delete;

public:
    explicit FDiv2Bissection(const int inDim, const int inHeight)
        : dim(inDim), height(inHeight),
          nbCells(0), cells(nullptr),
          nbLeaves(0), leaves(nullptr),
          totalNbBlocks(0){
        FAssertLF(FMath::pow2(height) <= inDim);

        cells   = new CellNode*[height-1];
        FSetToZeros(cells, height-1);
        nbCells = new int[height-1];
        FSetToZeros(nbCells, height-1);

        for(int idxLevel = 1 ; idxLevel < height-1 ; ++idxLevel){
            const int nbCellsAtLevel = FMath::pow2(idxLevel);
            cells[idxLevel]   = new CellNode[nbCellsAtLevel];
            nbCells[idxLevel] = nbCellsAtLevel;
            totalNbBlocks += nbCellsAtLevel;

            const int nbCellsInDirection = nbCells[idxLevel];
            for(int idxCell = 0 ; idxCell < nbCellsAtLevel ; ++idxCell){
                const int rowCellNumber = (idxCell&1? idxCell-1 : idxCell+1);
                const int colCellNumner = idxCell;
                cells[idxLevel][idxCell].infos.row = ((rowCellNumber*dim)/nbCellsInDirection);
                cells[idxLevel][idxCell].infos.col = ((colCellNumner*dim)/nbCellsInDirection);
                cells[idxLevel][idxCell].infos.nbRows = (((rowCellNumber+1)*dim)/nbCellsInDirection)
                                                        - cells[idxLevel][idxCell].infos.row;
                cells[idxLevel][idxCell].infos.nbCols =(((colCellNumner+1)*dim)/nbCellsInDirection)
                                                        - cells[idxLevel][idxCell].infos.col;
                cells[idxLevel][idxCell].infos.level = idxLevel;
            }
        }

        nbLeaves = FMath::pow2(height);
        leaves   = new LeafNode[nbLeaves];
        totalNbBlocks += nbLeaves;
        {
            const int nbLeavesInDirection = (nbLeaves/2);
            for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
                const int rowLeafNumber = ((idxLeaf/4)*2) + (idxLeaf&1?1:0);
                const int colLeafNumber = ((idxLeaf/4)*2) + (idxLeaf&2?1:0);
                leaves[idxLeaf].infos.row = ((rowLeafNumber*dim)/nbLeavesInDirection);
                leaves[idxLeaf].infos.col = ((colLeafNumber*dim)/nbLeavesInDirection);
                leaves[idxLeaf].infos.nbRows = (((rowLeafNumber+1)*dim)/nbLeavesInDirection)
                                                - leaves[idxLeaf].infos.row;
                leaves[idxLeaf].infos.nbCols = (((colLeafNumber+1)*dim)/nbLeavesInDirection)
                                                - leaves[idxLeaf].infos.col;
                leaves[idxLeaf].infos.level = height-1;
            }
        }
    }

    ~FDiv2Bissection(){
        for(int idxLevel = 0 ; idxLevel < height-1 ; ++idxLevel){
            delete[] cells[idxLevel];
        }
        delete[] cells;
        delete[] nbCells;
        delete[] leaves;
    }

    int getNbBlocks() const {
        return totalNbBlocks;
    }

    void forAllBlocksDescriptor(std::function<void(const FBlockDescriptor&)> callback){
        for(int idxLevel = 1 ; idxLevel < height-1 ; ++idxLevel){
            for(int idxCell = 0 ; idxCell < nbCells[idxLevel] ; ++idxCell){
                callback(cells[idxLevel][idxCell].infos);
            }
        }
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            callback(leaves[idxLeaf].infos);
        }
    }

    void forAllCellBlocks(std::function<void(const FBlockDescriptor&, CellClass&)> callback){
        for(int idxLevel = 1 ; idxLevel < height-1 ; ++idxLevel){
            for(int idxCell = 0 ; idxCell < nbCells[idxLevel] ; ++idxCell){
                callback(cells[idxLevel][idxCell].infos,
                         cells[idxLevel][idxCell].cell);
            }
        }
    }

    void forAllLeafBlocks(std::function<void(const FBlockDescriptor&, LeafClass&)> callback){
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            callback(leaves[idxLeaf].infos,
                     leaves[idxLeaf].leaf);
        }
    }

    template <class MatrixClass>
    void fillBlocks(MatrixClass& matrix){
        for(int idxLevel = 1 ; idxLevel < height-1 ; ++idxLevel){
            for(int idxCell = 0 ; idxCell < nbCells[idxLevel] ; ++idxCell){
                cells[idxLevel][idxCell].fill(matrix.getBlock(
                                                  cells[idxLevel][idxCell].infos.row,
                                                  cells[idxLevel][idxCell].infos.col,
                                                  cells[idxLevel][idxCell].infos.nbRows,
                                                  cells[idxLevel][idxCell].infos.nbCols
                                                  ));
            }
        }
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            leaves[idxLeaf].fill(matrix.getBlock(
                                      leaves[idxLeaf].infos.row,
                                      leaves[idxLeaf].infos.col,
                                      leaves[idxLeaf].infos.nbRows,
                                      leaves[idxLeaf].infos.nbCols
                                      ));
        }
    }

    template <class MatrixClass>
    void fillBlocks(const MatrixClass& matrix){
        for(int idxLevel = 1 ; idxLevel < height-1 ; ++idxLevel){
            for(int idxCell = 0 ; idxCell < nbCells[idxLevel] ; ++idxCell){
                cells[idxLevel][idxCell].fill(matrix.getBlock(
                                                  cells[idxLevel][idxCell].infos.row,
                                                  cells[idxLevel][idxCell].infos.col,
                                                  cells[idxLevel][idxCell].infos.nbRows,
                                                  cells[idxLevel][idxCell].infos.nbCols
                                                  ));
            }
        }
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            leaves[idxLeaf].fill(matrix.getBlock(
                                      leaves[idxLeaf].infos.row,
                                      leaves[idxLeaf].infos.col,
                                      leaves[idxLeaf].infos.nbRows,
                                      leaves[idxLeaf].infos.nbCols
                                      ));
        }
    }

    void gemv(FReal res[], const FReal vec[]) const {
        for(int idxLevel = 1 ; idxLevel < height-1 ; ++idxLevel){
            for(int idxCell = 0 ; idxCell < nbCells[idxLevel] ; ++idxCell){
                cells[idxLevel][idxCell].gemv(&res[cells[idxLevel][idxCell].infos.col],
                                              &vec[cells[idxLevel][idxCell].infos.row]);
            }
        }
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            leaves[idxLeaf].gemv(&res[leaves[idxLeaf].infos.col],
                                          &vec[leaves[idxLeaf].infos.row]);
        }
    }

    void gemm(FReal res[], const FReal mat[], const int nbRhs) const {
        for(int idxLevel = 1 ; idxLevel < height-1 ; ++idxLevel){
            for(int idxCell = 0 ; idxCell < nbCells[idxLevel] ; ++idxCell){
                cells[idxLevel][idxCell].gemm(&res[cells[idxLevel][idxCell].infos.col],
                                              &mat[cells[idxLevel][idxCell].infos.row],
                                              nbRhs, dim);
            }
        }
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            leaves[idxLeaf].gemm(&res[leaves[idxLeaf].infos.col],
                                 &mat[leaves[idxLeaf].infos.row],
                                 nbRhs, dim);
        }
    }
};


#endif // FDIV2BISSECTION_HPP

