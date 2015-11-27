// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FPARTITIONSMAPPING_HPP
#define FPARTITIONSMAPPING_HPP


// @SCALFMM_PRIVATE

#include "Utils/FGlobal.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FAssert.hpp"

#include "../Utils/FHUtils.hpp"

#include <functional>
#include <memory>


template <class FReal, class LeafClass, class CellClass >
class FPartitionsMapping {
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

    int nbCells;
    CellNode* cells;

    int nbLeaves;
    LeafNode* leaves;

    int totalNbBlocks;

    FPartitionsMapping(const FPartitionsMapping&) = delete;
    FPartitionsMapping& operator=(const FPartitionsMapping&) = delete;

public:
    explicit FPartitionsMapping(const int inDim, const int partitions[], const int nbPartitions,
                                const int ratioForLeaf = 0)
        : dim(inDim),
          nbCells(0), cells(nullptr),
          nbLeaves(0), leaves(nullptr),
          totalNbBlocks(0){
        FAssertLF(nbPartitions <= inDim);
        FAssertLF(1 <= nbPartitions);

        for(int idxPartRow = 0 ; idxPartRow < nbPartitions ; ++idxPartRow){
            for(int idxPartCol = 0 ; idxPartCol < nbPartitions ; ++idxPartCol){
                if(idxPartRow == idxPartCol
                    || partitions[idxPartRow]*partitions[idxPartCol] < ratioForLeaf){
                    nbLeaves += 1;
                }
                else{
                    nbCells += 1;
                }
                totalNbBlocks += 1;
            }
        }

        leaves   = new LeafNode[nbLeaves];
        cells    = new CellNode[nbCells];

        int idxLeaf = 0;
        int idxCell = 0;
        int offsetRows = 0;
        for(int idxPartRow = 0 ; idxPartRow < nbPartitions ; ++idxPartRow){
            int offsetCols = 0;
            for(int idxPartCol = 0 ; idxPartCol < nbPartitions ; ++idxPartCol){
                if(idxPartRow == idxPartCol
                    || partitions[idxPartRow]*partitions[idxPartCol] < ratioForLeaf){
                    leaves[idxLeaf].infos.row = offsetRows;
                    leaves[idxLeaf].infos.col = offsetCols;
                    leaves[idxLeaf].infos.nbRows = partitions[idxPartRow];
                    leaves[idxLeaf].infos.nbCols = partitions[idxPartCol];
                    leaves[idxLeaf].infos.level = 0;
                    idxLeaf += 1;
                }
                else{
                    cells[idxCell].infos.row = offsetRows;
                    cells[idxCell].infos.col = offsetCols;
                    cells[idxCell].infos.nbRows = partitions[idxPartRow];
                    cells[idxCell].infos.nbCols = partitions[idxPartCol];
                    cells[idxCell].infos.level = 1;
                    idxCell += 1;
                }
                offsetCols += partitions[idxPartCol];
            }
            FAssertLF(offsetCols == dim);
            offsetRows += partitions[idxPartRow];
        }
        FAssertLF(offsetRows == dim);
    }

    ~FPartitionsMapping(){
        delete[] cells;
        delete[] leaves;
    }

    int getNbBlocks() const {
        return totalNbBlocks;
    }

    void forAllBlocksDescriptor(std::function<void(const FBlockDescriptor&)> callback){
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            callback(cells[idxCell].infos);
        }
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            callback(leaves[idxLeaf].infos);
        }
    }

    void forAllCellBlocks(std::function<void(const FBlockDescriptor&, CellClass&)> callback){
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            callback(cells[idxCell].infos,
                     cells[idxCell].cell);
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
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            cells[idxCell].cell.fill(matrix.getBlock(
                                              cells[idxCell].infos.row,
                                              cells[idxCell].infos.col,
                                              cells[idxCell].infos.nbRows,
                                              cells[idxCell].infos.nbCols
                                              ),
                                               cells[idxCell].infos.level);
        }
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            leaves[idxLeaf].leaf.fill(matrix.getBlock(
                                      leaves[idxLeaf].infos.row,
                                      leaves[idxLeaf].infos.col,
                                      leaves[idxLeaf].infos.nbRows,
                                      leaves[idxLeaf].infos.nbCols
                                      ),
                                      leaves[idxLeaf].infos.level);
        }
    }

    template <class MatrixClass>
    void fillBlocks(const MatrixClass& matrix){
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            cells[idxCell].cell.fill(matrix.getBlock(
                                              cells[idxCell].infos.row,
                                              cells[idxCell].infos.col,
                                              cells[idxCell].infos.nbRows,
                                              cells[idxCell].infos.nbCols
                                              ),
                                          cells[idxCell].infos.level);
        }
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            leaves[idxLeaf].leaf.fill(matrix.getBlock(
                                      leaves[idxLeaf].infos.row,
                                      leaves[idxLeaf].infos.col,
                                      leaves[idxLeaf].infos.nbRows,
                                      leaves[idxLeaf].infos.nbCols
                                      ),
                                 leaves[idxLeaf].infos.level);
        }
    }

    void gemv(FReal res[], const FReal vec[]) const {
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            cells[idxCell].cell.gemv(&res[cells[idxCell].infos.row],
                                          &vec[cells[idxCell].infos.col]);
        }
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            leaves[idxLeaf].leaf.gemv(&res[leaves[idxLeaf].infos.row],
                                          &vec[leaves[idxLeaf].infos.col]);
        }
    }

    void gemm(FReal res[], const FReal mat[], const int nbRhs) const {
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            cells[idxCell].cell.gemm(&res[cells[idxCell].infos.col],
                                          &mat[cells[idxCell].infos.row],
                                          nbRhs, dim);
        }
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            leaves[idxLeaf].leaf.gemm(&res[leaves[idxLeaf].infos.col],
                                 &mat[leaves[idxLeaf].infos.row],
                                 nbRhs, dim);
        }
    }
};

#endif // FPARTITIONSMAPPING_HPP

