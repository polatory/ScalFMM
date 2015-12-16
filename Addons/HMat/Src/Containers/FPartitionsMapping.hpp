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


template <class FReal, class CellClass >
class FPartitionsMapping {
protected:
    struct CellNode {
        FBlockDescriptor infos;
        CellClass cell;
    };

    const int dim;
    const int nbPartitions;
    const int nbCells;

    CellNode* cells;

    FPartitionsMapping(const FPartitionsMapping&) = delete;
    FPartitionsMapping& operator=(const FPartitionsMapping&) = delete;

public:
    explicit FPartitionsMapping(const int inDim, const int partitions[], const int inNbPartitions)
        : dim(inDim),
          nbPartitions(inNbPartitions),
          nbCells(inNbPartitions*inNbPartitions),
          cells(nullptr){
        FAssertLF(nbPartitions <= inDim);
        FAssertLF(1 <= nbPartitions);

        std::unique_ptr<int[]> partitionsOffset(new int[nbPartitions]);
        partitionsOffset[0] = 0;
        for(int idxPart = 1 ; idxPart < nbPartitions ; ++idxPart){
            partitionsOffset[idxPart] = partitionsOffset[idxPart-1] + partitions[idxPart-1];
        }

        cells    = new CellNode[nbCells];

        for(int idxPartCol = 0 ; idxPartCol < nbPartitions ; ++idxPartCol){
            for(int idxPartRow = 0 ; idxPartRow < nbPartitions ; ++idxPartRow){
                cells[idxPartCol*nbPartitions + idxPartRow].infos.row = partitionsOffset[idxPartRow];
                cells[idxPartCol*nbPartitions + idxPartRow].infos.col = partitionsOffset[idxPartCol];
                cells[idxPartCol*nbPartitions + idxPartRow].infos.nbRows = partitions[idxPartRow];
                cells[idxPartCol*nbPartitions + idxPartRow].infos.nbCols = partitions[idxPartCol];
                cells[idxPartCol*nbPartitions + idxPartRow].infos.level = 0;
            }
        }
    }

    ~FPartitionsMapping(){
        delete[] cells;
    }

    int getNbBlocks() const {
        return nbCells;
    }

    CellClass& getCell(const int idxRowPart, const int idxColPart){
        return cells[idxColPart*nbPartitions + idxRowPart].cell;
    }

    const CellClass& getCell(const int idxRowPart, const int idxColPart) const {
        return cells[idxColPart*nbPartitions + idxRowPart].cell;
    }

    const FBlockDescriptor& getCellInfo(const int idxRowPart, const int idxColPart) const {
        return cells[idxColPart*nbPartitions + idxRowPart].infos;
    }

    void forAllBlocksDescriptor(std::function<void(const FBlockDescriptor&)> callback){
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            callback(cells[idxCell].infos);
        }
    }

    void forAllCellBlocks(std::function<void(const FBlockDescriptor&, CellClass&)> callback){
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            callback(cells[idxCell].infos,
                     cells[idxCell].cell);
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
    }

    void gemv(FReal res[], const FReal vec[]) const {
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            cells[idxCell].cell.gemv(&res[cells[idxCell].infos.row],
                                          &vec[cells[idxCell].infos.col]);
        }
    }

    void gemm(FReal res[], const FReal mat[], const int nbRhs) const {
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            cells[idxCell].cell.gemm(&res[cells[idxCell].infos.col],
                                          &mat[cells[idxCell].infos.row],
                                          nbRhs, dim);
        }
    }
};

#endif // FPARTITIONSMAPPING_HPP

