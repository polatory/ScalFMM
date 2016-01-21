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
#ifndef FBLOCKPMAPPING_HPP
#define FBLOCKPMAPPING_HPP


// @SCALFMM_PRIVATE

#include "Utils/FGlobal.hpp"
#include "Utils/FMath.hpp"
#include "Utils/FAssert.hpp"

#include "../Utils/FHUtils.hpp"

#include <functional>
#include <memory>


template <class FReal, class RowBlockClass, class ColBlockClass, class CoreCellClass >
class FBlockPMapping {
protected:
    struct CellNode {
        FBlockDescriptor infos;
        CoreCellClass cell;
    };

    struct RowNode {
        FBlockDescriptor infos;
        RowBlockClass cell;
    };

    struct ColNode {
        FBlockDescriptor infos;
        ColBlockClass cell;
    };

    const int dim;
    const int nbPartitions;
    const int nbCells;

    CellNode* cells;
    RowNode* rowBlocks;
    ColNode* colBlocks;

    FBlockPMapping(const FBlockPMapping&) = delete;
    FBlockPMapping& operator=(const FBlockPMapping&) = delete;

public:
    explicit FBlockPMapping(const int inDim, const int partitions[], const int inNbPartitions)
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

        rowBlocks = new RowNode[nbPartitions];
        for(int idxPartRow = 0 ; idxPartRow < nbPartitions ; ++idxPartRow){
            rowBlocks[idxPartRow].infos.row = partitionsOffset[idxPartRow];
            rowBlocks[idxPartRow].infos.col = 0;
            rowBlocks[idxPartRow].infos.nbRows = partitions[idxPartRow];
            rowBlocks[idxPartRow].infos.nbCols = dim;
            rowBlocks[idxPartRow].infos.level = 0;
        }

        colBlocks = new ColNode[nbPartitions];
        for(int idxPartCol = 0 ; idxPartCol < nbPartitions ; ++idxPartCol){
            colBlocks[idxPartRow].infos.row = 0;
            colBlocks[idxPartRow].infos.col = partitionsOffset[idxPartCol];
            colBlocks[idxPartRow].infos.nbRows = dim;
            colBlocks[idxPartRow].infos.nbCols = partitions[idxPartRow];
            colBlocks[idxPartRow].infos.level = 0;
        }
    }

    ~FBlockPMapping(){
        delete[] cells;
        delete[] rowBlocks;
        delete[] colBlocks;
    }

    int getNbBlocks() const {
        return nbCells;
    }

    // Iterate blocks

    CoreCellClass& getCell(const int idxRowPart, const int idxColPart){
        return cells[idxColPart*nbPartitions + idxRowPart].cell;
    }

    const CoreCellClass& getCell(const int idxRowPart, const int idxColPart) const {
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

    void forAllCellBlocks(std::function<void(const FBlockDescriptor&,
                          RowBlockClass&, CoreCellClass&, ColBlockClass&)> callback){
        for(int idxPartCol = 0 ; idxPartCol < nbPartitions ; ++idxPartCol){
            for(int idxPartRow = 0 ; idxPartRow < nbPartitions ; ++idxPartRow){
                callback(cells[idxPartCol*nbPartitions + idxPartRow].infos,
                     cells[idxPartCol*nbPartitions + idxPartRow].cell,
                     rowBlocks[idxPartRow],
                     colBlocks[idxPartCol]);
            }
        }
    }

    // Iterate row blocks

    RowBlockClass& getRowCell(const int idxRowPart){
        return rowBlocks[idxRowPart].cell;
    }

    const RowBlockClass& getRowCell(const int idxRowPart) const {
        return rowBlocks[idxRowPart].cell;
    }

    const FBlockDescriptor& getRowCellInfo(const int idxRowPart) const {
        return rowBlocks[idxRowPart].infos;
    }


    // Iterate col blocks

    ColBlockClass& getColCell(const int idxColPart){
        return colBlocks[idxColPart].cell;
    }

    const ColBlockClass& getColCell(const int idxColPart) const {
        return colBlocks[idxColPart].cell;
    }

    const ColBlockClass& getColCellInfo(const int idxColPart) const {
        return colBlocks[idxColPart].infos;
    }

    // Operations
};

#endif // FBLOCKPMAPPING_HPP


