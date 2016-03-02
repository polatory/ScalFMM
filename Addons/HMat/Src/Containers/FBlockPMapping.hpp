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


template <class FReal, class CellClass >
class FBlockPMapping {
protected:
    struct CellCNode {
        FBlockDescriptor infos;
        CellClass cell;
    };

    struct RowUNode {
        FBlockDescriptor infos;
        CellClass cell;
    };

    struct ColVNode {
        FBlockDescriptor infos;
        CellClass cell;
    };

    const int dim;
    const int nbPartitions;
    const int nbCells;

    CellCNode* cBlocks;
    RowUNode* uRowBlocks;
    ColVNode* vColBlocks;

    FBlockPMapping(const FBlockPMapping&) = delete;
    FBlockPMapping& operator=(const FBlockPMapping&) = delete;

public:
    explicit FBlockPMapping(const int inDim, const int partitions[], const int inNbPartitions)
        : dim(inDim),
          nbPartitions(inNbPartitions),
          nbCells(inNbPartitions*inNbPartitions),
          cBlocks(nullptr){
        FAssertLF(nbPartitions <= inDim);
        FAssertLF(1 <= nbPartitions);

        std::unique_ptr<int[]> partitionsOffset(new int[nbPartitions]);
        partitionsOffset[0] = 0;
        for(int idxPart = 1 ; idxPart < nbPartitions ; ++idxPart){
            partitionsOffset[idxPart] = partitionsOffset[idxPart-1] + partitions[idxPart-1];
        }

        cBlocks    = new CellCNode[nbCells];

        for(int idxPartCol = 0 ; idxPartCol < nbPartitions ; ++idxPartCol){
            for(int idxPartRow = 0 ; idxPartRow < nbPartitions ; ++idxPartRow){
                cBlocks[idxPartCol*nbPartitions + idxPartRow].infos.row = partitionsOffset[idxPartRow];
                cBlocks[idxPartCol*nbPartitions + idxPartRow].infos.col = partitionsOffset[idxPartCol];
                cBlocks[idxPartCol*nbPartitions + idxPartRow].infos.nbRows = partitions[idxPartRow];
                cBlocks[idxPartCol*nbPartitions + idxPartRow].infos.nbCols = partitions[idxPartCol];
                cBlocks[idxPartCol*nbPartitions + idxPartRow].infos.level = 0;
            }
        }

        uRowBlocks = new RowUNode[nbPartitions];
        for(int idxPartRow = 0 ; idxPartRow < nbPartitions ; ++idxPartRow){
            uRowBlocks[idxPartRow].infos.row = partitionsOffset[idxPartRow];
            uRowBlocks[idxPartRow].infos.col = 0;
            uRowBlocks[idxPartRow].infos.nbRows = partitions[idxPartRow];
            uRowBlocks[idxPartRow].infos.nbCols = dim;
            uRowBlocks[idxPartRow].infos.level = 0;
        }

        vColBlocks = new ColVNode[nbPartitions];
        for(int idxPartCol = 0 ; idxPartCol < nbPartitions ; ++idxPartCol){
            vColBlocks[idxPartCol].infos.row = 0;
            vColBlocks[idxPartCol].infos.col = partitionsOffset[idxPartCol];
            vColBlocks[idxPartCol].infos.nbRows = dim;
            vColBlocks[idxPartCol].infos.nbCols = partitions[idxPartCol];
            vColBlocks[idxPartCol].infos.level = 0;
        }
    }

    ~FBlockPMapping(){
        delete[] cBlocks;
        delete[] uRowBlocks;
        delete[] vColBlocks;
    }

    int getNbBlocks() const {
        return nbCells;
    }

    // Iterate blocks

    CellClass& getCBlock(const int idxRowPart, const int idxColPart){
        return cBlocks[idxColPart*nbPartitions + idxRowPart].cell;
    }

    const CellClass& getCBlock(const int idxRowPart, const int idxColPart) const {
        return cBlocks[idxColPart*nbPartitions + idxRowPart].cell;
    }

    const FBlockDescriptor& getCBlockInfo(const int idxRowPart, const int idxColPart) const {
        return cBlocks[idxColPart*nbPartitions + idxRowPart].infos;
    }

    void forAllCBlocksDescriptor(std::function<void(const FBlockDescriptor&)> callback){
        for(int idxCell = 0 ; idxCell < nbCells ; ++idxCell){
            callback(cBlocks[idxCell].infos);
        }
    }

    void forAllBlocks(std::function<void(const FBlockDescriptor&,
                          CellClass&, CellClass&, CellClass&)> callback){
        for(int idxPartCol = 0 ; idxPartCol < nbPartitions ; ++idxPartCol){
            for(int idxPartRow = 0 ; idxPartRow < nbPartitions ; ++idxPartRow){
                callback(cBlocks[idxPartCol*nbPartitions + idxPartRow].infos,
                     cBlocks[idxPartCol*nbPartitions + idxPartRow].cell,
                     uRowBlocks[idxPartRow].cell,
                     vColBlocks[idxPartCol].cell);
            }
        }
    }

    // Iterate row blocks

    CellClass& getUBlock(const int idxRowPart){
        return uRowBlocks[idxRowPart].cell;
    }

    const CellClass& getUBlock(const int idxRowPart) const {
        return uRowBlocks[idxRowPart].cell;
    }

    const FBlockDescriptor& getUBlockInfo(const int idxRowPart) const {
        return uRowBlocks[idxRowPart].infos;
    }


    // Iterate col blocks

    CellClass& getVBlock(const int idxColPart){
        return vColBlocks[idxColPart].cell;
    }

    const CellClass& getVBlock(const int idxColPart) const {
        return vColBlocks[idxColPart].cell;
    }

    const FBlockDescriptor& getVBlockInfo(const int idxColPart) const {
        return vColBlocks[idxColPart].infos;
    }

    // Operations
    void gemv(FReal res[], const FReal vec[]) const {
        for(int idxPartCol = 0 ; idxPartCol < nbPartitions ; ++idxPartCol){
            for(int idxPartRow = 0 ; idxPartRow < nbPartitions ; ++idxPartRow){
//                &res[cBlocks[idxPartCol*nbPartitions + idxPartRow].infos.row],
//                &vec[cBlocks[idxPartCol*nbPartitions + idxPartRow].infos.col])
//                cBlocks[idxPartCol*nbPartitions + idxPartRow].cell,
//                uRowBlocks[idxPartRow].cell,
//                vColBlocks[idxPartCol].cell;
            }
        }
    }

    void gemm(FReal res[], const FReal mat[], const int nbRhs) const {
        for(int idxPartCol = 0 ; idxPartCol < nbPartitions ; ++idxPartCol){
            for(int idxPartRow = 0 ; idxPartRow < nbPartitions ; ++idxPartRow){
//                &res[cBlocks[idxPartCol*nbPartitions + idxPartRow].infos.row],
//                &vec[cBlocks[idxPartCol*nbPartitions + idxPartRow].infos.col])
//                cBlocks[idxPartCol*nbPartitions + idxPartRow].cell,
//                uRowBlocks[idxPartRow].cell,
//                vColBlocks[idxPartCol].cell;
//                nbRhs, dim
            }
        }
    }
};

#endif // FBLOCKPMAPPING_HPP


