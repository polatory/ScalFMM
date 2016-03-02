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
#ifndef FMATDENSE_HPP
#define FMATDENSE_HPP

// @SCALFMM_PRIVATE

#include "Utils/FGlobal.hpp"
#include "Utils/FAssert.hpp"

#include "FDenseBlockWrapper.hpp"
#include "../Utils/FHUtils.hpp"
#include "../Utils/FMatrixIO.hpp"

template <class FReal>
class FMatDense {
protected:
    int matDim;
    FReal* values;

    FMatDense(const FMatDense&) = delete;
    FMatDense& operator=(const FMatDense&) = delete;
public:
    using BlockDescriptor = FDenseBlockWrapper<FReal>;

    explicit FMatDense(const int inDim, const FReal* inValues = nullptr)
        : matDim(inDim), values(nullptr){
        values = new FReal[matDim*matDim];
        if(inValues){
            for(int idxVal = 0 ; idxVal < matDim*matDim ; ++idxVal){
                values[idxVal] = inValues[idxVal];
            }
        }
        else{
            FSetToZeros(values, matDim*matDim);
        }
    }

    explicit FMatDense(const char inFilename[])
        : matDim(0), values(nullptr){
        int readNbRows = 0;
        int readNbCols = 0;
        FAssertLF(FMatrixIO::read(inFilename, &values, &readNbRows, &readNbCols));
        FAssertLF(readNbRows == readNbCols);
        matDim = readNbRows;
    }

    ~FMatDense(){
        delete[] values;
    }

    int getDim() const{
        return matDim;
    }

    const FReal& getVal(const int idxRow , const int idxCol) const{
        return values[idxCol*matDim+idxRow];
    }

    FReal& getVal(const int idxRow , const int idxCol) {
        return values[idxCol*matDim+idxRow];
    }

    void setVal(const int idxRow , const int idxCol, const FReal& val) {
        values[idxCol*matDim+idxRow] = val;
    }

    FDenseBlockWrapper<FReal> getBlock(const int rowIdx, const int colIdx, const int nbRows, const int nbCols) const {
        // static_assert(std::is_move_constructible<BlockClass>::value, "The block class must be movable");
        // static_assert(std::is_move_assignable<BlockClass>::value, "The block class must be movable");
        FAssertLF(0 < nbRows);
        FAssertLF(0 < nbCols);
        FAssertLF(rowIdx + nbRows <= matDim);
        FAssertLF(colIdx + nbRows <= matDim);
        return FDenseBlockWrapper<FReal>(&values[colIdx*matDim+rowIdx], nbRows, nbCols, matDim);
    }

    FDenseBlockWrapper<FReal> getBlock(const FBlockDescriptor& info) const {
        return getBlock(info.row, info.col, info.nbRows, info.nbCols);
    }

};

#endif // FMATDENSE_HPP

