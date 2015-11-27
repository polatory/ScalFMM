#ifndef FMATGRID_HPP
#define FMATGRID_HPP

// @SCALFMM_PRIVATE

#include "Utils/FGlobal.hpp"
#include "Utils/FAssert.hpp"

#include "FDenseBlockWrapper.hpp"

template <class FReal>
class FMatGrid {
protected:
    int matDim;
    FReal* values;

    FMatGrid(const FMatGrid&) = delete;
    FMatGrid& operator=(const FMatGrid&) = delete;
public:
    using BlockDescriptor = FDenseBlockWrapper<FReal>;

    FMatGrid(const int inDim, const FReal* inValues)
        : matDim(inDim), values(nullptr){
        values = new FReal[matDim*matDim];
        for(int idxVal = 0 ; idxVal < matDim*matDim ; ++idxVal){
            values[idxVal] = inValues[idxVal];
        }
    }

    ~FMatGrid(){
        delete[] values;
    }

    FDenseBlockWrapper<FReal> getBlock(const int rowIdx, const int colIdx, const int nbRows, const int nbCols) const {
        // static_assert(std::is_move_constructible<BlockClass>::value, "The block class must be movable");
        // static_assert(std::is_move_assignable<BlockClass>::value, "The block class must be movable");
        FAssertLF(0 < rowIdx);
        FAssertLF(0 < colIdx);
        FAssertLF(rowIdx + nbRows <= matDim);
        FAssertLF(colIdx + nbRows <= matDim);
        return FDenseBlockWrapper<FReal>(&values[colIdx*matDim+rowIdx], nbRows, nbCols, matDim);
    }
};

#endif // FMATGRID_HPP

