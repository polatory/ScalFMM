#ifndef FMATDENSE_HPP
#define FMATDENSE_HPP

// @SCALFMM_PRIVATE

#include "Utils/FGlobal.hpp"
#include "Utils/FAssert.hpp"

#include "FDenseBlockWrapper.hpp"
#include "../Utils/FHUtils.hpp"

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

    ~FMatDense(){
        delete[] values;
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
};

#endif // FMATDENSE_HPP

