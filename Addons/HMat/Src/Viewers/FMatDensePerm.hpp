#ifndef FMATDENSEPERM_HPP
#define FMATDENSEPERM_HPP


// @SCALFMM_PRIVATE

#include "Utils/FGlobal.hpp"
#include "Utils/FAssert.hpp"

#include "FDenseBlockPermWrapper.hpp"
#include "../Utils/FHUtils.hpp"
#include "../Utils/FMatrixIO.hpp"

template <class FReal>
class FMatDense {
protected:
    int matDim;
    FReal* values;
    int* permOrigToNew;

    FMatDense(const FMatDense&) = delete;
    FMatDense& operator=(const FMatDense&) = delete;
public:
    using BlockDescriptor = FDenseBlockPermWrapper<FReal>;

    explicit FMatDense(const int inDim, const FReal* inValues = nullptr, const int* inPermOrigToNew = nullptr)
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

        permOrigToNew = new int[matDim];
        if(inPermOrigToNew){
            memcpy(permOrigToNew, inPermOrigToNew, sizeof(int)*matDim);
        }
        else {
            for(int idxVal = 0 ; idxVal < matDim ; ++idxVal){
                permOrigToNew[idxVal] = idxVal;
            }
        }
    }

    explicit FMatDense(const char inFilename[], const int* inPermOrigToNew = nullptr)
        : matDim(0), values(nullptr){
        int readNbRows = 0;
        int readNbCols = 0;
        FAssertLF(FMatrixIO::read(inFilename, &values, &readNbRows, &readNbCols));
        FAssertLF(readNbRows == readNbCols);
        matDim = readNbRows;

        permOrigToNew = new int[matDim];
        if(inPermOrigToNew){
            memcpy(permOrigToNew, inPermOrigToNew, sizeof(int)*matDim);
        }
        else {
            for(int idxVal = 0 ; idxVal < matDim ; ++idxVal){
                permOrigToNew[idxVal] = idxVal;
            }
        }
    }

    ~FMatDense(){
        delete[] values;
        delete[] permOrigToNew;
    }

    void setPermutOrigToNew(const int* inPermOrigToNew){
        memcpy(permOrigToNew, inPermOrigToNew, sizeof(int)*matDim);
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

    FDenseBlockPermWrapper<FReal> getBlock(const int rowIdx, const int colIdx, const int nbRows, const int nbCols) const {
        // static_assert(std::is_move_constructible<BlockClass>::value, "The block class must be movable");
        // static_assert(std::is_move_assignable<BlockClass>::value, "The block class must be movable");
        FAssertLF(0 < nbRows);
        FAssertLF(0 < nbCols);
        FAssertLF(rowIdx + nbRows <= matDim);
        FAssertLF(colIdx + nbRows <= matDim);
        return FDenseBlockPermWrapper<FReal>(&values[colIdx*matDim+rowIdx], nbRows, nbCols, matDim, permOrigToNew);
    }
};

#endif // FMATDENSEPERM_HPP

