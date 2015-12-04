#ifndef FDENSEBLOCKPERMWRAPPER_HPP
#define FDENSEBLOCKPERMWRAPPER_HPP

// @SCALFMM_PRIVATE

template <class FReal>
class FDenseBlockPermWrapper{
protected:
    const FReal* values;
    const int nbRows;
    const int nbCols;
    const int leadingDim;
    const int* permuteVector;

public:
    FDenseBlockPermWrapper(const FReal* inValues, const int inNbRows, const int inNbCols,
                           const int inLeading, const int* inPermuteVector)
        : values(inValues), nbRows(inNbRows), nbCols(inNbCols), leadingDim(inLeading), permuteVector(inPermuteVector){
    }

    int getNbRows() const {
        return nbRows;
    }

    int getNbCols() const {
        return nbCols;
    }

    FReal getValue(const int rowIdx, const int colIdx) const{
        FAssertLF(rowIdx < nbRows);
        FAssertLF(colIdx < nbCols);
        return values[permuteVector[colIdx]*leadingDim + permuteVector[rowIdx]];
    }

    constexpr bool existsForReal() const{
        return true;
    }
};

#endif // FDENSEBLOCKPERMWRAPPER_HPP

