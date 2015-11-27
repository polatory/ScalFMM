#ifndef FDENSEBLOCKWRAPPER_HPP
#define FDENSEBLOCKWRAPPER_HPP

// @SCALFMM_PRIVATE

template <class FReal>
class FDenseBlockWrapper{
protected:
    const FReal* values;
    const int nbRows;
    const int nbCols;
    const int leadingDim;

public:
    FDenseBlockWrapper(const FReal* inValues, const int inNbRows, const int inNbCols, const int inLeading)
        : values(inValues), nbRows(inNbRows), nbCols(inNbCols), leadingDim(inLeading){
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
        return values[colIdx*leadingDim + rowIdx];
    }

    constexpr bool existsForReal() const{
        return true;
    }
};

#endif // FDENSEBLOCKWRAPPER_HPP

