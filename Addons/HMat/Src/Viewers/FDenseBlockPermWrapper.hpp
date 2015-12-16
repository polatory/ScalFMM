#ifndef FDENSEBLOCKPERMWRAPPER_HPP
#define FDENSEBLOCKPERMWRAPPER_HPP

// @SCALFMM_PRIVATE

template <class FReal, class SrcMatrixClass >
class FDenseBlockPermWrapper{
protected:
    const SrcMatrixClass& matrix;
    const int row;
    const int col;
    const int nbRows;
    const int nbCols;

public:
    FDenseBlockPermWrapper(const SrcMatrixClass& inMatrix, const int inRow, const int inCol,
                           const int inNbRows, const int inNbCols)
        : matrix(inMatrix), row(inRow), col(inCol),
          nbRows(inNbRows), nbCols(inNbCols) {
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
        return matrix.getVal(row+rowIdx, col+colIdx);
    }

    constexpr bool existsForReal() const{
        return true;
    }
};

#endif // FDENSEBLOCKPERMWRAPPER_HPP

