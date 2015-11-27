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

