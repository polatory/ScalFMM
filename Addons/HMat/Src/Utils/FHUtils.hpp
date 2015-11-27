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
#ifndef FHUTILS_HPP
#define FHUTILS_HPP

// @SCALFMM_PRIVATE

#include "Utils/FGlobal.hpp"

#include <cstring>

template <class Type>
void FSetToZeros(Type array[], const int length){
    memset(array, 0, length*sizeof(Type));
}


struct FBlockDescriptor {
    int row, col, nbRows, nbCols, level;
};

#endif // FHUTILS_HPP
