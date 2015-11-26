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
// 
// @SCALFMM_PRIVATE
// 
#ifndef FDENSEMATRIX_HPP
#define FDENSEMATRIX_HPP

#include "Utils/FBlas.hpp"

template <class FReal>
class FDenseMatrix
{

protected:
    // members
    FReal* block;
    const int height;
    const int width;
    const int leading_dimension;
    const int level;

public:

    // ctor
    explicit FDenseMatrix(const FReal* in_block, const int in_height, const int in_width, const int in_leading_dimension /*h<N*/, const int in_level)
    : height(in_height), width(in_width), leading_dimension(in_leading_dimension), level(in_level)
    {

        // Allocate memory
        block = new FReal[height*width];
        // Copy
        FBlas::copy(height*width,in_block,block);

    };
    // dtor
    ~FDenseMatrix(){

        // Free memory
        delete[] block;

    };

    void MV(FReal res[], const FReal vec[], const FReal scale = FReal(1.)) {

        FBlas::gemv(height, width, scale, block, vec, res);

    }

    
};

#endif // FDENSEMATRIX_HPP

// [--END--]
