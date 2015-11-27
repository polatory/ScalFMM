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
class FDenseMatrix{
protected:
    // members
    FReal* block;
    int nbRows;
    int nbCols;
    int level;

    FDenseMatrix(const FDenseMatrix&) = delete;
    FDenseMatrix& operator=(const FDenseMatrix&) = delete;

public:
    FDenseMatrix()
        : block(nullptr), nbRows(0), nbCols(0),  level(0) {
    }

    // ctor
    template <class ViewerClass>
    void fill(const ViewerClass& viewer, const int inLevel){
        clear();
        // Allocate memory
        level  = inLevel;
        nbRows = viewer.getNbRows();
        nbCols = viewer.getNbCols();
        block  = new FReal[nbRows*nbCols];

        for(int idxRow = 0 ; idxRow < nbRows ; ++idxRow){
            for(int idxCol = 0 ; idxCol < nbCols ; ++idxCol){
                block[idxCol*nbRows+idxRow] = viewer.getValue(idxRow,idxCol);
            }
        }
    };

    // dtor
    ~FDenseMatrix(){
        // Free memory
        clear();
    };

    void clear(){
        nbRows = 0;
        nbCols = 0;
        nbCols = 0;
        delete[] block;
        block = 0;
    }

    void gemv(FReal res[], const FReal vec[], const FReal scale = FReal(1.)) const {
        FBlas::gemva(nbRows, nbCols, scale, const_cast<FReal*>(block), const_cast<FReal*>(vec), res);
    }

    void gemm(FReal res[], const FReal vec[], const int nbrhs, const FReal scale = FReal(1.)) const {
        FBlas::gemva(nbRows, nbCols, scale, const_cast<FReal*>(block), const_cast<FReal*>(vec), res);
    }
};

#endif // FDENSEMATRIX_HPP

// [--END--]
