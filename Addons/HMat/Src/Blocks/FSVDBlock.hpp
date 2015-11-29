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
#ifndef FSVDBLOCK_HPP
#define FSVDBLOCK_HPP

#include "Utils/FBlas.hpp"

/*
 * Compute SVD $A=USV'$ and return $V'$, $S$ and $U$.
 * \param A contains input M x N matrix to be decomposed
 * \param S contains singular values $S$
 * \param U contains $U$
 * \param VT contains $V'$
 */
template<class FReal>
static void computeSVD(const FSize nbRows, const FSize nbCols, const FReal* A, FReal* S, FReal* U, FReal* VT){
    // verbose
    const bool verbose = false;
    // copy A
    //is_int(size*size);
    FBlas::copy(int(nbRows*nbCols),A,U);
    // init SVD
    const FSize minMN = std::min(nbRows,nbCols);
    const FSize maxMN = std::max(nbRows,nbCols);
    //const FSize LWORK = 2*4*minMN; // for square matrices
    const FSize LWORK = 2*std::max(3*minMN+maxMN, 5*minMN);
    FReal *const WORK = new FReal [LWORK];
    // singular value decomposition
    if(verbose) std::cout << "\nPerform SVD...";
    // SO means that first min(m,n) lines of U overwritten on VT and V' on U (A=VTSU)
    // AA means all lines
    // nothing means OS (the opposite of SO, A=USVT)
    //is_int(size); is_int(LWORK);
    const unsigned int INFOSVD
    = FBlas::gesvd(int(nbRows), int(nbCols), U, S, VT, int(minMN),
                   int(LWORK), WORK);
    if(verbose) {
        if(INFOSVD!=0) {std::cout << " failed!" << std::endl;}
        else {std::cout << " succeed!" << std::endl;}
    }
    // free memory
    delete[] WORK;
}

template <class FReal>
class FSVDBlock{
protected:
    // members
    FReal* block;
    FReal* U;
    FReal* S;
    FReal* VT;
    int nbRows;
    int nbCols;
    int level;
    int rank;
    FSVDBlock(const FSVDBlock&) = delete;
    FSVDBlock& operator=(const FSVDBlock&) = delete;

public:
    FSVDBlock()
        : block(nullptr), U(nullptr), S(nullptr), VT(nullptr), nbRows(0), nbCols(0),  level(0), rank(0) {
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

        // SVD specific (col major)
        rank = std::min(nbRows,nbCols);
        U  = new FReal[nbRows*nbCols]; // Call to computeSVD() copies block into U
        S  = new FReal[rank];
        VT = new FReal[rank*nbCols];
        FBlas::setzero(int(rank), S);        
        FBlas::setzero(int(rank*nbCols),VT);
        // Perform decomposition of rectangular block (jobu=O, jobvt=S => only first min(M,N) cols/rows of U/VT are stored)
        computeSVD(nbRows, nbCols, block, S, U ,VT);

        // TODO Determine numerical rank using prescribed accuracy
        // ...

    };

    // dtor
    ~FSVDBlock(){
        // Free memory
        clear();
    };

    void clear(){
        nbRows = 0;
        nbCols = 0;
        level = 0;
        rank = 0;
        delete[] block;
        block = 0;
        delete[] U;
        U = 0;
        delete[] S;
        S = 0;
        delete[] VT;
        VT = 0;
    }

    void gemv(FReal res[], const FReal vec[], const FReal scale = FReal(1.)) const {
        
        //// Apply (dense) block
        //FBlas::gemva(nbRows, nbCols, scale, const_cast<FReal*>(block), const_cast<FReal*>(vec), res);
        
        // Apply low-rank block
        FReal* VTvec = new FReal[nbCols];
        // Apply VT
        FBlas::gemv(rank, nbCols, scale, const_cast<FReal*>(VT), const_cast<FReal*>(vec), VTvec);
        // Apply S
        for(int idxS = 0 ; idxS < rank ; ++idxS)
            VTvec[idxS]*=S[idxS];
        // Apply U
        FBlas::gemva(nbRows, rank, scale, const_cast<FReal*>(U), const_cast<FReal*>(VTvec), res);

    }

    void gemm(FReal res[], const FReal mat[], const int nbRhs, const FReal scale = FReal(1.)) const {
        
        //// Apply (dense) block
        //FBlas::gemma(nbRows, nbCols, nbRhs, scale, const_cast<FReal*>(block), nbRows, const_cast<FReal*>(mat), nbCols, res, nbRows);

        // Apply low-rank block
        FReal* VTmat = new FReal[nbCols*nbRhs];
        // Apply VT
        FBlas::gemm(rank, nbCols, nbRhs, scale, const_cast<FReal*>(VT), rank, const_cast<FReal*>(mat), nbCols, VTmat, rank);
        // Apply S
        for(int idxRow = 0 ; idxRow < rank ; ++idxRow)
            for(int idxRhs = 0 ; idxRhs < nbRhs ; ++idxRhs)
                VTmat[idxRhs*rank+idxRow]*=S[idxRow];
        // Apply U
        FBlas::gemma(nbRows, rank, nbRhs, scale, const_cast<FReal*>(U), nbRows, const_cast<FReal*>(VTmat), rank, res, nbRows);

    }


};

#endif // FSVDBLOCK_HPP

// [--END--]
