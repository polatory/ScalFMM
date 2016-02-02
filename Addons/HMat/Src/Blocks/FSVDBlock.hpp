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
    = FBlas::gesvd(int(nbRows), int(nbCols), U, S, VT, int(minMN)/*ldVT*/,
                   int(LWORK), WORK);
    if(verbose) {
        if(INFOSVD!=0) {std::cout << " failed!" << std::endl;}
        else {std::cout << " succeed!" << std::endl;}
    }
    // free memory
    delete[] WORK;
}


/*
 * Compute SVD $A=USV'$ and return $V'$, $S$ and $U$.
 * \param 
 */
template<class FReal>
static void computeNumericalRank(int &rank, const FReal* S, const FReal epsilon){
    // verbose
    const bool verbose = false;

    // init
    const FSize maxRank = rank;
    FReal sumSigma2 = FReal(0.0);
    for(int idxRow = 0 ; idxRow < rank ; ++idxRow)
        sumSigma2+=S[idxRow]*S[idxRow];
    FReal SqrtSumSigma2 = std::sqrt(sumSigma2);

    // set rank to 1
    rank = 1;
    // increase
    FReal sumSigma2r = S[0]*S[0];
    while(std::sqrt(sumSigma2r)<(FReal(1.)-epsilon)*SqrtSumSigma2 && rank<maxRank){
        sumSigma2r+=S[rank]*S[rank];
        rank++;
    }
    //std::cout << "std::sqrt(sumSigma2r)=" << std::sqrt(sumSigma2r) << std::endl;
    //std::cout << "std::sqrt(sumSigma2)=" << std::sqrt(sumSigma2) << std::endl;
    //std::cout << "R/S=" << (std::sqrt(sumSigma2)-std::sqrt(sumSigma2r))/std::sqrt(sumSigma2) << std::endl;

}




template <class FReal, int ORDER = 14>
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
    FReal accuracy;
    FSVDBlock(const FSVDBlock&) = delete;
    FSVDBlock& operator=(const FSVDBlock&) = delete;

public:
    FSVDBlock()
        : block(nullptr), U(nullptr), S(nullptr), VT(nullptr), nbRows(0), nbCols(0),  level(0), rank(0), accuracy(FMath::pow(FReal(10.0),static_cast<FReal>(-ORDER))) {
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
        S  = new FReal[rank];
        FReal* _U  = new FReal[nbRows*nbCols]; // Call to computeSVD() copies block MxN into _U and stores first min(M,N) cols of U into _U
        FReal* _VT = new FReal[rank*nbCols];
        FBlas::setzero(int(rank), S);        
        FBlas::setzero(int(rank*nbCols),_VT);
        // Perform decomposition of rectangular block (jobu=O, jobvt=S => only first min(M,N) cols/rows of U/VT are stored)
        computeSVD(nbRows, nbCols, block, S, _U ,_VT);

        // Determine numerical rank using prescribed accuracy
        computeNumericalRank(rank, S, accuracy);

        //// Display singular values
        //std::cout << "S = [";
        //for(int idxRow = 0 ; idxRow < rank ; ++idxRow)
        //    std::cout << S[idxRow] << " " ;
        //std::cout << "]" << std::endl;

        //// display rank
        std::cout << "rank SVD =" << rank << " (" << nbRows << "," << nbCols << ")" << std::endl;

        // Resize U and VT
        U  = new FReal[nbRows*rank]; // Call to computeSVD() copies block into U
        VT = new FReal[rank*nbCols];

        for(int idxRow = 0 ; idxRow < nbRows ; ++idxRow)
            for(int idxCol = 0 ; idxCol < rank ; ++idxCol)
                U[idxCol*nbRows+idxRow] = _U[idxCol*nbRows+idxRow];
        for(int idxRow = 0 ; idxRow < rank ; ++idxRow)
            for(int idxCol = 0 ; idxCol < nbCols ; ++idxCol)
                VT[idxCol*rank+idxRow] = _VT[idxCol*std::min(nbRows,nbCols)+idxRow]; 
        delete [] _U;           
        delete [] _VT;

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
        //FReal* res_dense = new FReal[nbRows];
        //FBlas::copy(nbRows,res,res_dense);
        //FBlas::gemva(nbRows, nbCols, scale, const_cast<FReal*>(block), const_cast<FReal*>(vec), res_dense);
        
        // Apply low-rank block
        FReal* VTvec = new FReal[rank];
        FBlas::setzero(rank,VTvec);

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

    int getRank() const{
        return rank;
    }
};

#endif // FSVDBLOCK_HPP

// [--END--]
