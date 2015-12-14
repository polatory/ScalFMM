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
#ifndef FACABLOCK_HPP
#define FACABLOCK_HPP

#include "Utils/FBlas.hpp"
#include "Utils/FAca.hpp"


/*!  Choose either \a FULLY_PIVOTED_ACA or \a PARTIALLY_PIVOTED_ACA */
#define FULLY_PIVOTED_ACA
//#define PARTIALLY_PIVOTED_ACA

/*!  Choose \a RECOMPRESSED_ACA */
//#define RECOMPRESSED_ACA


template <class FReal, int ORDER = 14>
class FACABlock{
protected:
    // members
    FReal* block;
    FReal* U;
    FReal* VT;
    int nbRows;
    int nbCols;
    int level;
    unsigned int rank;
    FReal accuracy;
    FACABlock(const FACABlock&) = delete;
    FACABlock& operator=(const FACABlock&) = delete;

public:
    FACABlock()
        : block(nullptr), U(nullptr), VT(nullptr), nbRows(0), nbCols(0),  level(0), rank(0), accuracy(FMath::pow(10.0,static_cast<FReal>(-ORDER))) {
    }

    // ctor
    template <class ViewerClass>
    void fill(const ViewerClass& viewer, const int inLevel){
        clear();
        // Allocate memory
        level  = inLevel;
        nbRows = viewer.getNbRows();
        nbCols = viewer.getNbCols();


#if (defined FULLY_PIVOTED_ACA)
        block  = new FReal[nbRows*nbCols];

        for(int idxRow = 0 ; idxRow < nbRows ; ++idxRow){
            for(int idxCol = 0 ; idxCol < nbCols ; ++idxCol){
                block[idxCol*nbRows+idxRow] = viewer.getValue(idxRow,idxCol);
            }
        }
#endif


        // SVD specific (col major)
        rank = std::min(nbRows,nbCols);
        //S  = new FReal[rank];
        FReal* _U ;// = new FReal[nbRows*nbCols]; // Call to computeSVD() copies block into U
        FReal* _V;// = new FReal[rank*nbCols];

#if (defined FULLY_PIVOTED_ACA)

        // Perform fully pivoted ACA
        FAca::fACA<FReal>(block, nbRows, nbCols, accuracy, _U, _V, rank);

#elif (defined PARTIAL_PIVOTED_ACA)

        // TODO
        // Perform partially pivoted ACA
        //FAca::fACA<FReal>(viewer, nbRows, nbCols, accuracy, _U, _V, rank);

#endif
        // display rank
        std::cout << "rank after ACA=" << rank << " (" << nbRows << "," << nbCols << ")" << std::endl;

#if (defined RECOMPRESSED_ACA)

        // Recompression by QR+SVD
        FAca::recompress(_U,_V,nbRows,nbCols,accuracy,U,VT,rank);
        // display rank after recompression
        std::cout << "  rank after QR+SVD=" << rank << std::endl;

#else

        // Resize U and VT
        U  = new FReal[nbRows*rank]; 
        for(int idxRow = 0 ; idxRow < nbRows ; ++idxRow)
            for(int idxCol = 0 ; idxCol < rank ; ++idxCol)
                U[idxCol*nbRows+idxRow] = _U[idxCol*nbRows+idxRow];
  
        VT = new FReal[rank*nbCols];
        for(int idxRow = 0 ; idxRow < rank ; ++idxRow)
            for(int idxCol = 0 ; idxCol < nbCols ; ++idxCol)
                VT[idxCol*rank+idxRow] = _V[idxRow*nbCols+idxCol]; // transposition

#endif

        // Free memory
        delete [] _U;
        delete [] _V;

    };

    // dtor
    ~FACABlock(){
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
        // Apply U
        FBlas::gemma(nbRows, rank, nbRhs, scale, const_cast<FReal*>(U), nbRows, const_cast<FReal*>(VTmat), rank, res, nbRows);

    }

    int getRank() const{
        return rank;
    }

};

#endif // FACABLOCK_HPP

// [--END--]
