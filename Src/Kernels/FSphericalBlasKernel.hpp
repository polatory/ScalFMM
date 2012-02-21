// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FSPHERICALBLASKERNEL_HPP
#define FSPHERICALBLASKERNEL_HPP

#include "FAbstractSphericalKernel.hpp"

#include "../Utils/FMemUtils.hpp"
#include "../Utils/FBlas.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* This class is a spherical harmonic kernels using blas
*/
template< class ParticleClass, class CellClass, class ContainerClass>
class FSphericalBlasKernel : public FAbstractSphericalKernel<ParticleClass,CellClass,ContainerClass> {
protected:
    typedef FAbstractSphericalKernel<ParticleClass,CellClass,ContainerClass> Parent;

    const int FF_MATRIX_ROW_DIM;     //< The blas matrix number of rows
    const int FF_MATRIX_COLUMN_DIM;  //< The blas matrix number of columns
    const int FF_MATRIX_SIZE;        //< The blas matrix size

    FComplexe* temporaryMultiSource;               //< To perform the M2L without allocating at each call
    FSmartPointer<FComplexe*> preM2LTransitions;   //< The pre-computation for the M2L based on the level and the 189 possibilities

    /** To access te precomputed M2L transfer matrixes */
    int indexM2LTransition(const int idxX,const int idxY,const int idxZ) const {
        return (((((idxX+3) * 7) + (idxY+3)) * 7 ) + (idxZ+3)) * FF_MATRIX_SIZE;
    }

    /** Alloc and init pre-vectors*/
    void allocAndInit(){
        FHarmonic blasHarmonic(Parent::devP * 2);

        // M2L transfer, there is a maximum of 3 neighbors in each direction,
        // so 6 in each dimension
        FReal treeWidthAtLevel = Parent::boxWidth * FReal( 1 << Parent::periodicLevels);
        preM2LTransitions = new FComplexe*[Parent::treeHeight + Parent::periodicLevels];
        memset(preM2LTransitions.getPtr(), 0, sizeof(FComplexe*) * (Parent::treeHeight + Parent::periodicLevels));

        for(int idxLevel = -Parent::periodicLevels ; idxLevel < Parent::treeHeight ; ++idxLevel ){
            preM2LTransitions[idxLevel + Parent::periodicLevels] = new FComplexe[(7 * 7 * 7) * FF_MATRIX_SIZE];

            for(int idxX = -3 ; idxX <= 3 ; ++idxX ){
                for(int idxY = -3 ; idxY <= 3 ; ++idxY ){
                    for(int idxZ = -3 ; idxZ <= 3 ; ++idxZ ){
                        if(FMath::Abs(idxX) > 1 || FMath::Abs(idxY) > 1 || FMath::Abs(idxZ) > 1){
                            const F3DPosition relativePos( FReal(idxX) * treeWidthAtLevel , FReal(idxY) * treeWidthAtLevel , FReal(idxZ) * treeWidthAtLevel );
                            blasHarmonic.computeOuter(FSpherical(relativePos));

                            FComplexe* FRestrict fillTransfer = &preM2LTransitions[idxLevel + Parent::periodicLevels][indexM2LTransition(idxX,idxY,idxZ)];

                            for(int M = 0 ; M <= Parent::devP ; ++M){
                                for (int m = 0 ;  m <= M ; ++m){
                                    for (int N = 0 ; N <= Parent::devP ; ++N){
                                        for (int n = 0 ; n <= 2*N ;  ++n, ++fillTransfer){
                                            const int k = N-n-m;
                                            if (k < 0){
                                                const FReal pow_of_minus_1 = FReal((k&1) ? -1 : 1);
                                                fillTransfer->setReal( pow_of_minus_1 * blasHarmonic.result()[blasHarmonic.getPreExpRedirJ(M+N)-k].getReal());
                                                fillTransfer->setImag((-pow_of_minus_1) * blasHarmonic.result()[blasHarmonic.getPreExpRedirJ(M+N)-k].getImag());
                                            }
                                            else{
                                                (*fillTransfer) = blasHarmonic.result()[blasHarmonic.getPreExpRedirJ(M+N)+k];
                                            }

                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            treeWidthAtLevel /= 2;
        }
    }


public:
    /** Constructor
      * @param inDevP the polynomial degree
      * @param inThreeHeight the height of the tree
      * @param inBoxWidth the size of the simulation box
      * @param inPeriodicLevel the number of level upper to 0 that will be requiried
      */
    FSphericalBlasKernel(const int inDevP, const int inTreeHeight, const FReal inBoxWidth, const F3DPosition& inBoxCenter, const int inPeriodicLevel = 0)
        : Parent(inDevP, inTreeHeight, inBoxWidth, inBoxCenter, inPeriodicLevel),
          FF_MATRIX_ROW_DIM(Parent::harmonic.getExpSize()), FF_MATRIX_COLUMN_DIM(Parent::harmonic.getNExpSize()),
          FF_MATRIX_SIZE(FF_MATRIX_ROW_DIM * FF_MATRIX_COLUMN_DIM),
          temporaryMultiSource(new FComplexe[FF_MATRIX_COLUMN_DIM]),
          preM2LTransitions(0){

        allocAndInit();
    }

    /** Copy constructor */
    FSphericalBlasKernel(const FSphericalBlasKernel& other)
        : Parent(other),
          FF_MATRIX_ROW_DIM(other.FF_MATRIX_ROW_DIM), FF_MATRIX_COLUMN_DIM(other.FF_MATRIX_COLUMN_DIM),
          FF_MATRIX_SIZE(other.FF_MATRIX_SIZE),
          temporaryMultiSource(new FComplexe[FF_MATRIX_COLUMN_DIM]),
          preM2LTransitions(other.preM2LTransitions) {
    }

    /** Destructor */
    ~FSphericalBlasKernel(){
        delete[] temporaryMultiSource;
        if(preM2LTransitions.isLast()){
            FMemUtils::DeleteAll(preM2LTransitions.getPtr(), Parent::treeHeight + Parent::periodicLevels);
        }
    }

    /** M2L with a cell and all the existing neighbors */
    void M2L(CellClass* const FRestrict pole, const CellClass* distantNeighbors[343],
             const int /*size*/, const int inLevel) {
        // For all neighbors compute M2L
        for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
            if( distantNeighbors[idxNeigh] ){
                const FComplexe* const transitionVector = &preM2LTransitions[inLevel + Parent::periodicLevels][idxNeigh * FF_MATRIX_SIZE];
                multipoleToLocal(pole->getLocal(), distantNeighbors[idxNeigh]->getMultipole(), transitionVector);
            }
        }
    }


    /** preExpNExp
      * @param exp an exponent vector to create an computable vector
      */
    void preExpNExp(FComplexe* const exp) const {
        for(int j = Parent::devP; j>= 0 ; --j){
            // Position in 'exp':  (j*(j+1)*0.5) + k
            // Position in 'nexp':  j*(j+1)      + k
            const int j_j1       = j*(j+1);
            const int j_j1_div_2 = int(j_j1 * 0.5);

            // Positive (or null) orders:
            for(int k = j ; k >= 0; --k){
                exp[j_j1 + k] = exp[j_j1_div_2 + k];
            }

            // Negative orders:
            FReal minus_1_pow_k = FReal( j&1 ? -1 : 1);
            for(int k = -j ; k < 0 ; ++k ){
                exp[j_j1 + k].setReal(minus_1_pow_k * exp[j_j1 + (-k)].getReal());
                exp[j_j1 + k].setImag((-minus_1_pow_k) * exp[j_j1 + (-k)].getImag());
                minus_1_pow_k = -minus_1_pow_k;
            }
        }
    }

    /** M2L
    *We compute the conversion of multipole_exp_src in *p_center_of_exp_src to
    *a local expansion in *p_center_of_exp_target, and add the result to local_exp_target.
    *
    *O_n^l (with n=0..P, l=-n..n) being the former multipole expansion terms
    *(whose center is *p_center_of_multipole_exp_src) we have for the new local
    *expansion terms (whose center is *p_center_of_local_exp_target):
    *
    *L_j^k = sum{n=0..+}
    *sum{l=-n..n}
    *O_n^l Outer_{j+n}^{-k-l}(rho, alpha, beta)
    *
    *where (rho, alpha, beta) are the spherical coordinates of the vector :
    *p_center_of_local_exp_src - *p_center_of_multipole_exp_target
    *
    *Remark: here we have always j+n >= |-k-l|
    *
    */
    void multipoleToLocal(FComplexe*const FRestrict local_exp, const FComplexe* const FRestrict multipole_exp_src,
                          const FComplexe* const FRestrict M2L_Outer_transfer){
        // Copy original vector and compute exp2nexp
        FMemUtils::copyall<FComplexe>(temporaryMultiSource, multipole_exp_src, CellClass::GetPoleSize());
        // Get a computable vector
        preExpNExp(temporaryMultiSource);

        //CblasTrans
        FBlas::c_gemtva(
                    static_cast<unsigned int>(FF_MATRIX_COLUMN_DIM),
                    static_cast<unsigned int>(FF_MATRIX_ROW_DIM),
                    FReal(1.0),
                    (FReal*)M2L_Outer_transfer,
                    (FReal*)temporaryMultiSource,
                    (FReal*)local_exp);
        // TODO delete
        static int count = 0;
        if( FMath::IsNan(local_exp[CellClass::GetLocalSize()-1].getReal()) && count != -1){
            std::cout << "count is " << count << std::endl;
            count = -1;
        }
        if(count != -1) count++;
    }
};

#endif // FSPHERICALBLASKERNEL_HPP
