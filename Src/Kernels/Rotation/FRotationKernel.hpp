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
#ifndef FROTATIONKERNEL_HPP
#define FROTATIONKERNEL_HPP

#include "../../Components/FAbstractKernels.hpp"
#include "../../Utils/FSmartPointer.hpp"
#include "../../Utils/FComplexe.hpp"
#include "../../Utils/FMemUtils.hpp"
#include "../../Utils/FSpherical.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FRotationKernel
* @brief
*
* This kernels is a complete rotation based kernel with spherical
* harmonic.
*
* Here is the optimizated kernel, please refer to FRotationOriginalKernel
* to see the non optimized easy to understand kernel.
*/
template< class ParticleClass, class CellClass, class ContainerClass, int P>
class FRotationKernel : public FAbstractKernels<ParticleClass,CellClass,ContainerClass> {    
    //< Size of the data array computed using a suite relation
    static const int SizeArray = ((P+2)*(P+1))/2;
    static const int P2 = P*2;

    ///////////////////////////////////////////////////////
    // Object attributes
    ///////////////////////////////////////////////////////

    const FReal boxWidth;               //< the box width at leaf level
    const int   treeHeight;             //< The height of the tree
    const FReal widthAtLeafLevel;       //< width of box at leaf level
    const FReal widthAtLeafLevelDiv2;   //< width of box at leaf leve div 2
    const FPoint boxCorner;             //< position of the box corner


    FSmartPointer<FSpherical> childrenPosition;     //< The transfers between level
    FSmartPointer<FSpherical> interactionsPosition; //< The transfers at a same level

    FReal factorials[P2+1];

    FSmartPointer<FReal[P+1]>      M2MTranslationCoef;
    FSmartPointer<FReal[343][P+1]> M2LTranslationCoef;
    FSmartPointer<FReal[P+1]>      L2LTranslationCoef;

    FComplexe rotationExpMinusImPhi[8][SizeArray];
    FComplexe rotationExpImPhi[8][SizeArray];

    FSmartPointer<FComplexe[343][SizeArray]> M2LAroundZCoef;

    ///////////////////////////////////////////////////////
    // d_lmk
    ///////////////////////////////////////////////////////

    /** Return position in the array of the l/m couple */
    int atLm(const int l, const int m){
        // summation series over l + m => (l*(l+1))/2 + m
        return ((l*(l+1))>>1) + m;
    }

    /** Return the sign of a value */
    template < class T >
    T Sgn(const T a){
        if(a < 0) return T(-1);
        else if(a > 0) return T(1);
        return T(0);
    }

    /** Compute the factorial from 0 to P*2 */
    void precomputeFactorials(){
        factorials[0] = 1;
        FReal fidx = 1;
        for(int idx = 1 ; idx <= P2 ; ++idx, ++fidx){
            factorials[idx] = fidx * factorials[idx-1];
        }
    }

    void precomputeTranslationCoef(){
        M2MTranslationCoef = new FReal[treeHeight-1][P+1];
        L2LTranslationCoef = new FReal[treeHeight-1][P+1];

        FReal widthAtLevel = boxWidth/4;
        for( int idxLevel = 0 ; idxLevel < treeHeight - 1 ; ++idxLevel){
            FReal b = FMath::Sqrt(widthAtLevel*widthAtLevel*3);
            FReal bPowIdx = 1.0;
            FReal minus_1_pow_idx = 1.0;
            for(int idx = 0 ; idx <= P ; ++idx){
                // coef m2m = (-b)^j/j!
                M2MTranslationCoef[idxLevel][idx] = minus_1_pow_idx * bPowIdx / factorials[idx];
                // coef l2l = b^j/j!
                L2LTranslationCoef[idxLevel][idx] = bPowIdx / factorials[idx];
                // increase
                bPowIdx *= b;
                minus_1_pow_idx = -minus_1_pow_idx;
            }
            widthAtLevel /= 2;
        }

        M2LTranslationCoef = new FReal[treeHeight][343][P+1];
        for( int idxLevel = 0 ; idxLevel < treeHeight ; ++idxLevel){
            for(int idxInteraction = 0 ; idxInteraction < 343 ; ++idxInteraction){
                int x = (idxInteraction/7*7)-3;
                int y = ((idxInteraction - (x+3)*7*7)/7) - 3;
                int z = idxInteraction - ((x+3)*7 + (y+3))*7;//idxInteraction%7;
                if( x < -1 || 1 < x || y < -1 || 1 < y || z < -1 || 1 < 1 ){
                    const FReal b = getSphericalInteraction(idxLevel, idxInteraction).getR();
                    FReal bPowIdx1 = b;
                    for(int idx = 0 ; idx <= P ; ++idx){
                        // factorials[j+l] / FMath::pow(b,j+l+1)
                        M2LTranslationCoef[idxLevel][idxInteraction][idx] = factorials[idx] / bPowIdx1;
                        bPowIdx1 *= b;
                    }
                }
            }
        }
    }

    ///////////////////////////////////////////////////////
    // Position precomputation
    ///////////////////////////////////////////////////////

    /** This function precompute the position between cells */
    void preComputePosition(){
        // Compute the parent/children relation for M2M L2L
        childrenPosition = new FSpherical[8 * (treeHeight-1)];
        {
            FReal subBoxWidth = widthAtLeafLevelDiv2;
            // For each
            for(int idxLevel = treeHeight-2 ; idxLevel > 0 ; --idxLevel){
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild ){
                    // coord from child to parent
                    const FReal x = FReal((idxChild&4)? -subBoxWidth : subBoxWidth);
                    const FReal y = FReal((idxChild&2)? -subBoxWidth : subBoxWidth);
                    const FReal z = FReal((idxChild&1)? -subBoxWidth : subBoxWidth);
                    const FPoint relativePosition( x , y , z );

                    childrenPosition[(idxLevel-1)*8 + idxChild] = FSpherical(relativePosition);
                }
                subBoxWidth *= FReal(2.0);
            }
        }
        // Compute the interaction relations for M2L
        interactionsPosition = new FSpherical[343 * (treeHeight-1)];
        {
            FReal boxWidthAtLevel = widthAtLeafLevel;
            for(int idxLevel = treeHeight-1 ; idxLevel > 0 ; --idxLevel){
                for(int idxX = -3 ; idxX <= 3 ; ++idxX ){
                    for(int idxY = -3 ; idxY <= 3 ; ++idxY ){
                        for(int idxZ = -3 ; idxZ <= 3 ; ++idxZ ){
                            if( idxX != 0 || idxY != 0 || idxZ != 0 ){
                                const FPoint relativePosition( -FReal(idxX)*boxWidthAtLevel,
                                                               -FReal(idxY)*boxWidthAtLevel,
                                                               -FReal(idxZ)*boxWidthAtLevel);
                                const int position = ((( (idxX+3) * 7) + (idxY+3))) * 7 + idxZ + 3;
                                interactionsPosition[(idxLevel-1)*343 + position] = FSpherical(relativePosition);
                            }
                        }
                    }
                }
                boxWidthAtLevel *= FReal(2.0);
            }
        }
    }

    /** This function rotate a multipole vector by an angle azimuth phi
      * The formula used is present in several paper, but we refer to
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1
      * \f[
      * O_{l,m}( \alpha, \beta + \phi ) = e^{-i \phi m} O_{l,m}( \alpha, \beta )
      * \f]
      * The computation is simply a multiplication per a complex number \f$ e^{-i \phi m} \f$
      * Phi should be in [0,2pi]
      */
    /** This function rotate a local vector by an angle azimuth phi
      * The formula used is present in several paper, but we refer to
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1
      * \f[
      * M_{l,m}( \alpha, \beta + \phi ) = e^{i \phi m} M_{l,m}( \alpha, \beta )
      * \f]
      * The computation is simply a multiplication per a complex number \f$ e^{i \phi m} \f$
      * Phi should be in [0,2pi]
      */

    void precomputeRotationVectors(){
        const int index_P0 = atLm(P,0);
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            const FReal x = FReal((idxChild&4)? -boxWidth : boxWidth);
            const FReal y = FReal((idxChild&2)? -boxWidth : boxWidth);
            const FReal z = FReal((idxChild&1)? -boxWidth : boxWidth);
            const FPoint relativePosition( x , y , z );

            FSpherical sph(relativePosition);

            // compute the last part with l == P
            {
                int index_lm = index_P0;
                for(int m = 0 ; m <= P ; ++m, ++index_lm ){
                    const FReal mphi = (sph.getAzimuth() + FMath::FPiDiv2) * FReal(m);
                    // O_{l,m}( \alpha, \beta + \phi ) = e^{-i \phi m} O_{l,m}( \alpha, \beta )
                    rotationExpMinusImPhi[idxChild][index_lm].setRealImag(FMath::Cos(-mphi), FMath::Sin(-mphi));
                    // M_{l,m}( \alpha, \beta + \phi ) = e^{i \phi m} M_{l,m}( \alpha, \beta )
                    rotationExpImPhi[idxChild][index_lm].setRealImag(FMath::Cos(mphi), FMath::Sin(mphi));
                }
            }
            // Then copy
            {
                int index_lm = 0;
                for(int l = 0 ; l < P ; ++l){
                    FMemUtils::copyall(rotationExpMinusImPhi[idxChild] + index_lm,
                                       rotationExpMinusImPhi[idxChild] + index_P0,
                                       l + 1);
                    FMemUtils::copyall(rotationExpImPhi[idxChild] + index_lm,
                                       rotationExpImPhi[idxChild] + index_P0,
                                       l + 1);
                    index_lm += l + 1;
                }
            }
        }
    }

    /** To get the right spherical object from level and child position */
    FSpherical getSphericalChild(const int idxLevel, const int position) const {
        return childrenPosition[(idxLevel-1)*8 + position];
    }

    /** To get the right spherical object from level and interaction position */
    FSpherical getSphericalInteraction(const int idxLevel, const int position) const {
        return interactionsPosition[(idxLevel-1)*343 + position];
    }

    /** Return the position of a leaf from its tree coordinate */
    FPoint getLeafCenter(const FTreeCoordinate coordinate) const {
        return FPoint(
                    FReal(coordinate.getX()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getX(),
                    FReal(coordinate.getY()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getY(),
                    FReal(coordinate.getZ()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getZ());
    }


    /** Return the combine of a paire of number */    
    FReal combin(const int& a, const int& b){
        if(a-b<0) printf("Error combi negative!! a=%d b=%d\n",a,b);
        return factorials[a] / (factorials[b]*factorials[a-b]);
    }

    /** To access the rotation matrix value with an analytical computation on the fly
      * These analytical formula has been taken from two papers:
      * Fast and accurate determination of the wigner rotation matrices in the fast multipole method
      * Formulas 13, 14, 15, 16
      * Parallelization of the Fast Multipole Method
      * Formula 19
      *  \f[
      * P^{l}_{m,k} = \frac{1}{2} \sqrt{ \frac{(l-m)!(l+m)!}{(l-k)!(l+k)!} } (1+sgn(k)cos( \theta))^{|k|} sin( \theta )^{m-|k|}
      *              \sum_{n=max(-(m+k),0)}^{min(l-m,l-k)}{ (-1)^{l-m-n} {{l-k} \choose {n}} {{l+k} \choose {l-m-n}}
      *              (1+cos( \theta ))^n (1+cos( \theta ))^{l-m-n} } , l \geq 0, -l \leq k \leq l, |k| \leq m \leq l
      * P^{l}_{m,k} = (-1)^{m+k} P^{l}_{k,m} , l \ge 0, -l \leq m \le 0, |m| \leq k \leq l
      * P^{l}_{m,k} = (-1)^{m+k} P^{l}_{k,m} , l \ge 0, 0 \leq m \le l, m \le k \leq l
      * \f]
      */
    FReal d_lmk_analytical(const FReal cosTheta, const FReal sinTheta, const int l , const int m, const int k){
        if( l >= 0 && -l <= k && k <= l && FMath::Abs(k) <= m && m <= l ){
            FReal sum = 0;
            for(int n = FMath::Max(-(m+k),0) ; n <= FMath::Min(l-m,l-k) ; ++n){
                sum += FMath::pow(FReal(-1.0), l-m-n) * combin(l-k, n) * combin(l+k, l-m-n) * FMath::pow(FReal(1.0)+cosTheta,n) * FMath::pow(FReal(1.0)-cosTheta, l-m-n);
            }
            return (FReal(1.0)/FMath::pow(FReal(2.0),l)) * FMath::Sqrt((factorials[l-m]*factorials[l+m])/(factorials[l-k]*factorials[l+k]))
                    * FMath::pow(FReal(1.0) + FReal(Sgn(k))*cosTheta, FMath::Abs(k)) * FMath::pow(sinTheta, m - FMath::Abs(k)) * sum;
        }
        // P^{l}_{m,k} = (-1)^{m+k} P^{l}_{k,m} , l \ge 0, -l \leq m \le 0, |m| \leq k \leq l
        else if( (l > 0 && -l <= m && m < 0 && FMath::Abs(m) <= k && k <= l)
                 ||  (l > 0 && 0 <= m && m < l && m < k && k <= l )) {
            return FMath::pow(FReal(-1.0), m+k) * d_lmk_analytical(cosTheta,sinTheta, l, k ,m);
        }
        // P^{l}_{m,k} = (-1)^{m+k} P^{l}_{k,m} , l \ge 0, 0 \leq m \le l, m \le k \leq l
        else if( l > 0 && -l <= m && m < l && -l <= k && k < -m ){
            return FMath::pow(FReal(-1.0), m+k) * d_lmk_analytical(cosTheta,sinTheta, l, -m, -k);
        }
        else{
            printf("Error that should not be possible!\n");
            return FReal(0.0);
        }
    }

    ///////////////////////////////////////////////////////
    // legendre
    ///////////////////////////////////////////////////////

    /** Compute the legendre polynomial from {0,0} to {P,P}
      * the computation is made by recurence (P cannot be equal to 0)
      * @todo use pointer
      *
      * The formula has been taken from:
      * Fast and accurate determination of the wigner rotation matrices in the fast multipole method
      * Formula number (22)
      * \f[
      * P_{0,0} = 1
      * P_{l,l} = (2l-1) sin( \theta ) P_{l-1,l-1} ,l \ge 0
      * P_{l,l-1} = (2l-1) cos( \theta ) P_{l-1,l-1} ,l \ge 0
      * P_{l,m} = \frac{(2l-1) cos( \theta ) P_{l-1,m} - (l+m-1) P_{l-2,m}x}{(l-k)} ,l \ge 1, 0 \leq m \le l-1
      * \f]
      */
    void computeLegendre(FReal legendre[], const FReal inCosTheta, const FReal inSinTheta){
        const FReal invSinTheta = -inSinTheta;

        legendre[0] = 1.0;             // P_0,0(1) = 1

        legendre[1] = inCosTheta;      // P_1,0 = cos(theta)
        legendre[2] = invSinTheta;     // P_1,1 = -sin(theta)

        // work with pointers
        FReal* FRestrict legendre_l1_m1 = legendre;     // P{l-2,m} starts with P_{0,0}
        FReal* FRestrict legendre_l1_m  = legendre + 1; // P{l-1,m} starts with P_{1,0}
        FReal* FRestrict legendre_lm  = legendre + 3;   // P{l,m} starts with P_{2,0}

        // Compute using recurrence
        FReal l2_minus_1 = 3; // 2 * l - 1
        FReal fl = FReal(2.0);// To get 'l' as a float
        for(int l = 2; l <= P ; ++l, ++fl ){
            FReal lm_minus_1 = fl - FReal(1.0); // l + m - 1
            FReal l_minus_m = fl;               // l - m
            for( int m = 0; m < l - 1 ; ++m ){
                // P_{l,m} = \frac{(2l-1) cos( \theta ) P_{l-1,m} - (l+m-1) P_{l-2,m}x}{(l-m)}
                *(legendre_lm++) = (l2_minus_1 * inCosTheta * (*legendre_l1_m++) - (lm_minus_1++) * (*legendre_l1_m1++) )
                                                        / (l_minus_m--);
            }
            // P_{l,l-1} = (2l-1) cos( \theta ) P_{l-1,l-1}
            *(legendre_lm++) = l2_minus_1 * inCosTheta * (*legendre_l1_m);
            // P_{l,l} = (2l-1) sin( \theta ) P_{l-1,l-1}
            *(legendre_lm++) = l2_minus_1 * invSinTheta * (*legendre_l1_m);
            // goto P_{l-1,0}
            ++legendre_l1_m;
            l2_minus_1 += FReal(2.0); // 2 * l - 1 => progress by two
        }
    }

    ///////////////////////////////////////////////////////
    // Real rotation
    ///////////////////////////////////////////////////////


    /** This function rotate a multipole vector by an angle inclination \theta
      * The formula used is present in several paper, but we refer to
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1
      * \f[
      * O_{l,m}( \alpha + \theta, \beta ) = \sum_{k=-l}^l{ \sqrt{ \frac{(l-k)!(l+k)!}{(l-|m|)!(l+|m|)!} }
      *                                     d^l_{km}( \theta ) O_{l,k}( \alpha, \beta ) }
      * \f]
      * Because we store only P_lm for l >= 0 and m >= 0 we use the relation of symetrie as:
      * \f$ O_{l,-m} = \bar{ O_{l,m} } (-1)^m \f$
      * Theta should be in [0,pi]
      */
    void rotateMultipoleAroundY(FComplexe vec[], const FReal theta){
        FComplexe cell_rotate[SizeArray];
        for(int l = 0 ; l <= P ; ++l){
            for(int m = 0 ; m <= l ; ++m ){
                FReal w_lkm_real = 0.0;
                FReal w_lkm_imag = 0.0;

                for(int k = -l ; k < 0 ; ++k){
                    // O_{l,-m} = \bar{ O_{l,m} } (-1)^m
                    const FReal d_lmk = d_lmk_analytical(FMath::Cos(theta),FMath::Sin(theta),l,m,k);
                    // \sqrt{ \frac{(l-k)!(l+k)!}{(l-|m|)!(l+|m|)!} }
                    const FReal factor = FMath::Sqrt((factorials[l-k]*factorials[l+k])/(factorials[l-abs(m)]*factorials[l+abs(m)]));
                    w_lkm_real += FMath::pow(FReal(-1.0),-k) * factor * d_lmk * vec[atLm(l,-k)].getReal(); // k<0 => Conjugate * -1^k
                    w_lkm_imag -= FMath::pow(FReal(-1.0),-k) * factor * d_lmk * vec[atLm(l,-k)].getImag(); // k<0 => Conjugate * -1^k
                }
                for(int k = 0 ; k <= l ; ++k){
                    const FReal d_lmk = d_lmk_analytical(FMath::Cos(theta),FMath::Sin(theta),l,m,k);
                    // \sqrt{ \frac{(l-k)!(l+k)!}{(l-|m|)!(l+|m|)!} }
                    const FReal factor = FMath::Sqrt((factorials[l-k]*factorials[l+k])/(factorials[l-abs(m)]*factorials[l+abs(m)]));
                    w_lkm_real += factor * d_lmk * vec[atLm(l,k)].getReal();
                    w_lkm_imag += factor * d_lmk * vec[atLm(l,k)].getImag();
                }
                cell_rotate[atLm(l,m)].setRealImag(w_lkm_real, w_lkm_imag);
            }
        }
        FMemUtils::copyall(vec,cell_rotate,SizeArray);
    }

    /** This function rotate a local vector by an angle inclination \theta
      * The formula used is present in several paper, but we refer to
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1
      * \f[
      * M_{l,m}( \alpha + \theta, \beta ) = \sum_{k=-l}^l{ \sqrt{ \frac{(l-|m|)!(l+|m|)!}{(l-k)!(l+k)!} }
      *                                     d^l_{km}( \theta ) M_{l,k}( \alpha, \beta ) }
      * \f]
      * Because we store only P_lm for l >= 0 and m >= 0 we use the relation of symetrie as:
      * \f$ M_{l,-m} = \bar{ M_{l,m} } (-1)^m \f$
      * Theta should be in [0,pi]
      */
    void rotateTaylorAroundY(FComplexe vec[], const FReal theta){
        FComplexe cell_rotate[SizeArray];
        for(int l = 0 ; l <= P ; ++l){
            for(int m = 0 ; m <= l ; ++m ){
                FReal w_lkm_real = 0.0;
                FReal w_lkm_imag = 0.0;

                for(int k = -l ; k < 0 ; ++k){
                    // M_{l,-m} = \bar{ M_{l,m} } (-1)^m
                    const FReal d_lmk = d_lmk_analytical(FMath::Cos(theta),FMath::Sin(theta),l,m,k);
                    // \sqrt{ \frac{(l-|m|)!(l+|m|)!}{(l-k)!(l+k)!} }
                    const FReal factor = FMath::Sqrt((factorials[l-abs(m)]*factorials[l+abs(m)])/(factorials[l-k]*factorials[l+k]));
                    w_lkm_real += FMath::pow(FReal(-1.0),-k) * factor * d_lmk * vec[atLm(l,-k)].getReal(); // k<0 => Conjugate * -1^k
                    w_lkm_imag -= FMath::pow(FReal(-1.0),-k) * factor * d_lmk * vec[atLm(l,-k)].getImag(); // k<0 => Conjugate * -1^k
                }
                for(int k = 0 ; k <= l ; ++k){
                    const FReal d_lmk = d_lmk_analytical(FMath::Cos(theta),FMath::Sin(theta),l,m,k);
                    // \sqrt{ \frac{(l-|m|)!(l+|m|)!}{(l-k)!(l+k)!} }
                    const FReal factor = FMath::Sqrt((factorials[l-abs(m)]*factorials[l+abs(m)])/(factorials[l-k]*factorials[l+k]));
                    w_lkm_real += factor * d_lmk * vec[atLm(l,k)].getReal();
                    w_lkm_imag += factor * d_lmk * vec[atLm(l,k)].getImag();
                }
                cell_rotate[atLm(l,m)].setRealImag(w_lkm_real, w_lkm_imag);
            }
        }
        FMemUtils::copyall(vec,cell_rotate,SizeArray);
    }

    /** This function rotate a multipole vector by a angles inclination & azimuth
      * The formula used is present in several paper, but we refer to
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1 as the forward rotation
      *
      * Rotation are not commutative so we have to do it in the right order
      */
    void rotateMultipole(FComplexe vec[], const FReal azimuth, const FReal inclination){
        //rotateMultipoleAroundZ(vec,(FMath::FPiDiv2 + azimuth));
        rotateMultipoleAroundY(vec,inclination);
    }
    /** This function rotate a multipole vector by a angles inclination & azimuth
      * The formula used is present in several paper, but we refer to
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1 as the forward rotation
      *
      * Rotation are not commutative so we have to do it in the right order
      */
    void deRotateMultipole(FComplexe vec[], const FReal azimuth, const FReal inclination){
        rotateMultipoleAroundY(vec,-inclination);
        //rotateMultipoleAroundZ(vec,-(FMath::FPiDiv2 + azimuth));
    }

    /** This function rotate a local vector by a angles inclination & azimuth
      * The formula used is present in several paper, but we refer to
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1 as the forward rotation
      *
      * Rotation are not commutative so we have to do it in the right order
      */
    void rotateTaylor(FComplexe vec[], const FReal azimuth, const FReal inclination){
        //rotateTaylorAroundZ(vec,(FMath::FPiDiv2 + azimuth));
        rotateTaylorAroundY(vec,inclination);
    }
    /** This function rotate a local vector by a angles inclination & azimuth
      * The formula used is present in several paper, but we refer to
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1 as the forward rotation
      *
      * Rotation are not commutative so we have to do it in the right order
      */
    void deRotateTaylor(FComplexe vec[], const FReal azimuth, const FReal inclination){
        rotateTaylorAroundY(vec,-inclination);
        //rotateTaylorAroundZ(vec,-(FMath::FPiDiv2 + azimuth));
    }

    static void ComplexeArrayMulEqual(FComplexe dest[], const FComplexe src[], const int sizeOfArray){
        for(int idx = 0 ; idx < sizeOfArray ; ++idx) {
            dest[idx] *= src[idx];
        }
    }

public:

    /** Constructor, needs system information */
    FRotationKernel( const int inTreeHeight, const FReal inBoxWidth, const FPoint& inBoxCenter) :
        boxWidth(inBoxWidth),
        treeHeight(inTreeHeight),
        widthAtLeafLevel(inBoxWidth/FReal(1 << (inTreeHeight-1))),
        widthAtLeafLevelDiv2(widthAtLeafLevel/2),
        boxCorner(inBoxCenter.getX()-(inBoxWidth/2),inBoxCenter.getY()-(inBoxWidth/2),inBoxCenter.getZ()-(inBoxWidth/2))
        {
        // simply does the precomputation
        precomputeFactorials();
        preComputePosition();
        precomputeTranslationCoef();
        precomputeRotationVectors();
    }

    /** Default destructor */
    virtual ~FRotationKernel(){
    }

    /** P2M
      * The computation is based on the paper :
      * Parallelization of the fast multipole method
      * Formula number 10, page 3
      * \f[
      * \omega (q,a) = q \frac{a^{l}}{(l+|m|)!} P_{lm}(cos( \alpha ) )e^{-im \beta}
      * \f]
      */
    void P2M(CellClass* const inPole, const ContainerClass* const inParticles ) {
        const FReal i_pow_m[4] = {0, FMath::FPiDiv2, FMath::FPi, -FMath::FPiDiv2};
        // w is the multipole moment
        FComplexe* FRestrict const w = inPole->getMultipole();

        // Copying the position is faster than using cell position
        const FPoint cellPosition = getLeafCenter(inPole->getCoordinate());

        // We need a legendre array
        FReal legendre[SizeArray];

        // For all particles in the leaf box
        typename ContainerClass::ConstBasicIterator iterParticle(*inParticles);
        while( iterParticle.hasNotFinished()){
            // P2M
            const ParticleClass& particle = iterParticle.data();
            const FSpherical sph(particle.getPosition() - cellPosition);

            // The physical value (charge, mass)
            const FReal q = particle.getPhysicalValue();
            // The distance between the SH and the particle
            const FReal a = sph.getR();

            // Compute the legendre polynomial
            computeLegendre(legendre, sph.getCosTheta(), sph.getSinTheta());

            // w{l,m}(q,a) = q a^l/(l+|m|)! P{l,m}(cos(alpha)) exp(-i m Beta)
            FReal q_aPowL = q; // To consutrct q*a^l continously
            int index_l_m = 0; // To construct the index of (l,m) continously
            for(int l = 0 ; l <= P ; ++l ){
                FReal fm = 0.0; // To have "m" has a float
                for(int m = 0 ; m <= l ; ++m, ++index_l_m, ++fm){
                    const FReal magnitude = q_aPowL * legendre[index_l_m] / factorials[l+m];
                    w[index_l_m].incReal(magnitude * FMath::Cos(fm * sph.getPhi() + i_pow_m[m & 0x3]));
                    w[index_l_m].incImag(magnitude * FMath::Sin(fm * sph.getPhi() + i_pow_m[m & 0x3]));
                }
                q_aPowL *= a;
            }

            // Goto next particle
            iterParticle.gotoNext();
        }
    }

    /** M2M
      * The operator A has been taken from :
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1 as the operator A
      * \f[
      * O_{l,m}(a+b') = \sum_{j=|m|}^l{ \frac{ b^{l-j} }{ (l-j)! } O_{j,m}(a) }
      * \f]
      * As describe in the paper, when need first to rotate the SH
      * then transfer using the formula
      * and finaly rotate back.
      */
    void M2M(CellClass* const FRestrict inPole, const CellClass*const FRestrict *const FRestrict inChildren, const int inLevel) {
        // Get the translation coef for this level (same for all chidl)
        const FReal*const coef = M2MTranslationCoef[inLevel];
        // A buffer to copy the source w allocated once
        FComplexe source_w[SizeArray];
        // For all children
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            // if child exists
            if(inChildren[idxChild]){
                // Copy the source
                FMemUtils::copyall(source_w, inChildren[idxChild]->getMultipole(), SizeArray);

                // rotate it forward
                const FSpherical sph = getSphericalChild(inLevel, idxChild);
                ComplexeArrayMulEqual(source_w,rotationExpMinusImPhi[idxChild],SizeArray);
                rotateMultipole(source_w, sph.getAzimuth(), sph.getInclination());

                //const FReal b = -sph.getR();
                // Translate it
                FComplexe target_w[SizeArray];
                int index_lm = 0;
                for(int l = 0 ; l <= P ; ++l ){
                    for(int m = 0 ; m <= l ; ++m, ++index_lm ){
                        // w{l,m}(a+b) = sum(j=m:l, b^(l-j)/(l-j)! w{j,m}(a)
                        FReal w_lm_real = 0.0;
                        FReal w_lm_imag = 0.0;
                        int index_jm = atLm(m,m);   // get atLm(l,m)
                        int index_l_minus_j = l-m;  // get l-j continously
                        for(int j = m ; j <= l ; ++j, --index_l_minus_j, index_jm += j ){
                            //const coef = (b^l-j) / (l-j)!;
                            w_lm_real += coef[index_l_minus_j] * source_w[index_jm].getReal();
                            w_lm_imag += coef[index_l_minus_j] * source_w[index_jm].getImag();
                        }
                        target_w[index_lm].setRealImag(w_lm_real,w_lm_imag);
                    }
                }

                // Rotate it back
                deRotateMultipole(target_w, sph.getAzimuth(), sph.getInclination());
                ComplexeArrayMulEqual(target_w,rotationExpImPhi[idxChild],SizeArray);

                // Sum the result
                FMemUtils::addall( inPole->getMultipole(), target_w, SizeArray);
            }
        }
    }

    /** M2L
      * The operator B has been taken from :
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1 as the operator B
      * \f[
      * M_{l,m}(a-b') = \sum_{j=|m|}^{\infty}{ \frac{ (j+l)! } { b^{j+l+1} } O_{j,-m}(a) } , \textrm{j bounded by P-l}
      * \f]
      * As describe in the paper, when need first to rotate the SH
      * then transfer using the formula
      * and finaly rotate back.
      */
    void M2L(CellClass* const FRestrict inLocal, const CellClass* inInteractions[], const int /*inSize*/, const int inLevel) {
        // To copy the multipole data allocated once
        FComplexe source_w[SizeArray];
        // For all children
        for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
            // if interaction exits
            if(inInteractions[idxNeigh]){
                const FReal*const coef = M2LTranslationCoef[inLevel][idxNeigh];
                // Copy multipole data into buffer
                FMemUtils::copyall(source_w, inInteractions[idxNeigh]->getMultipole(), SizeArray);

                // Rotate
                const FSpherical sph = getSphericalInteraction(inLevel, idxNeigh);
                rotateMultipole(source_w, sph.getAzimuth(), sph.getInclination());

                // Transfer to u
                FComplexe target_u[SizeArray];
                int index_lm = 0;
                for(int l = 0 ; l <= P ; ++l ){
                    FReal minus_1_pow_m = 1.0;
                    for(int m = 0 ; m <= l ; ++m, ++index_lm ){
                        // u{l,m}(a-b) = sum(j=|m|:P-l, (j+l)!/b^(j+l+1) w{j,-m}(a)
                        FReal u_lm_real = 0.0;
                        FReal u_lm_imag = 0.0;
                        int index_jl = m + l;       // get j+l
                        int index_jm = atLm(m,m);   // get atLm(l,m)
                        for(int j = m ; j <= P-l ; ++j, ++index_jl, index_jm += j ){
                            // coef = (j+l)!/b^(j+l+1)
                            // because {l,-m} => {l,m} conjugate -1^m with -i
                            u_lm_real += minus_1_pow_m * coef[index_jl] * source_w[index_jm].getReal();
                            u_lm_imag -= minus_1_pow_m * coef[index_jl] * source_w[index_jm].getImag();
                        }
                        target_u[index_lm].setRealImag(u_lm_real,u_lm_imag);
                        minus_1_pow_m = -minus_1_pow_m;
                    }
                }

                // Rotate it back
                deRotateTaylor(target_u, sph.getAzimuth(), sph.getInclination());

                // Sum
                FMemUtils::addall(inLocal->getLocal(), target_u, SizeArray);
            }
        }
    }

    /** L2L
      * The operator C has been taken from :
      * Implementation of rotation-based operators for Fast Multipole Method in X10
      * At page 5 .1 as the operator C
      * \f[
      * M_{l,m}(a-b') = \sum_{j=l}^{\infty}{ \frac{ b^{j-l} }{ (j-l)! } M_{j,m}(a) } , \textrm{j bounded by P}
      * \f]
      * As describe in the paper, when need first to rotate the SH
      * then transfer using the formula
      * and finaly rotate back.
      */
    void L2L(const CellClass* const FRestrict inLocal, CellClass* FRestrict *const FRestrict  inChildren, const int inLevel) {
        // Get the translation coef for this level (same for all chidl)
        const FReal*const coef = L2LTranslationCoef[inLevel];
        // To copy the source local allocated once
        FComplexe source_u[SizeArray];
        // For all children
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            // if child exists
            if(inChildren[idxChild]){
                // Copy the local data into the buffer
                FMemUtils::copyall(source_u, inLocal->getLocal(), SizeArray);

                // Rotate
                const FSpherical sph = getSphericalChild(inLevel, idxChild);
                ComplexeArrayMulEqual(source_u,rotationExpImPhi[idxChild],SizeArray);
                rotateTaylor(source_u, sph.getAzimuth(), sph.getInclination());

                // Translate
                FComplexe target_u[SizeArray];
                for(int l = 0 ; l <= P ; ++l ){
                    for(int m = 0 ; m <= l ; ++m ){
                        // u{l,m}(r-b) = sum(j=0:P, b^(j-l)/(j-l)! u{j,m}(r);
                        FReal u_lm_real = 0.0;
                        FReal u_lm_imag = 0.0;
                        int index_jm = atLm(l,m);   // get atLm(j,m)
                        int index_j_minus_l = 0;    // get l-j continously
                        for(int j = l ; j <= P ; ++j, ++index_j_minus_l, index_jm += j){
                            // coef = b^j-l/j-l!
                            u_lm_real += coef[index_j_minus_l] * source_u[index_jm].getReal();
                            u_lm_imag += coef[index_j_minus_l] * source_u[index_jm].getImag();
                        }
                        target_u[atLm(l,m)].setRealImag(u_lm_real,u_lm_imag);
                    }
                }

                // Rotate
                deRotateTaylor(target_u, sph.getAzimuth(), sph.getInclination());
                ComplexeArrayMulEqual(source_u,rotationExpMinusImPhi[idxChild],SizeArray);

                // Sum in child
                FMemUtils::addall(inChildren[idxChild]->getLocal(), target_u, SizeArray);
            }
        }
    }

    /** L2P
      * Equation are coming from the PhD report of Pierre Fortin.
      * We have two different computations, one for the potential (end of function)
      * the other for the forces.
      *
      * The potential use the fallowing formula, page 36, formula 2.14 + 1:
      * \f[
      *  \Phi = \sum_{j=0}^P{\left( u_{j,0} I_{j,0}(r, \theta, \phi) + \sum_{k=1}^j{2 Re(u_{j,k} I_{j,k}(r, \theta, \phi))} \right)},
      *  \textrm{since } u_{l,-m} = (-1)^m \overline{ u_{l,m} }
      * \f]
      *
      * The forces are coming form the formulas, page 37, formulas 2.14 + 3:
      * \f[
      * F_r = -\frac{1}{r} \left( \sum_{j=1}^P{j u_{j,0} I_{j,0}(r, \theta, \phi) } + \sum_{k=1}^j{2 j Re(u_{j,k} I_{j,k}(r, \theta, \phi))} \right)
      * F_{ \theta } = -\frac{1}{r} \left( \sum_{j=0}^P{j u_{j,0} \frac{ \partial I_{j,0}(r, \theta, \phi) }{ \partial \theta } } + \sum_{k=1}^j{2 Re(u_{j,k} \frac{ \partial I_{j,k}(r, \theta, \phi) }{ \partial \theta })} \right)
      * F_{ \phi } = -\frac{1}{r sin \phi} \sum_{j=0}^P \sum_{k=1}^j{(-2k) Im(u_{j,k} I_{j,k}(r, \theta, \phi)) }}
      * \f]
      */
    void L2P(const CellClass* const inLocal, ContainerClass* const inParticles){
        const FReal i_pow_m[4] = {0, FMath::FPiDiv2, FMath::FPi, -FMath::FPiDiv2};
        // Take the local value from the cell
        const FComplexe* FRestrict const u = inLocal->getLocal();

        // Copying the position is faster than using cell position
        const FPoint cellPosition = getLeafCenter(inLocal->getCoordinate());

        // For all particles in the leaf box
        typename ContainerClass::BasicIterator iterParticle(*inParticles);
        while( iterParticle.hasNotFinished()){
            // L2P
            ParticleClass& particle = iterParticle.data();
            const FSpherical sph(particle.getPosition() - cellPosition);

            // The distance between the SH and the particle
            const FReal r = sph.getR();

            // Compute the legendre polynomial
            FReal legendre[SizeArray];
            computeLegendre(legendre, sph.getCosTheta(), sph.getSinTheta());

            // pre compute what is used more than once
            FReal minus_r_pow_l_div_fact_lm[SizeArray];
            FReal minus_r_pow_l_legendre_div_fact_lm[SizeArray];
            {
                int index_lm = 0;
                FReal minus_r_pow_l = 1.0;  // To get (-1*r)^l
                for(int l = 0 ; l <= P ; ++l){
                    for(int m = 0 ; m <= l ; ++m, ++index_lm){
                        minus_r_pow_l_div_fact_lm[index_lm] = minus_r_pow_l / factorials[l+m];
                        minus_r_pow_l_legendre_div_fact_lm[index_lm] = minus_r_pow_l_div_fact_lm[index_lm] * legendre[index_lm];
                    }
                    minus_r_pow_l *= -r;
                }
            }
            // pre compute what is use more than once
            FReal cos_m_phi_i_pow_m[P+1];
            FReal sin_m_phi_i_pow_m[P+1];
            {
                for(int m = 0 ; m <= P ; ++m){
                    const FReal m_phi_i_pow_m = FReal(m) * sph.getPhi() + i_pow_m[m & 0x3];
                    cos_m_phi_i_pow_m[m] = FMath::Cos(m_phi_i_pow_m);
                    sin_m_phi_i_pow_m[m] = FMath::Sin(m_phi_i_pow_m);
                }
            }

            // compute the forces
            {
                FReal Fr = 0;
                FReal FO = 0;
                FReal Fp = 0;

                int index_lm = 1;          // To get atLm(l,m), warning starts with l = 1
                FReal fl = 1.0;            // To get "l" as a float

                for(int l = 1 ; l <= P ; ++l, ++fl){
                    // first m == 0
                    {
                        Fr += fl * u[index_lm].getReal() * minus_r_pow_l_legendre_div_fact_lm[index_lm];
                    }
                    {
                        const FReal coef = minus_r_pow_l_div_fact_lm[index_lm] * (fl * (sph.getCosTheta()*legendre[index_lm]
                                                  - legendre[index_lm-l]) / sph.getSinTheta());
                        const FReal dI_real = coef;
                        // F(O) += 2 * Real(L dI/dO)
                        FO += u[index_lm].getReal() * dI_real;
                    }
                    ++index_lm;
                    // then 0 < m
                    for(int m = 1 ; m <= l ; ++m, ++index_lm){
                        {
                            const FReal coef = minus_r_pow_l_legendre_div_fact_lm[index_lm];
                            const FReal I_real = coef * cos_m_phi_i_pow_m[m];
                            const FReal I_imag = coef * sin_m_phi_i_pow_m[m];
                            // F(r) += 2 x l x Real(LI)
                            Fr += 2 * fl * (u[index_lm].getReal() * I_real - u[index_lm].getImag() * I_imag);
                            // F(p) += -2 x m x Imag(LI)
                            Fp -= 2 * FReal(m) * (u[index_lm].getReal() * I_imag + u[index_lm].getImag() * I_real);
                        }
                        {
                            const FReal legendre_l_minus_1 = (m == l) ? FReal(0.0) : FReal(l+m)*legendre[index_lm-l];
                            const FReal coef = minus_r_pow_l_div_fact_lm[index_lm] * ((fl * sph.getCosTheta()*legendre[index_lm]
                                                      - legendre_l_minus_1) / sph.getSinTheta());
                            const FReal dI_real = coef * cos_m_phi_i_pow_m[m];
                            const FReal dI_imag = coef * sin_m_phi_i_pow_m[m];
                            // F(O) += 2 * Real(L dI/dO)
                            FO += FReal(2.0) * (u[index_lm].getReal() * dI_real - u[index_lm].getImag() * dI_imag);
                        }
                    }
                }
                // div by r
                Fr /= sph.getR();
                FO /= sph.getR();
                Fp /= sph.getR() * sph.getSinTheta();

                // copy variable from spherical position
                const FReal cosPhi     = FMath::Cos(sph.getPhi());
                const FReal sinPhi     = FMath::Sin(sph.getPhi());
                const FReal physicalValue = particle.getPhysicalValue();

                // compute forces
                const FReal forceX = (
                        cosPhi * sph.getSinTheta() * Fr  +
                        cosPhi * sph.getCosTheta() * FO +
                        (-sinPhi) * Fp) * physicalValue;

                const FReal forceY = (
                        sinPhi * sph.getSinTheta() * Fr  +
                        sinPhi * sph.getCosTheta() * FO +
                        cosPhi * Fp) * physicalValue;

                const FReal forceZ = (
                        sph.getCosTheta() * Fr +
                        (-sph.getSinTheta()) * FO) * physicalValue;

                // inc particles forces
                particle.incForces( forceX, forceY, forceZ );
            }
            // compute the potential
            {
                FReal magnitude = 0;
                // E = sum( l = 0:P, sum(m = -l:l, u{l,m} ))
                int index_lm = 0;
                for(int l = 0 ; l <= P ; ++l ){
                    {//for m == 0
                        // (l-|m|)! * P{l,0} / r^(l+1)
                        magnitude += u[index_lm].getReal() * minus_r_pow_l_legendre_div_fact_lm[index_lm];
                        ++index_lm;
                    }
                    for(int m = 1 ; m <= l ; ++m, ++index_lm ){
                        const FReal coef = minus_r_pow_l_legendre_div_fact_lm[index_lm];
                        const FReal I_real = coef * cos_m_phi_i_pow_m[m];
                        const FReal I_imag = coef * sin_m_phi_i_pow_m[m];
                        magnitude += FReal(2.0) * ( u[index_lm].getReal() * I_real - u[index_lm].getImag() * I_imag );
                    }
                }
                // inc potential
                particle.incPotential(magnitude);
            }
            // progress
            iterParticle.gotoNext();
        }
    }


    /** P2P
      * This function proceed the P2P using particlesMutualInteraction
      * The computation is done for interactions with an index <= 13.
      * (13 means current leaf (x;y;z) = (0;0;0)).
      * Calling this method in multi thread should be done carrefully.
      */
    void P2P(const FTreeCoordinate& /*inPosition*/,
                     ContainerClass* const FRestrict inTargets, const ContainerClass* const FRestrict /*inSources*/,
                     ContainerClass* const inNeighbors[27], const int /*inSize*/){

        {
            typename ContainerClass::BasicIterator iterTarget(*inTargets);
            while( iterTarget.hasNotFinished() ){
                // We copy the target particle to work with a particle in the heap
                ParticleClass target( iterTarget.data() );

                // For all particles after the current one
                typename ContainerClass::BasicIterator iterSameBox = iterTarget;
                iterSameBox.gotoNext();
                while( iterSameBox.hasNotFinished() ){
                    particlesMutualInteraction(&target, &iterSameBox.data());
                    iterSameBox.gotoNext();
                }
                // Set data and progress
                iterTarget.setData(target);
                iterTarget.gotoNext();
            }
        }
        // For all the neigbors leaves
        for(int idxDirectNeighbors = 0 ; idxDirectNeighbors <= 13 ; ++idxDirectNeighbors){
            if( inNeighbors[idxDirectNeighbors] ){
                // For all particles in current leaf
                typename ContainerClass::BasicIterator iterTarget(*inTargets);
                while( iterTarget.hasNotFinished() ){
                    ParticleClass target( iterTarget.data() );
                    // For all the particles in the other leaf
                    typename ContainerClass::BasicIterator iterSource(*inNeighbors[idxDirectNeighbors]);
                    while( iterSource.hasNotFinished() ){
                        particlesMutualInteraction(&target, &iterSource.data());
                        iterSource.gotoNext();
                    }
                    // Set data and progress
                    iterTarget.setData(target);
                    iterTarget.gotoNext();
                }
            }
        }
    }


    /** Use mutual even if it not useful and call particlesMutualInteraction */
    void P2PRemote(const FTreeCoordinate& inPosition,
                   ContainerClass* const FRestrict inTargets, const ContainerClass* const FRestrict inSources,
                   ContainerClass* const inNeighbors[27], const int inSize){
        for(int idxDirectNeighbors = 0 ; idxDirectNeighbors < 27 ; ++idxDirectNeighbors){
            if( inNeighbors[idxDirectNeighbors] ){
                // For all particles in current leaf
                typename ContainerClass::BasicIterator iterTarget(*inTargets);
                while( iterTarget.hasNotFinished() ){
                    ParticleClass target( iterTarget.data() );
                    // For all the particles in the other leaf
                    typename ContainerClass::BasicIterator iterSource(*inNeighbors[idxDirectNeighbors]);
                    while( iterSource.hasNotFinished() ){
                        particlesMutualInteraction(&target, &iterSource.data());
                        iterSource.gotoNext();
                    }
                    // Set data and progress
                    iterTarget.setData(target);
                    iterTarget.gotoNext();
                }
            }
        }
    }

    /** P2P mutual interaction,
      * this function computes the interaction for 2 particles.
      *
      * Formulas are:
      * \f[
      * F = q_1 * q_2 / r^2
      * P_1 = q_2 / r ; P_2 = q_1 / r
      * \f]
      * In details :
      * \f$ F(x) = \frac{ \Delta_x * q_1 * q_2 }{ r^2 } = \Delta_x * F \f$
      */
    void particlesMutualInteraction(ParticleClass*const FRestrict target, ParticleClass*const FRestrict source) const {

        FReal dx = source->getPosition().getX() - target->getPosition().getX();
        FReal dy = source->getPosition().getY() - target->getPosition().getY();
        FReal dz = source->getPosition().getZ() - target->getPosition().getZ();

        FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
        FReal inv_distance = FMath::Sqrt(inv_square_distance);

        inv_square_distance *= inv_distance;
        inv_square_distance *= target->getPhysicalValue() * source->getPhysicalValue();

        dx *= inv_square_distance;
        dy *= inv_square_distance;
        dz *= inv_square_distance;

        target->incForces( dx, dy, dz);
        target->incPotential( inv_distance * source->getPhysicalValue() );

        source->incForces( (-dx), (-dy), (-dz));
        source->incPotential( inv_distance * target->getPhysicalValue() );

    }
};


#endif // FROTATIONKERNEL_HPP
