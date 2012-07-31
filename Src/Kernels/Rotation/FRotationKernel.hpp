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

#include "../../Components/FAbstractSerializable.hpp"
#include "../../Components/FAbstractSendable.hpp"
#include "../../Utils/FComplexe.hpp"
#include "../../Utils/FMemUtils.hpp"

#include "../../Extensions/FExtendCellType.hpp"

#include "../../Components/FBasicCell.hpp"

#include "../../Containers/FBufferWriter.hpp"
#include "../../Containers/FBufferReader.hpp"

template <int P>
class FRotationCell : public FBasicCell {
protected:
    static const int MultipoleSize = ((P+2)*(P+1))/2; // Artimethique suite (n+1)*n/2
    static const int LocalSize = ((P+2)*(P+1))/2;     // Artimethique suite (n+1)*n/2

    FComplexe multipole_exp[MultipoleSize]; //< For multipole extenssion
    FComplexe local_exp[LocalSize];         //< For local extenssion

public:
    static int GetLocalSize(){
        return LocalSize;
    }

    static int GetPoleSize(){
        return MultipoleSize;
    }

    /** Default constructor */
    FRotationCell(){
    }

    /** Constructor */
    FRotationCell(const FRotationCell& other){
        (*this) = other;
    }

    /** Default destructor */
    virtual ~FRotationCell(){
    }

    /** Copy constructor */
    FRotationCell& operator=(const FRotationCell& other) {
        FMemUtils::copyall(multipole_exp, other.multipole_exp, MultipoleSize);
        FMemUtils::copyall(local_exp, other.local_exp, LocalSize);
        return *this;
    }

    /** Get Multipole */
    const FComplexe* getMultipole() const {
        return multipole_exp;
    }
    /** Get Local */
    const FComplexe* getLocal() const {
        return local_exp;
    }

    /** Get Multipole */
    FComplexe* getMultipole() {
        return multipole_exp;
    }
    /** Get Local */
    FComplexe* getLocal() {
        return local_exp;
    }

    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    void serializeUp(FBufferWriter& buffer) const{
        buffer.write(multipole_exp, MultipoleSize);
    }
    void deserializeUp(FBufferReader& buffer){
        buffer.fillArray(multipole_exp, MultipoleSize);
    }

    void serializeDown(FBufferWriter& buffer) const{
        buffer.write(local_exp, LocalSize);
    }
    void deserializeDown(FBufferReader& buffer){
        buffer.fillArray(local_exp, LocalSize);
    }

    ///////////////////////////////////////////////////////
    // to extend Serializable
    ///////////////////////////////////////////////////////
    void save(FBufferWriter& buffer) const{
        FBasicCell::save(buffer);
        buffer.write(multipole_exp, MultipoleSize);
        buffer.write(local_exp, LocalSize);
    }
    void restore(FBufferReader& buffer){
        FBasicCell::restore(buffer);
        buffer.fillArray(multipole_exp, MultipoleSize);
        buffer.fillArray(local_exp, LocalSize);
    }
};



#include "../../Components/FAbstractKernels.hpp"
#include "../../Utils/FSmartPointer.hpp"
#include "../../Utils/FComplexe.hpp"
#include "../../Utils/FSpherical.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FRotationKernel
* @brief
*
* This kernels is empty for now and does nothing.
*/
template< class ParticleClass, class CellClass, class ContainerClass, int P>
class FRotationKernel : public FAbstractKernels<ParticleClass,CellClass,ContainerClass> {    
    static const int SizeArray = ((P+2)*(P+1))/2;

    /** Return position in the array of the l/m couple */
    int atLm(const int l, const int m){
        // summation series over l + m => (l*(l+1))/2 + m
        return ((l*(l+1))>>1) + m;
    }

    template < class T >
    T Sgn(const T a){
        if(a < 0) return T(-1);
        else if(a > 0) return T(1);
        return T(0);
    }

    template < class T >
    T fact(const T a){
        if(a<0) printf("Error factorial negative!! a=%d\n",a);
        int result = 1;
        for(int i = 1; i<= a ; ++i){
            result *= i;
        }
        return T(result);
    }

    template < class T >
    FReal combin(const T& a, const T& b){
        if(a-b<0) printf("Error combi negative!! a=%d b=%d\n",a,b);
        return FReal(fact(a)) / FReal(fact(b)*fact(a-b));
    }

    int Min(const int a , const int b){
        return (a<b?a:b);
    }

    int Max(const int a , const int b){
        return (a>b?a:b);
    }

    double analytical(const double cosTheta, const double sinTheta, const int l , const int m, const int k){
        if( l >= 0 && -l <= k && k <= l && FMath::Abs(k) <= m && m <= l ){
            FReal sum = 0;
            for(int n = Max(-(m+k),0) ; n <= Min(l-m,l-k) ; ++n){
                sum += FMath::pow(-1.0, l-m-n) * combin(l-k, n) * combin(l+k, l-m-n) * FMath::pow(1.0+cosTheta,n) * FMath::pow(1.0-cosTheta, l-m-n);
            }
            return (1.0/FMath::pow(2.0,l)) * FMath::Sqrt(FReal(fact(l-m)*fact(l+m))/FReal(fact(l-k)*fact(l+k)))
                    * FMath::pow(1.0 + Sgn(k)*cosTheta, FMath::Abs(k)) * FMath::pow(sinTheta, m - FMath::Abs(k)) * sum;
        }
        else if( (l > 0 && -l <= m && m < 0 && FMath::Abs(m) <= k && k <= l)
                 ||  (l > 0 && 0 <= m && m < l && m < k && k <= l )) {
            return FMath::pow(-1.0, m+k) * analytical(cosTheta,sinTheta, l, k ,m);
        }
        else if( l > 0 && -l <= m && m < l && -l <= k && k < -m ){
            return FMath::pow(-1.0, m+k) * analytical(cosTheta,sinTheta, l, -m, -k);
        }
        else{
            return -999999.99999;
        }
    }



    /** @todo use pointer */
    void computeLegendre(FReal legendre[], const FReal inCosTheta, const FReal inSinTheta){
        const FReal invSinTheta = inSinTheta;
        legendre[atLm(0,0)] = 1.0;        // P_0,0(1) = 1

        legendre[atLm(1,0)] = inCosTheta;
        legendre[atLm(1,1)] = invSinTheta;

        for(int l = 2; l <= P ; ++l ){
            for( int m = 0; m < l - 1 ; ++m ){
                legendre[atLm(l,m)] = (FReal(2*l-1) * inCosTheta * legendre[atLm(l-1,m)] - FReal( l + m - 1 ) * legendre[atLm(l-2,m)] )
                                                                    / FReal( l - m );
            }
            legendre[atLm(l,l-1)] = FReal(2*l-1) * inCosTheta * legendre[atLm(l-1,l-1)];
            legendre[atLm(l,l)] = FReal(2*l-1) * invSinTheta * legendre[atLm(l-1,l-1)];
        }
    }


    FPoint getLeafCenter(const FTreeCoordinate coordinate) const {
        return FPoint(
                    FReal(coordinate.getX()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getX(),
                    FReal(coordinate.getY()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getY(),
                    FReal(coordinate.getZ()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getZ());
    }

    void rotateMultipoleFull(FComplexe vec[], const FReal theta, const FReal phi){
        FComplexe cell_rotate[SizeArray];
        for(int l = 0 ; l <= P ; ++l){
            for(int m = 0 ; m <= l ; ++m ){
                for(int k = -l ; k < 0 ; ++k){
                    // Equ 25
                    // Y{l,m} = SUM D{l,m,k} Y{l,k}
                    // z = p exp(i*O) = p(cos(O) + i sin(O))
                    const FReal factor = FMath::pow(-1,-k) * FMath::Sqrt(FReal(fact(l-k)*fact(l+k))/FReal(fact(l-m)*fact(l+m)));
                    const FReal d_lmk = analytical(FMath::Cos(theta),
                                                   FMath::Sin(theta),l,m,k);
                    const FReal minus_k = FReal(-k);
                    const FComplexe Dlmk(factor * d_lmk * FMath::Cos(phi * minus_k),
                                         factor * d_lmk * FMath::Sin(phi * minus_k));
                    cell_rotate[atLm(l,m)].addMul( Dlmk , vec[atLm(l,-k)].conjugate());
                }
                for(int k = 0 ; k <= l ; ++k){
                    // Equ 25
                    // Y{l,m} = SUM D{l,m,k} Y{l,k}
                    // z = p exp(i*O) = p(cos(O) + i sin(O))
                    const FReal factor = FMath::Sqrt(FReal(fact(l-k)*fact(l+k))/FReal(fact(l-m)*fact(l+m)));
                    const FReal d_lmk = analytical(FMath::Cos(theta),
                                                   FMath::Sin(theta),l,m,k);
                    const FReal minus_k = FReal(-k);
                    const FComplexe Dlmk(factor * d_lmk * FMath::Cos(phi * minus_k),
                                         factor * d_lmk * FMath::Sin(phi * minus_k));
                    cell_rotate[atLm(l,m)].addMul( Dlmk , vec[atLm(l,k)]);
                }
            }
        }
        FMemUtils::copyall(vec,cell_rotate,SizeArray);
    }

    void rotateMultipolePhi(FComplexe vec[], const FReal phi){
        FComplexe cell_rotate[SizeArray];
        for(int l = 0 ; l <= P ; ++l){
            for(int m = 0 ; m <= l ; ++m ){
                const FReal k = FReal(m);
                const FComplexe exp_ikphi(FMath::Cos(phi * k), FMath::Sin(phi * k));
                cell_rotate[atLm(l,m)].equalMul(exp_ikphi , vec[atLm(l,m)]);
            }
        }
        FMemUtils::copyall(vec,cell_rotate,SizeArray);
    }

    void rotateMultipolePhi2(FComplexe vec[], const FReal phi){
        FComplexe cell_rotate[SizeArray];
        for(int l = 0 ; l <= P ; ++l){
            for(int m = 0 ; m <= l ; ++m ){
                const FReal minus_m = FReal(-m);
                const FComplexe exp_minus_imphi(FMath::Cos(phi * minus_m), FMath::Sin(phi * minus_m));
                cell_rotate[atLm(l,m)].equalMul(exp_minus_imphi , vec[atLm(l,m)]);
            }
        }
        FMemUtils::copyall(vec,cell_rotate,SizeArray);
    }

    void rotateTaylorPhi(FComplexe vec[], const FReal phi){
        FComplexe cell_rotate[SizeArray];
        for(int l = 0 ; l <= P ; ++l){
            for(int m = 0 ; m <= l ; ++m ){
                const FReal k = FReal(m);
                const FComplexe exp_ikphi(FMath::Cos(phi * k), FMath::Sin(phi * k));
                cell_rotate[atLm(l,m)].equalMul(exp_ikphi , vec[atLm(l,m)]);
            }
        }
        FMemUtils::copyall(vec,cell_rotate,SizeArray);
    }

    void rotateMultipoleTheta(FComplexe vec[], const FReal theta){
        FComplexe cell_rotate[SizeArray];
        for(int l = 0 ; l <= P ; ++l){
            for(int m = 0 ; m <= l ; ++m ){
                FReal w_lkm_real = 0.0;
                FReal w_lkm_imag = 0.0;

                for(int k = -l ; k < 0 ; ++k){
                    const FReal d_lmk = analytical(FMath::Cos(theta),FMath::Sin(theta),l,m,k);
                    w_lkm_real += FMath::pow(-1,-k) * d_lmk * vec[atLm(l,-k)].getReal(); // k<0 => Conjugate * -1^k
                    w_lkm_imag -= FMath::pow(-1,-k) * d_lmk * vec[atLm(l,-k)].getImag(); // k<0 => Conjugate * -1^k
                    //printf("l%d m%d k%d, dlkm = %f, vec = %f\n", l,m,k,d_lmk,vec[atLm(l,-k)].getReal());
                }
                for(int k = 0 ; k <= l ; ++k){
                    const FReal d_lmk = analytical(FMath::Cos(theta),FMath::Sin(theta),l,m,k);
                    w_lkm_real += d_lmk * vec[atLm(l,k)].getReal();
                    w_lkm_imag += d_lmk * vec[atLm(l,k)].getImag();
                    //printf("l%d m%d k%d, dlkm = %f, vec = %f\n", l,m,k,d_lmk,vec[atLm(l,k)].getReal());
                }
                //printf("l%d m%d res = %f\n", l,m,w_lkm_real);
                cell_rotate[atLm(l,m)].setRealImag(w_lkm_real, w_lkm_imag);
            }
        }
        FMemUtils::copyall(vec,cell_rotate,SizeArray);
    }

    void rotateMultipoleTheta2(FComplexe vec[], const FReal theta){
        FComplexe cell_rotate[SizeArray];
        for(int l = 0 ; l <= P ; ++l){
            for(int m = 0 ; m <= l ; ++m ){
                FReal w_lkm_real = 0.0;
                FReal w_lkm_imag = 0.0;

                for(int k = -l ; k < 0 ; ++k){
                    const FReal d_lmk = analytical(FMath::Cos(theta),FMath::Sin(theta),l,m,k);
                    const FReal factor = FMath::Sqrt(FReal(fact(l-k)*fact(l+k))/FReal(fact(l-abs(m))*fact(l+abs(m))));
                    w_lkm_real += FMath::pow(-1,-k) * factor * d_lmk * vec[atLm(l,-k)].getReal(); // k<0 => Conjugate * -1^k
                    w_lkm_imag -= FMath::pow(-1,-k) * factor * d_lmk * vec[atLm(l,-k)].getImag(); // k<0 => Conjugate * -1^k
                }
                for(int k = 0 ; k <= l ; ++k){
                    const FReal d_lmk = analytical(FMath::Cos(theta),FMath::Sin(theta),l,m,k);
                    const FReal factor = FMath::Sqrt(FReal(fact(l-k)*fact(l+k))/FReal(fact(l-abs(m))*fact(l+abs(m))));
                    w_lkm_real += factor * d_lmk * vec[atLm(l,k)].getReal();
                    w_lkm_imag += factor * d_lmk * vec[atLm(l,k)].getImag();
                }
                //printf("l%d m%d res = %f\n", l,m,w_lkm_real);
                cell_rotate[atLm(l,m)].setRealImag(w_lkm_real, w_lkm_imag);
            }
        }
        FMemUtils::copyall(vec,cell_rotate,SizeArray);
    }

    void rotateMultipole(FComplexe vec[], FReal theta, FReal phi){
        rotateMultipolePhi2(vec,-phi);
        rotateMultipoleTheta2(vec,-theta);
    }
    void deRotateMultipole(FComplexe vec[], const FReal theta, const FReal phi){
        rotateMultipoleTheta2(vec,theta);
        rotateMultipolePhi2(vec,phi);
    }

    void rotateTaylor(FComplexe vec[], const FReal theta, const FReal phi){
        rotateTaylorPhi(vec,phi);
        rotateMultipoleTheta(vec,theta);
    }
    void deRotateTaylor(FComplexe vec[], const FReal theta, const FReal phi){
        rotateMultipoleTheta(vec,FMath::FPi*2 - theta);
        rotateTaylorPhi(vec,FMath::FPi*2 - phi);
    }

    /*void rotateTaylor(FComplexe vec[], const FReal theta, const FReal phi){
        FComplexe cell_rotate[SizeArray];
        for(int l = 0 ; l <= P ; ++l){
            for(int m = 0 ; m <= l ; ++m ){
                for(int k = -l ; k < 0 ; ++k){
                    // Equ 25
                    // Y{l,m} = SUM D{l,m,k} Y{l,k}
                    // z = p exp(i*O) = p(cos(O) + i sin(O))
                    const FReal d_lmk = analytical(FMath::Cos(theta),FMath::Sin(theta),l,m,k);
                    const FReal factor = FMath::Sqrt( FReal(fact(l-m)*fact(l+m)) / FReal(fact(l-k)*fact(l+k)) );
                    const FComplexe Dlmk( factor * d_lmk * FMath::Cos(phi * FReal(k)),
                                          factor * d_lmk * FMath::Sin(phi * FReal(k)));
                    cell_rotate[atLm(l,m)].addMul( Dlmk , vec[atLm(l,-k)].conjugate() );
                }
                for(int k = 0 ; k <= l ; ++k){
                    // Equ 25
                    // Y{l,m} = SUM D{l,m,k} Y{l,k}
                    // z = p exp(i*O) = p(cos(O) + i sin(O))
                    const FReal d_lmk = analytical(FMath::Cos(theta),FMath::Sin(theta),l,m,k);
                    const FReal factor = FMath::Sqrt( FReal(fact(l-m)*fact(l+m)) / FReal(fact(l-k)*fact(l+k)) );
                    const FComplexe Dlmk( factor * d_lmk * FMath::Cos(phi * FReal(k)),
                                          factor * d_lmk * FMath::Sin(phi * FReal(k)));
                    cell_rotate[atLm(l,m)].addMul( Dlmk , vec[atLm(l,k)] );
                }
            }
        }
        FMemUtils::copyall(vec,cell_rotate,SizeArray);
    }*/

    const FReal boxWidth;       //< the box width at leaf level
    const int   treeHeight;     //< The height of the tree
    const FReal widthAtLeafLevel;
    const FReal widthAtLeafLevelDiv2;
    const FPoint boxCorner;

    FSmartPointer<FSpherical> childrenPosition;
    FSmartPointer<FSpherical> interactionsPosition;

    void preComputePosition(){
        childrenPosition = new FSpherical[8 * (treeHeight-1)];
        {
            FReal subBoxWidth = widthAtLeafLevelDiv2;
            // For each
            for(int idxLevel = treeHeight-2 ; idxLevel > 0 ; --idxLevel){
                for(int idxChild = 0 ; idxChild < 8 ; ++idxChild ){
                    // coord from child to parent
                    const FReal x = FReal((idxChild&4)? -1.0 : 1.0) * subBoxWidth;
                    const FReal y = FReal((idxChild&2)? -1.0 : 1.0) * subBoxWidth;
                    const FReal z = FReal((idxChild&1)? -1.0 : 1.0) * subBoxWidth;
                    const FPoint relativePosition( x , y , z );

                    childrenPosition[(idxLevel-1)*8 + idxChild] = FSpherical(relativePosition);

                    /*printf("Level %d, child %d (%f,%f,%f), theta %f phi %f \n",
                           idxLevel, idxChild, x, y, z, childrenPosition[(idxLevel-1)*8 + idxChild].getTheta(),
                           childrenPosition[(idxLevel-1)*8 + idxChild].getPhi());*/
                }
                subBoxWidth *= FReal(2.0);
            }
        }

        interactionsPosition = new FSpherical[343 * (treeHeight-1)];
        {
            FReal boxWidthAtLevel = widthAtLeafLevel;
            for(int idxLevel = treeHeight-1 ; idxLevel > 0 ; --idxLevel){
                for(int idxX = -3 ; idxX <= 3 ; ++idxX ){
                    for(int idxY = -3 ; idxY <= 3 ; ++idxY ){
                        for(int idxZ = -3 ; idxZ <= 3 ; ++idxZ ){
                            if( idxX || idxY || idxZ ){
                                const FPoint relativePosition( -FReal(idxX)*boxWidthAtLevel, -FReal(idxY)*boxWidthAtLevel, -FReal(idxZ)*boxWidthAtLevel);
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

    FSpherical getSphericalChild(const int idxLevel, const int position) const {
        return childrenPosition[(idxLevel-1)*8 + position];
    }

    FSpherical getSphericalInteraction(const int idxLevel, const int position) const {
        return interactionsPosition[(idxLevel-1)*343 + position];
    }

public:


    FRotationKernel(const int inDevP, const int inTreeHeight, const FReal inBoxWidth, const FPoint& inBoxCenter) :
        boxWidth(inBoxWidth),
        treeHeight(inTreeHeight),
        widthAtLeafLevel(inBoxWidth/FReal(1 << (inTreeHeight-1))),
        widthAtLeafLevelDiv2(widthAtLeafLevel/2),
        boxCorner(inBoxCenter.getX()-(inBoxWidth/2),inBoxCenter.getY()-(inBoxWidth/2),inBoxCenter.getZ()-(inBoxWidth/2))
        {

        preComputePosition();
    }

    /** Default destructor */
    virtual ~FRotationKernel(){
    }

    /** P2M
      * From equation 10, fmmparadix.pdf
      */
    void P2M(CellClass* const inPole, const ContainerClass* const inParticles ) {
        //kernel.P2M(inPole, inParticles);

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

            // w{l,m} = q a^l/(l+|m|)! P{l,m} exp(-i m Betha)
            for(int l = 0 ; l <= P ; ++l ){
                for(int m = 0 ; m <= l ; ++m ){
                    const FReal magnitude = q * (FMath::pow( a , l )/fact(l+m))
                            * legendre[atLm(l,m)];
                    w[atLm(l,m)].setReal(magnitude * FMath::Cos(FReal(-m) * sph.getPhi()));
                    w[atLm(l,m)].setImag(magnitude * FMath::Sin(FReal(-m) * sph.getPhi()));
                }
            }
            // Goto next particle
            iterParticle.gotoNext();
        }
    }

    /** M2M
      * Equation 29, from fmm paradix
      */
    void M2M(CellClass* const FRestrict inPole, const CellClass*const FRestrict *const FRestrict inChildren, const int inLevel) {
        // A buffer to copy the source w
        FComplexe source_w[SizeArray];
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            // if child exists
            if(inChildren[idxChild]){
                // Copy the source
                FMemUtils::copyall(source_w, inChildren[idxChild]->getMultipole(), SizeArray);
                //delete me
                printf("----------------- child %d \n", idxChild);
                for(int l = 0 ; l <= P ; ++l ){
                    for(int m = 0 ; m <= l ; ++m ){
                        printf("0 - [%d][%d] %f i%f\t", l, m, source_w[atLm(l,m)].getReal(), source_w[atLm(l,m)].getImag());
                    }
                    printf("\n");
                }
                //end delete me
                // rotate it
                const FSpherical sph = getSphericalChild(inLevel, idxChild);
                rotateMultipole(source_w, sph.getTheta(), sph.getPhi());
                const FReal b = -sph.getR();

                //delete me
                printf("b %f , theta %f, Phi %f \n", b,  sph.getTheta(), sph.getPhi());
                for(int l = 0 ; l <= P ; ++l ){
                    for(int m = 0 ; m <= l ; ++m ){
                        printf("1 - [%d][%d] %f i%f\t", l, m, source_w[atLm(l,m)].getReal(), source_w[atLm(l,m)].getImag());
                    }
                    printf("\n");
                }
                //end delete me

                // Translate it
                FComplexe target_w[SizeArray];
                for(int l = 0 ; l <= P ; ++l ){
                    for(int m = 0 ; m <= l ; ++m ){
                        // w{l,m}(a+b) = sum(j=m:l, b^(l-j)/(l-j)! w{j,m}(a)
                        FReal w_lm_real = 0.0;
                        FReal w_lm_imag = 0.0;
                        printf("T %d %d >> ", l , m);
                        for(int j = m ; j <= l ; ++j ){
                            const FReal coef = (FMath::pow(b,l-j) / FReal(fact(l-j)));
                            w_lm_real += coef * source_w[atLm(j,m)].getReal();
                            w_lm_imag += coef * source_w[atLm(j,m)].getImag();
                            printf(" \t use : %d %d (l-j %d) coef %f at %d, ",j,m,l-j, coef, atLm(j,m));
                        }
                        printf("\n");
                        target_w[atLm(l,m)].setRealImag(w_lm_real,w_lm_imag);
                    }
                }
                //FMemUtils::copyall( target_w, source_w, SizeArray);
                //delete me
                for(int l = 0 ; l <= P ; ++l ){
                    for(int m = 0 ; m <= l ; ++m ){
                        printf("2 - [%d][%d] %f i%f\t", l, m, target_w[atLm(l,m)].getReal(), target_w[atLm(l,m)].getImag());
                    }
                    printf("\n");
                }
                //end delete me

                // Rotate it back
                deRotateMultipole(target_w, sph.getTheta(),sph.getPhi());
                //delete me
                for(int l = 0 ; l <= P ; ++l ){
                    for(int m = 0 ; m <= l ; ++m ){
                        printf("3 - [%d][%d] %f i%f\t", l, m, target_w[atLm(l,m)].getReal(), target_w[atLm(l,m)].getImag());
                    }
                    printf("\n");
                }
                //end delete me
                // Sum the result
                FMemUtils::addall( inPole->getMultipole(), target_w, SizeArray);

            }
        }
    }

    /** M2L
      * Equation 33, from fmmparadix.pdf
      */
    void M2L(CellClass* const FRestrict inLocal, const CellClass* inInteractions[], const int /*inSize*/, const int inLevel) {
        // To copy the multipole data
        FComplexe source_w[SizeArray];
        for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
            // if interaction exits
            if(inInteractions[idxNeigh]){
                // Copy multipole data into buffer
                FMemUtils::copyall(source_w, inInteractions[idxNeigh]->getMultipole(), SizeArray);

                // Rotate
                const FSpherical sph = getSphericalInteraction(inLevel, idxNeigh);
                rotateMultipole(source_w, sph.getTheta(), sph.getPhi());

                const FReal b = sph.getR();

                // Transfer to u
                FComplexe target_u[SizeArray];
                for(int l = 0 ; l <= P ; ++l ){
                    for(int m = 0 ; m <= l ; ++m ){
                        // u{l,m}(a-b) = sum(j=0:P, (j+l)!/b^(j+l+1) w{j,-m}(a)
                        FReal u_lm_real = 0.0;
                        FReal u_lm_imag = 0.0;
                        for(int j = 0 ; j <= P ; ++j ){
                            const FReal coef = (fact(j+l)/FMath::pow(b,j+l+1));
                            // conjugate because {l,-m} => {l,m} with -i
                            u_lm_real += coef * source_w[atLm(j,m)].getReal();
                            u_lm_imag -= coef * source_w[atLm(j,m)].getImag();
                        }
                        target_u[atLm(l,m)].setRealImag(u_lm_real,u_lm_imag);
                    }
                }

                // Rotate it back
                deRotateTaylor(target_u, FMath::FPi*2 - sph.getTheta(), FMath::FPi*2 - sph.getPhi());
                // Sum
                FMemUtils::addall(inLocal->getLocal(), target_u, SizeArray);
            }
        }
    }

    /** L2L
      * Equation 37, from fmmparadix.pdf
      */
    void L2L(const CellClass* const FRestrict inLocal, CellClass* FRestrict *const FRestrict  inChildren, const int inLevel) {
        // To copy the source local
        FComplexe source_u[SizeArray];
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            // if child exists
            if(inChildren[idxChild]){
                // Copy the local data into the buffer
                FMemUtils::copyall(source_u, inLocal->getLocal(), SizeArray);
                // Rotate
                const FSpherical sph = getSphericalChild(inLevel, idxChild);
                rotateTaylor(source_u, FMath::FPi*2 -sph.getTheta(), FMath::FPi*2 - sph.getPhi());
                const FReal b = sph.getR();

                // Translate
                FComplexe target_u[SizeArray];
                for(int l = 0 ; l <= P ; ++l ){
                    for(int m = 0 ; m <= l ; ++m ){
                        // u{l,m}(r-b) = sum(j=0:P, b^(j-l)/(j-l)! u{j,m}(r);
                        FReal u_lm_real = 0.0;
                        FReal u_lm_imag = 0.0;
                        for(int j = l ; j <= P ; ++j ){
                            const FReal coef = (FMath::pow(b,j-l) / fact(j-l));
                            u_lm_real += coef * source_u[atLm(j,m)].getReal();
                            u_lm_imag += coef * source_u[atLm(j,m)].getImag();
                        }
                        target_u[atLm(l,m)].setRealImag(u_lm_real,u_lm_imag);
                    }
                }
                // Rotate
                rotateTaylor(target_u, FMath::FPi*2 - sph.getTheta(), FMath::FPi*2 - sph.getPhi());
                // Sum in child
                FMemUtils::addall(inChildren[idxChild]->getLocal(), target_u, SizeArray);
            }
        }
    }

    /** L2P
      * Equation 13/14, from fmmparadix.pdf
      */
    void L2P(const CellClass* const inLocal, ContainerClass* const inParticles){
        // Take the local value from the cell
        const FComplexe* FRestrict const u = inLocal->getLocal();

        // For all particles in the leaf box
        typename ContainerClass::BasicIterator iterParticle(*inParticles);
        while( iterParticle.hasNotFinished()){
            // L2P
            ParticleClass& particle = iterParticle.data();
            FReal magnitude = 0;

            // E = sum( l = 0:P, sum(m = -l:l, u{l,m} ))
            for(int l = 0 ; l <= P ; ++l ){
                magnitude += u[atLm(l,0)].getReal();
                for(int m = 1 ; m <= l ; ++m ){
                    // we sum u{l,m} + u{l-m} = 2 * u{l,m} since
                    // for u{l,-m} = u*{l,m}
                    magnitude += 2 * u[atLm(l,m)].getReal();
                }
            }
            // inc potential
            particle.incPotential(magnitude);
            // progress
            iterParticle.gotoNext();
        }
    }


    /** Do nothing */
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

    /** P2P mutual interaction
      * F = q * q' / r²
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

    /** Do nothing */
    void P2PRemote(const FTreeCoordinate& inPosition,
                   ContainerClass* const FRestrict inTargets, const ContainerClass* const FRestrict inSources,
                   ContainerClass* const inNeighbors[27], const int inSize){
        //kernel.P2PRemote(inPosition, inTargets, inSources, inNeighbors,inSize);
    }

};


#endif // FROTATIONKERNEL_HPP
