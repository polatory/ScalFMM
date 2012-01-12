#ifndef FHARMONIC_HPP
#define FHARMONIC_HPP

#include "../Utils/FComplexe.hpp"

/** Todo move this class outside when kernels will be developed */
class FSpherical {
    FReal r, cosTheta, sinTheta, phi;

public:
    explicit FSpherical(const F3DPosition& inVector){
        const FReal x2y2 = (inVector.getX() * inVector.getX()) + (inVector.getY() * inVector.getY());
        this->r = FMath::Sqrt( x2y2 + (inVector.getZ() * inVector.getZ()));
        this->phi = FMath::Atan2(inVector.getY(),inVector.getX());
        this->cosTheta = inVector.getZ() / r;
        this->sinTheta = FMath::Sqrt(x2y2) / r;
    }

    FReal getR() const{
        return r;
    }

    FReal getCosTheta() const{
        return cosTheta;
    }

    FReal getSinTheta() const{
        return sinTheta;
    }

    FReal getPhi() const{
        return phi;
    }
};



class FHarmonic {
    const int devP;     //< P
    const int expSize;  //< Exponen Size

    FComplexe* harmonic;//< Harmonic Result
    FComplexe* cosSin;  //< Cos/Sin precomputed values
    FReal*     legendre;//< Legendre results

    FComplexe* thetaDerivatedResult; //< the theta derivated result
    FReal*     sphereHarmoInnerCoef; //< sphere innter pre computed coefficients
    FReal*     sphereHarmoOuterCoef; //< sphere outer pre computed coefficients

    int* preExpRedirJ; //< A indirection to acceess the good value in the preX complexes vector

    /** Allocate and init */
    void allocAndInit(){
        harmonic = new FComplexe[expSize];
        cosSin   = new FComplexe[2 * devP + 1];
        legendre = new FReal[expSize];

        thetaDerivatedResult = new FComplexe[expSize];
        sphereHarmoInnerCoef = new FReal[int(((devP*2)+1) * ((devP*2)+2) * 0.5)];
        sphereHarmoOuterCoef = new FReal[devP + 1];

        // Pre compute coef
        FReal factOuter = 1.0;
        for(int idxP = 0 ; idxP <= devP; factOuter *= FReal(++idxP) ){
            sphereHarmoOuterCoef[idxP] = factOuter;
        }

        // Pre compute coef
        FReal* currentInner = sphereHarmoInnerCoef;
        FReal factInner = 1.0;
        FReal powN1idxP = 1.0;
        for(int idxP = 0 ; idxP <= this->devP ; factInner *= FReal(++idxP), powN1idxP = -powN1idxP){
            for(int idxMP = 0, fact_l_m = int(factInner); idxMP <= idxP ; fact_l_m *= idxP+(++idxMP), ++currentInner){
                *currentInner = powN1idxP / FReal(fact_l_m);
            }
        }

        // Smart redirection
        preExpRedirJ = new int[2 * devP + 1];
        for( int h = 0; h <= (2 * devP) ; ++h ){
            preExpRedirJ[h] = static_cast<int>( h * ( h + 1 ) * 0.5 );
        }
    }

    /** Compute cos/sin from phi */
    void computeCosSin(const FReal inPhi, const FReal piArray[4]){
        for(int idxl = 0 ; idxl <= devP ; ++idxl){
            const FReal angle = FReal(idxl) * inPhi + piArray[idxl & 0x3];
            cosSin[idxl].setReal( FMath::Sin(angle + FMath::FPiDiv2) );
            cosSin[idxl].setImag( FMath::Sin(angle) );
        }
    }

    /** Compute legendre */
    void computeLegendre(const FReal inCosTheta, const FReal inSinTheta){
        const FReal invSinTheta = -inSinTheta;

        legendre[0] = 1.0;        // P_0,0(cosTheta) = 1
        legendre[1] = inCosTheta; // P_1,0(cosTheta) = cosTheta
        legendre[2] = invSinTheta;// P_1,1(cosTheta) = -sinTheta

        int idxCurrentLM  = 3; //current pointer on P_l,m
        int idxCurrentL1M = 1; //pointer on P_{l-1},m => P_1,0
        int idxCurrentL2M = 0; //pointer on P_{l-2},m => P_0,0
        FReal fact = 3.0;

        for(int idxl = 2; idxl <= devP ; ++idxl ){
            // m from 0 to l - 2
            for( int idxm = 0; idxm <= idxl - 2 ; ++idxm ){
                //         cosTheta x (2 x l - 1) x P_l-1_m - (l + m - 1) x P_l-2_m
                // P_l_m = --------------------------------------------------------
                //                               l - m
                legendre[idxCurrentLM] = (inCosTheta * FReal( 2 * idxl - 1 ) * legendre[idxCurrentL1M]
                                          - FReal( idxl + idxm - 1 ) * legendre[idxCurrentL2M] )
                                          / FReal( idxl - idxm );
                // progress
                ++idxCurrentLM;
                ++idxCurrentL1M;
                ++idxCurrentL2M;
            }

            // Compute P_l,{l-1}
            legendre[idxCurrentLM++] = inCosTheta * FReal( 2 * idxl - 1 ) * legendre[idxCurrentL1M];
            // Compute P_l,l
            legendre[idxCurrentLM++] = fact * invSinTheta * legendre[idxCurrentL1M];

            fact += FReal(2.0);
            ++idxCurrentL1M;
        }
    }

    /** Forbid copy operator */
    FHarmonic& operator=(const FHarmonic&){ return *this;}

public:
    /////////////////////////////////////////////////////////////////
    // Constructor & Accessor
    /////////////////////////////////////////////////////////////////

    explicit FHarmonic(const int inDevP)
        : devP(inDevP),expSize(int(((inDevP)+1) * ((inDevP)+2) * 0.5)),
          harmonic(0), cosSin(0), legendre(0), thetaDerivatedResult(0),
          sphereHarmoInnerCoef(0), sphereHarmoOuterCoef(0), preExpRedirJ(0)  {

        allocAndInit();
    }

    FHarmonic(const FHarmonic& other)
        : devP(other.devP),expSize(int(((other.devP)+1) * ((other.devP)+2) * 0.5)),
          harmonic(0), cosSin(0), legendre(0), thetaDerivatedResult(0),
          sphereHarmoInnerCoef(0), sphereHarmoOuterCoef(0), preExpRedirJ(0)  {

        allocAndInit();
    }

    ~FHarmonic(){
        delete[] harmonic;
        delete[] cosSin;
        delete[] legendre;
        delete[] thetaDerivatedResult;
        delete[] sphereHarmoInnerCoef;
        delete[] sphereHarmoOuterCoef;
        delete[] preExpRedirJ;
    }

    int getExpSize() const{
        return expSize;
    }

    FComplexe* result(){
        return harmonic;
    }

    const FComplexe* result() const {
        return harmonic;
    }

    FComplexe& result(const int index){
        return harmonic[index];
    }

    const FComplexe& result(const int index) const{
        return harmonic[index];
    }

    FComplexe* resultThetaDerivated(){
        return thetaDerivatedResult;
    }

    const FComplexe* resultThetaDerivated() const {
        return thetaDerivatedResult;
    }

    FComplexe& resultThetaDerivated(const int index){
        return thetaDerivatedResult[index];
    }

    const FComplexe& resultThetaDerivated(const int index) const{
        return thetaDerivatedResult[index];
    }

    int getPreExpRedirJ(const int index) const{
        return preExpRedirJ[index];
    }

    /////////////////////////////////////////////////////////////////
    // Computation
    /////////////////////////////////////////////////////////////////

    void computeInner(const FSpherical& inSphere){
        const FReal PiArrayInner[4] = {0, FMath::FPiDiv2, FMath::FPi, -FMath::FPiDiv2};
        computeCosSin(inSphere.getPhi(), PiArrayInner);

        // p_associated_Legendre_function_Array
        computeLegendre(inSphere.getCosTheta(), inSphere.getSinTheta());

        FComplexe* currentResult = harmonic;
        int idxLegendre = 0;//ptr_associated_Legendre_function_Array
        int idxSphereHarmoCoef = 0;
        FReal idxRl = 1.0 ;

        for(int idxl = 0; idxl <= devP ; ++idxl, idxRl *= inSphere.getR()){
            for(int idxm = 0 ; idxm <= idxl ; ++idxm, ++currentResult, ++idxSphereHarmoCoef, ++idxLegendre){
                const FReal magnitude = sphereHarmoInnerCoef[idxSphereHarmoCoef] * idxRl * legendre[idxLegendre];
                currentResult->setReal( magnitude * cosSin[idxm].getReal() );
                currentResult->setImag( magnitude * cosSin[idxm].getImag() );
            }
        }
    }

    void computeOuter(const FSpherical& inSphere){
        const FReal PiArrayOuter[4] = {0, -FMath::FPiDiv2, FMath::FPi, FMath::FPiDiv2};
        computeCosSin(inSphere.getPhi(), PiArrayOuter);

        // p_associated_Legendre_function_Array
        computeLegendre(inSphere.getCosTheta(), inSphere.getSinTheta());

        int idxLegendre = 0;
        FComplexe* currentResult = harmonic;
        const FReal invR = 1/inSphere.getR();
        FReal idxRl1 = invR;

        for(int idxl = 0 ; idxl <= devP ; ++idxl, idxRl1 *= invR){
            for(int idxm = 0 ; idxm <= idxl ; ++idxm, ++currentResult, ++idxLegendre){
                const FReal magnitude = this->sphereHarmoOuterCoef[idxl-idxm] * idxRl1 * legendre[idxLegendre];
                currentResult->setReal( magnitude * cosSin[idxm].getReal() );
                currentResult->setImag( magnitude * cosSin[idxm].getImag() );
            }
        }
    }

    /** spherical_harmonic_Inner_and_theta_derivated
        * Returns the value of the partial derivative of the spherical harmonic
        *relative to theta. We have for all m such that -(l-1) <= m <= l-1 :
        *
        *(d H_l^m(theta, phi))/(d theta)
        *= (-1)^m sqrt((l-|m|)!/(l+|m|)!) (d P_l^{|m|}(cos theta))/(d theta) e^{i.m.phi}
        *= (-1)^m sqrt((l-|m|)!/(l+|m|)!) 1/sqrt{1-cos^2 theta} [l cos theta P_l^{|m|}(cos theta) - (l+|m|) P_{l-1}^{|m|}(cos theta) ] e^{i.m.phi}
        *Since theta is in the range [0, Pi], we have: sin theta > 0 and therefore
        *sqrt{1-cos^2 theta} = sin theta. Thus:
        *
        *(d H_l^m(theta, phi))/(d theta)
        *= (-1)^m sqrt((l-|m|)!/(l+|m|)!) 1/(sin theta) [l cos theta P_l^{|m|}(cos theta) - (l+|m|) P_{l-1}^{|m|}(cos theta) ] e^{i.m.phi}
        *For |m|=l, we have~:
        *(d H_l^l(theta, phi))/(d theta)
        *= (-1)^m sqrt(1/(2l)!) 1/(sin theta) [l cos theta P_l^l(cos theta) ] e^{i.m.phi}
        *
        *Remark: for 0<m<=l:
        *(d H_l^{-m}(theta, phi))/(d theta) = [(d H_l^{-m}(theta, phi))/(d theta)]*
        *
        *
        *
        *Therefore, we have for (d Inner_l^m(r, theta, phi))/(d theta):
        *
        *|m|<l: (d Inner_l^m(r, theta, phi))/(d theta) =
        *(i^m (-1)^l / (l+|m|)!) 1/(sin theta) [l cos theta P_l^{|m|}(cos theta) - (l+|m|) P_{l-1}^{|m|}(cos theta) ] e^{i.m.phi} r^l
        *|m|=l: (d Inner_l^m(r, theta, phi))/(d theta) =
        *(i^m (-1)^l / (l+|m|)!) 1/(sin theta) [l cos theta P_l^l(cos theta) ] e^{i.m.phi} r^l
        *
        *
      */
    void computeInnerTheta(const FSpherical& inSphere){
        const FReal PiArrayInner[4] = {0, FMath::FPiDiv2, FMath::FPi, -FMath::FPiDiv2};
        computeCosSin(inSphere.getPhi(),PiArrayInner);

        // p_associated_Legendre_function_Array
        computeLegendre(inSphere.getCosTheta(), inSphere.getSinTheta());

        FComplexe *p_term = harmonic;
        FComplexe *p_theta_derivated_term = thetaDerivatedResult;
        FReal *p_spherical_harmonic_Inner_coefficients_array = sphereHarmoInnerCoef;
        FReal *ptr_associated_Legendre_function_Array = legendre;
        FReal *start_ptr_associated_Legendre_function_Array = ptr_associated_Legendre_function_Array;

        // r^l
        FReal r_l = 1.0;
        for (int l = 0 ; l <= FMB_Info_P ; ++l, r_l *= inSphere.getR()){
            FReal magnitude;
            // m<l:
            int m = 0;
            for(; m < l ; ++m, ++p_term, ++p_theta_derivated_term, ++p_spherical_harmonic_Inner_coefficients_array, ++ptr_associated_Legendre_function_Array){
                magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * r_l * (*ptr_associated_Legendre_function_Array);

                // Computation of Inner_l^m(r, theta, phi):
                p_term->setReal( magnitude * cosSin[m].getReal());
                p_term->setImag( magnitude * cosSin[m].getImag());

                // Computation of {\partial Inner_l^m(r, theta, phi)}/{\partial theta}:
                magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * r_l * ((FReal(l)*inSphere.getCosTheta()*(*ptr_associated_Legendre_function_Array)
                                                                                       - FReal(l+m)*(*(start_ptr_associated_Legendre_function_Array + getPreExpRedirJ(l-1) + m) )) / inSphere.getSinTheta());
                p_theta_derivated_term->setReal(magnitude * cosSin[m].getReal());
                p_theta_derivated_term->setImag(magnitude * cosSin[m].getImag());

            }

            // m=l:
            // Computation of Inner_m^m(r, theta, phi):
            magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * r_l * (*ptr_associated_Legendre_function_Array);
            p_term->setReal(magnitude * cosSin[m].getReal());
            p_term->setImag(magnitude * cosSin[m].getImag());

            // Computation of {\partial Inner_m^m(r, theta, phi)}/{\partial theta}:
            magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * FReal(r_l) * (FReal(m) * inSphere.getCosTheta() * (*ptr_associated_Legendre_function_Array) / inSphere.getSinTheta());
            p_theta_derivated_term->setReal(magnitude * cosSin[m].getReal());
            p_theta_derivated_term->setImag(magnitude * cosSin[m].getImag());

            ++p_term;
            ++p_theta_derivated_term;
            ++p_spherical_harmonic_Inner_coefficients_array;
            ++ptr_associated_Legendre_function_Array;
        }
    }

};



#endif // FHARMONIC_HPP
