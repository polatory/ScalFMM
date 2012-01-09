#ifndef FHARMONIC_HPP
#define FHARMONIC_HPP

#include "../Utils/FComplexe.hpp"

/** Todo move this class outside when kernels will be developed */
class FSpherical {
    FReal r, cosTheta, sinTheta, phi;
public:
    FSpherical(): r(0), cosTheta(0), sinTheta(0), phi(0){
    }

    explicit FSpherical(const F3DPosition& inVector){
        const FReal x2y2 = (inVector.getX() * inVector.getX()) + (inVector.getY() * inVector.getY());
        this->r = FMath::Sqrt( x2y2 + (inVector.getZ() * inVector.getZ()));
        this->phi = FMath::Atan2(inVector.getY(),inVector.getX());
        this->cosTheta = inVector.getZ() / outSphere->r; // cos_th = z/r
        this->sinTheta = FMath::Sqrt(x2y2) / outSphere->r; // sin_th = sqrt(x^2 + y^2)/r
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
    const int devP;
    const int expSize;

    FComplexe* harmonic;
    FComplexe* cosSin;
    FReal*     legendre;
    FComplexe* thetaDerivatedResult;
    FReal*     sphereHarmoInnerCoef;
    FReal*     sphereHarmoOuterCoef;


    void sphericalHarmonicInitialize(){
        FReal factOuter = 1.0;
        for(int idxP = 0 ; idxP <= devP; factOuter *= FReal(++idxP) ){
            sphereHarmoOuterCoef[idxP] = factOuter;
        }

        FReal* currentInner = sphereHarmoInnerCoef;
        FReal factInner = 1.0;
        FReal powN1idxP = 1.0;
        for(int idxP = 0 ; idxP <= this->devP ; factInner *= FReal(++idxP), powN1idxP = -powN1idxP){
            for(int idxMP = 0, fact_l_m = int(factInner); idxMP <= idxP ; fact_l_m *= idxP+(++idxMP), ++currentInner){
                *currentInner = powN1idxP / FReal(fact_l_m);
            }
        }
    }

    void computeCosSin(const FReal inPhi){
        for(int idxl = 0, idxlMod4 = 0; idxl <= devP ; ++idxl, ++idxlMod4){
            if(idxlMod4 == 4) idxlMod4 = 0;
            const FReal angle = FReal(idxl) * inPhi + this->PiArrayOuter[idxlMod4];

            cosSin[idxl].setReal( FMath::Sin(angle + FMath::FPiDiv2) );
            cosSin[idxl].setImag( FMath::Sin(angle) );
        }
    }

    void computeLegendre(const FReal inCosTheta, const FReal inSinTheta){
        int idxCurrent = 0;
        legendre[idxCurrent++] = 1.0; // P_0^0(cosTheta) = 1

        legendre[idxCurrent++] = inCosTheta; // P_1^{0} using (3)

        // Compute P_1^1 using (2 bis) and store it into results_array
        const FReal invSinTheta = -inSinTheta; // somx2 => -sinTheta
        legendre[idxCurrent++] = invSinTheta;

        // l>1:
        int idxCurrent1m = 1; //pointer on P_{l-1}^m P_1^0
        int idxCurrent2m = 0; //pointer on P_{l-2}^m P_0^0
        FReal fact = 3.0;

        // Remark: p_results_array_l_minus_1_m and p_results_array_l_minus_2_m
        // just need to be incremented at each iteration.
        for(int idxl = 2; idxl <= devP ; ++idxl ){
            for( int idxm = 0; idxm <= idxl - 2 ; ++idxm , ++idxCurrent , ++idxCurrent1m , ++idxCurrent2m ){
                // Compute P_l^m, l >= m+2, using (1) and store it into results_array:
                legendre[idxCurrent] = (inCosTheta * FReal( 2 * idxl - 1 ) * legendre[idxCurrent1m] - FReal( idxl + idxm - 1 )
                                          * legendre[idxCurrent2m] ) / FReal( idxl - idxm );
            }
            // p_results_array_l_minus_1_m now points on P_{l-1}^{l-1}

            // Compute P_l^{l-1} using (3) and store it into ptrResults:
            legendre[idxCurrent++] = inCosTheta * FReal( 2 * idxl - 1 ) * legendre[idxCurrent1m];

            // Compute P_l^l using (2 bis) and store it into results_array:
            legendre[idxCurrent++] = fact * invSinTheta * legendre[idxCurrent1m];

            fact += FReal(2.0);
            ++idxCurrent1m;
        }
    }

public:
    explicit FHarmonic(const int inDevP)
        : devP(inDevP),expSize(int(((inDevP)+1) * ((inDevP)+2) * 0.5)),
          harmonic(0), cosSin(0), legendre(0), thetaDerivatedResult(0),
          sphereHarmoInnerCoef(0), sphereHarmoOuterCoef(0) {

        harmonic = new FComplexe[expSize];
        cosSin   = new FComplexe[2 * devP + 1];
        legendre = new FReal[expSize];

        thetaDerivatedResult = new FComplexe[expSize];
        sphereHarmoInnerCoef = new FReal[expSize];
        sphereHarmoOuterCoef = new FReal[devP + 1];

        sphericalHarmonicInitialize();
    }

    ~FHarmonic(){
        delete[] harmonic;
        delete[] cosSin;
        delete[] legendre;
        delete[] thetaDerivatedResult;
        delete[] sphereHarmoInnerCoef;
        delete[] sphereHarmoOuterCoef;
    }

    FComplexe* data(){
        return harmonic;
    }

    const FComplexe* data() const {
        return harmonic;
    }

    FComplexe& result(const int index){
        return harmonic[index];
    }

    const FComplexe& result(const int index) const{
        return harmonic[index];
    }

    void computeInner(const FSpherical& inSphere){
        computeCosSin(inSphere.getPhi());

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
        computeCosSin(inSphere.getPhi());

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

    void harmonicInnerTheta(const FSpherical& inSphere){
        computeCosSin(inSphere.getPhi());

        // p_associated_Legendre_function_Array
        computeLegendre(inSphere.getCosTheta(), inSphere.getSinTheta());

        FComplexe *p_term = harmonic;
        FComplexe *p_theta_derivated_term = thetaDerivatedResult;
        FReal *p_spherical_harmonic_Inner_coefficients_array = sphereHarmoInnerCoef;
        FReal *ptr_associated_Legendre_function_Array = legendre;
        FReal *start_ptr_associated_Legendre_function_Array = ptr_associated_Legendre_function_Array;

        // r^l
        FReal r_l = 1.0;
        for (int l = 0 ; l <= FMB_Info_P ; ++l, r_l *= inSphere.r){
            FReal magnitude;
            // m<l:
            int m = 0;
            for(; m < l ; ++m, ++p_term, ++p_theta_derivated_term, ++p_spherical_harmonic_Inner_coefficients_array, ++ptr_associated_Legendre_function_Array){
                magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * r_l * (*ptr_associated_Legendre_function_Array);

                // Computation of Inner_l^m(r, theta, phi):
                p_term->setReal( magnitude * cosSin[m].getReal());
                p_term->setImag( magnitude * cosSin[m].getImag());

                // Computation of {\partial Inner_l^m(r, theta, phi)}/{\partial theta}:
                magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * r_l * ((FReal(l)*inSphere.cosTheta*(*ptr_associated_Legendre_function_Array)
                                                                                       - FReal(l+m)*(*(start_ptr_associated_Legendre_function_Array + expansion_Redirection_array_for_j[l-1] + m) )) / inSphere.sinTheta);
                p_theta_derivated_term->setReal(magnitude * cosSin[m].getReal());
                p_theta_derivated_term->setImag(magnitude * cosSin[m].getImag());

            }

            // m=l:
            // Computation of Inner_m^m(r, theta, phi):
            magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * r_l * (*ptr_associated_Legendre_function_Array);
            p_term->setReal(magnitude * cosSin[m].getReal());
            p_term->setImag(magnitude * cosSin[m].getImag());

            // Computation of {\partial Inner_m^m(r, theta, phi)}/{\partial theta}:
            magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * FReal(r_l) * (FReal(m) * inSphere.cosTheta * (*ptr_associated_Legendre_function_Array) / inSphere.sinTheta);
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
