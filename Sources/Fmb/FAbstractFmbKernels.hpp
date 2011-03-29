#ifndef FABSTRACTFMBKERNELS_HPP
#define FABSTRACTFMBKERNELS_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Core/FAbstractKernels.hpp"

#include "../Containers/FTreeCoordinate.hpp"
#include "../Containers/FList.hpp"

#include "../Utils/F3DPosition.hpp"
#include "../Utils/FComplexe.hpp"
#include "../Utils/FMath.hpp"

#include <iostream>


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractFmbKernels
* @brief
* Please read the license
*
*/
template< class ParticuleClass, class CellClass>
class FAbstractFmbKernels : public FAbstractKernels<ParticuleClass,CellClass> {
protected:
    // P is a input parameter
    static const int FMB_Info_P = 2;

    // _GRAVITATIONAL_
    static const int FMB_Info_eps_soft_square = 1;

    // LMax is used in several algorithm
    // it is just a copy of FMB_Info.P
    static const int LMax = FMB_Info_P;

    // Can be 2 * FMB_Info_P if user ask to
    static const int FMB_Info_M2L_P = FMB_Info_P;
    static const int FMB_Info_M2L_exp_size = ((FMB_Info_M2L_P)+1) * ((FMB_Info_M2L_P)+2) * 0.5;

    // Default value set in main
    static const int FMB_Info_ws = 1;

    // INTERACTION_LIST_SIZE_ALONG_1_DIM
    static const int size1Dim =  (2*(2*(FMB_Info_ws)+1) +1);
    // HALF_INTERACTION_LIST_SIZE_ALONG_1_DIM
    static const int halphSize1Dim =  (2*(FMB_Info_ws)+1);

    // EXPANSION_SIZE(FMB_Info.P)
    static const int FMB_Info_exp_size = ((FMB_Info_P)+1) * ((FMB_Info_P)+2) * 0.5;
    // NEXP_SIZE(FMB_Info.P)
    static const int FMB_Info_nexp_size = (FMB_Info_P + 1) * (FMB_Info_P + 1);

    // Level of the tree
    int treeHeight;
    // Width of the box at the root level
    double treeWidthAtRoot;

    // transfer_M2M_container
    FComplexe*** transitionM2M;
    // transfer_L2L_container
    FComplexe*** transitionL2L;

    // transfer_container
    FComplexe***** transferM2L;

    //[OK] spherical_harmonic_Outer_coefficients_array
    double sphereHarmoOuterCoef[FMB_Info_exp_size];
    //[OK] spherical_harmonic_Inner_coefficients_array
    double sphereHarmoInnerCoef[FMB_Info_exp_size];

    FComplexe current_thread_Y[FMB_Info_exp_size];

    // p_Y_theta_derivated
    FComplexe current_thread_Y_theta_derivated[FMB_Info_exp_size];

    // [OK?] p_precomputed_cos_and_sin_array
    FComplexe cosSin[FMB_Info_M2L_P + 1];
    // [OK?] p_associated_Legendre_function_Array
    double legendre[FMB_Info_M2L_exp_size];

    // pow_of_I_array
    static const double PiArrayInner[4];
    // pow_of_O_array
    static const double PiArrayOuter[4];

    // To store spherical position
    struct Spherical {
        double r, cosTheta, sinTheta, phi;
    };

    int expansion_Redirection_array_for_j[FMB_Info_M2L_P + 1 ];

    F3DPosition force_sum;
    double potential_sum;

    bool FMB_Info_up_to_P_in_M2L;

    //////////////////////////////////////////////////////////////////
    // Allocation
    //////////////////////////////////////////////////////////////////

    void expansion_Redirection_array_for_j_Initialize() {
            for( int h = 0; h <= FMB_Info_M2L_P ; ++h ){
                expansion_Redirection_array_for_j[h] = static_cast<int>( h * ( h + 1 ) * 0.5 );
            }
    }

    //spherical_harmonic_Outer_and_Inner_coefficients_array_Initialize
    void sphericalHarmonicInitialize(){
        // Outer coefficients:
        //std::cout << "sphereHarmoOuterCoef\n";
        double factOuter = 1.0;
        // in FMB code stoped at <= FMB_Info_M2L_P but this is not sufficient
        for(int idxP = 0 ; idxP < FMB_Info_exp_size; factOuter *= (++idxP) ){
            this->sphereHarmoOuterCoef[idxP] = factOuter;
            //std::cout << this->sphereHarmoOuterCoef[idxl] << "\n";
        }

        // Inner coefficients:
        double* currentInner = this->sphereHarmoInnerCoef;
        double factInner = 1.0;
        double powN1idxP = 1.0;
        //std::cout << "sphereHarmoInnerCoef\n";
        for(int idxP = 0 ; idxP <= this->FMB_Info_M2L_P ; factInner *= (++idxP), powN1idxP = -powN1idxP){
            for(int idxMP = 0, fact_l_m = factInner; idxMP <= idxP ; fact_l_m *= idxP+(++idxMP), ++currentInner){
                *currentInner = powN1idxP / fact_l_m;
                //std::cout << (*currentInner) << "\n";
            }
        }

        //for(int temp = 0 ; temp < 6 ; ++temp){
        //    std::cout << this->sphereHarmoInnerCoef[temp] << "\n";
        //}
    }


    // transfer_L2L_Allocate
    // transfer_M2M_Allocate
    // transfer_M2L_Allocate
    void transferAllocate(){
        // M2M L2L
        this->transitionM2M = new FComplexe**[this->treeHeight];
        this->transitionL2L = new FComplexe**[this->treeHeight];

        for(int idxLevel = 0 ; idxLevel < this->treeHeight ; ++idxLevel){

            this->transitionM2M[idxLevel] = new FComplexe*[8];
            this->transitionL2L[idxLevel] = new FComplexe*[8];

            for(long idxChild = 0; idxChild < 8; ++idxChild){
                this->transitionM2M[idxLevel][idxChild] = new FComplexe[FMB_Info_exp_size];
                this->transitionL2L[idxLevel][idxChild] = new FComplexe[FMB_Info_exp_size];
            }
        }
        // M2L
        this->transferM2L = new FComplexe****[this->treeHeight+1];

        for(int idxLevel = 0; idxLevel < this->treeHeight; ++idxLevel){
            this->transferM2L[idxLevel] = new FComplexe***[this->size1Dim];

            for(long idxD1 = 0 ; idxD1 < this->size1Dim; ++idxD1){
                this->transferM2L[idxLevel][idxD1] = new FComplexe**[this->size1Dim];

                for(long idxD2 = 0; idxD2 < this->size1Dim; ++idxD2){
                    this->transferM2L[idxLevel][idxD1][idxD2] = new FComplexe*[this->size1Dim];

                    for(long idxD3 = 0; idxD3 < this->size1Dim; ++idxD3){
                        const int x = idxD1 - this->halphSize1Dim;
                        const int y = idxD2 - this->halphSize1Dim;
                        const int z = idxD3 - this->halphSize1Dim;

                        if( ( x*x + y*y + z*z ) >= ( 3*this->FMB_Info_ws*this->FMB_Info_ws + 0.1 ) ){
                            //[Blas] this->transferM2L[idxLevel][idxD1][idxD2][idxD3] = new FComplexe[this->FMB_Info_exp_size * this->FMB_Info_nexp_size];
                            this->transferM2L[idxLevel][idxD1][idxD2][idxD3] = new FComplexe[this->FMB_Info_M2L_exp_size];
                        }
                        else {
                            this->transferM2L[idxLevel][idxD1][idxD2][idxD3] = NULL;
                        }
                    }
                }
            }
        }
    }
    // transfer_L2L_free
    // transfer_M2M_free
    // transfer_M2L_free
    void transferDeallocate(){
        // M2M L2L
        for(long idxLevel = 0 ; idxLevel < this->treeHeight ; ++idxLevel){
            for(long idxChild = 0; idxChild < 8; ++idxChild){
                delete [] this->transitionM2M[idxLevel][idxChild];
                delete [] this->transitionL2L[idxLevel][idxChild];
            }
            delete [] this->transitionM2M[idxLevel];
            delete [] transitionL2L[idxLevel];
        }
        delete [] this->transitionM2M;
        delete [] this->transitionL2L;
        // M2L
        for(int idxLevel = 0 ; idxLevel < this->treeHeight; ++idxLevel){
            for(long idxD1 = 0 ; idxD1 < this->size1Dim ; ++idxD1){
                for(long idxD2 = 0 ; idxD2 < this->size1Dim ; ++idxD2){
                    for(long idxD3 = 0 ; idxD3 < this->size1Dim; ++idxD3){
                        delete [] this->transferM2L[idxLevel][idxD1][idxD2][idxD3];
                    }
                    delete [] this->transferM2L[idxLevel][idxD1][idxD2];
                }
                delete [] this->transferM2L[idxLevel][idxD1];
            }
            delete [] this->transferM2L[idxLevel];
        }
        delete [] this->transferM2L;
    }

    //////////////////////////////////////////////////////////////////
    // Utils
    //////////////////////////////////////////////////////////////////

    // position_2_r_cos_th_sin_th_ph
    void positionToSphere(const F3DPosition& inVector, Spherical* const outSphere ){
        const double x2y2 = (inVector.getX() * inVector.getX()) + (inVector.getY() * inVector.getY());

        outSphere->r = FMath::Sqrt( x2y2 + (inVector.getZ() * inVector.getZ()));
        outSphere->phi = FMath::Atan2(inVector.getY(),inVector.getX());
        outSphere->cosTheta = inVector.getZ() / outSphere->r; // cos_th = z/r
        outSphere->sinTheta = FMath::Sqrt(x2y2) / outSphere->r; // sin_th = sqrt(x^2 + y^2)/r
    }

    // associated_Legendre_function_Fill_complete_array_of_values_for_cos
    void legendreFunction( const double inCosTheta, const double inSinTheta, double* const outResults ){
        // l=0:         results[current++] = 1.0; // P_0^0(cosTheta) = 1
        int idxCurrent = 0;
        outResults[idxCurrent++] = 1.0;

        // l=1:
        // Compute P_1^{0} using (3) and store it into results_array:
        outResults[idxCurrent++] = inCosTheta;

        // Compute P_1^1 using (2 bis) and store it into results_array
        const double invSinTheta = -inSinTheta; // somx2 => -sinTheta
        outResults[idxCurrent++] = invSinTheta;

        // l>1:
        int idxCurrent1m = 1; //pointer on P_{l-1}^m P_1^0
        int idxCurrent2m = 0; //pointer on P_{l-2}^m P_0^0
        double fact = 3.0;

        // Remark: p_results_array_l_minus_1_m and p_results_array_l_minus_2_m
        // just need to be incremented at each iteration.
        for(int idxl = 2; idxl <= this->LMax ; ++idxl ){
                for( int idxm = 0; idxm <= idxl - 2 ; ++idxm , ++idxCurrent , ++idxCurrent1m , ++idxCurrent2m ){
                        // Compute P_l^m, l >= m+2, using (1) and store it into results_array:
                        outResults[idxCurrent] = (inCosTheta * ( 2 * idxl - 1 ) * outResults[idxCurrent1m] - ( idxl + idxm - 1 )
                                        * outResults[idxCurrent2m] ) / ( idxl - idxm );
                }
                // p_results_array_l_minus_1_m now points on P_{l-1}^{l-1}

                // Compute P_l^{l-1} using (3) and store it into ptrResults:
                outResults[idxCurrent++] = inCosTheta * ( 2 * idxl - 1 ) * outResults[idxCurrent1m];

                // Compute P_l^l using (2 bis) and store it into results_array:
                outResults[idxCurrent++] = fact * invSinTheta * outResults[idxCurrent1m];

                fact += 2.0;
                ++idxCurrent1m;
        }
    }

    // spherical_harmonic_Inner
    void harmonicInner(const Spherical& inSphere, FComplexe* const outResults){

        FComplexe* ptrCosSin = this->cosSin;
        for(int idxl = 0 , idxlMod4 = 0; idxl <= LMax ; ++idxl, ++idxlMod4, ++ptrCosSin){
            if(idxlMod4 == 4) idxlMod4 = 0;
            const double angle = idxl * inSphere.phi + this->PiArrayInner[idxlMod4];

            ptrCosSin->setReal( FMath::Sin(angle + FMath::FPiDiv2) );
            ptrCosSin->setImag( FMath::Sin(angle) );
            //std::cout<< idxl << "=" << ptrCosSin->getReal() << "/" << ptrCosSin->getImag() << " (" << idxl << "/" << inSphere.phi << "/" << PiArray[lMod4] << ")\n";
        }

        legendreFunction(inSphere.cosTheta, inSphere.sinTheta, this->legendre);
        //printf("FMB_Info_M2L_exp_size=%d\n",FMB_Info_M2L_exp_size);
        //for(int temp = 0 ; temp < FMB_Info_M2L_exp_size ; ++temp){
        //    printf("%f\n",this->legendre[temp]);
        //}

        FComplexe* currentResult = outResults;
        int idxLegendre = 0;//ptr_associated_Legendre_function_Array
        int idxSphereHarmoCoef = 0;
        double idxRl = 1.0 ;

        for(int idxl = 0; idxl <= this->LMax ; ++idxl, idxRl *= inSphere.r){
            for(int idxm = 0 ; idxm <= idxl ; ++idxm, ++currentResult, ++idxSphereHarmoCoef, ++idxLegendre){
                const double magnitude = this->sphereHarmoInnerCoef[idxSphereHarmoCoef] * idxRl * legendre[idxLegendre];
                currentResult->setReal( magnitude * this->cosSin[idxm].getReal() );
                currentResult->setImag( magnitude * this->cosSin[idxm].getImag() );
                //printf("magnitude=%f idxRl=%f sphereHarmoInnerCoef=%f real=%f imag=%f\n",magnitude,idxRl,this->sphereHarmoInnerCoef[idxSphereHarmoCoef],currentResult->getReal(),currentResult->getImag());
            }
        }

    }
    // spherical_harmonic_Outer
    void harmonicOuter(const Spherical& inSphere, FComplexe* const outResults){

        FComplexe* ptrCosSin = this->cosSin;
        for(int idxl = 0, idxlMod4 = 0; idxl <= LMax ; ++idxl, ++idxlMod4, ++ptrCosSin){
            if(idxlMod4 == 4) idxlMod4 = 0;
            const double angle = idxl * inSphere.phi + this->PiArrayOuter[idxlMod4];

            ptrCosSin->setReal( FMath::Sin(angle + FMath::FPiDiv2) );
            ptrCosSin->setImag( FMath::Sin(angle) );

            //printf("l=%d \t inSphere.phi=%f \t this->PiArrayOuter[idxlMod4]=%f \t angle=%f \t FMath::Sin(angle + FMath::FPiDiv2)=%f \t FMath::Sin(angle)=%f\n",
            //        idxl, inSphere.phi, this->PiArrayOuter[idxlMod4], angle, FMath::Sin(angle + FMath::FPiDiv2) , FMath::Sin(angle));
        }

        legendreFunction(inSphere.cosTheta, inSphere.sinTheta, this->legendre);

        int idxLegendre = 0;
        FComplexe* currentResult = outResults;

        const double invR = 1/inSphere.r;
        double idxRl1 = invR;
        for(int idxl = 0 ; idxl <= this->LMax ; ++idxl, idxRl1 *= invR){
            for(int idxm = 0 ; idxm <= idxl ; ++idxm, ++currentResult, ++idxLegendre){
                const double magnitude = this->sphereHarmoOuterCoef[idxl-idxm] * idxRl1 * this->legendre[idxLegendre];
                currentResult->setReal( magnitude * this->cosSin[idxm].getReal() );
                currentResult->setImag( magnitude * this->cosSin[idxm].getImag() );
                //printf("l=%d\t m=%d\t idxRl1=%f\t magnitude=%f\n",idxl,idxm,idxRl1,magnitude);
                //printf("l=%d\t m=%d\t this->cosSin[idxm].getReal()=%f\t this->cosSin[idxm].getImag()=%f\n",
                //       idxl,idxm,this->cosSin[idxm].getReal(),this->cosSin[idxm].getImag());
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
    void harmonicInnerThetaDerivated(
                         const Spherical& inSphere,
                         FComplexe * results_array,
                         FComplexe * theta_derivated_results_array
                         ){
        FComplexe *p_term = results_array;
        FComplexe *p_theta_derivated_term = theta_derivated_results_array;
        double *p_spherical_harmonic_Inner_coefficients_array = sphereHarmoInnerCoef;
        double *ptr_associated_Legendre_function_Array = legendre;
        double *start_ptr_associated_Legendre_function_Array = ptr_associated_Legendre_function_Array;


        // Initialization of precomputed_cos_and_sin_array:
        for(int idxm = 0 , idxmMod4 = 0; idxm <= FMB_Info_P ; ++idxm, ++idxmMod4){
            if(idxmMod4 == 4) idxmMod4 = 0;
            const double angle = idxm*inSphere.phi + PiArrayInner[idxmMod4];
            (cosSin + idxm)->setReal(FMath::Sin(angle + FMath::FPiDiv2));
            (cosSin + idxm)->setImag(FMath::Sin(angle));

            //printf("l=%d \t inSphere.phi=%f \t this->PiArrayOuter[idxlMod4]=%f \t angle=%f \t FMath::Sin(angle + FMath::FPiDiv2)=%f \t FMath::Sin(angle)=%f\n",
            //        idxm, inSphere.phi, this->PiArrayInner[idxmMod4], angle, FMath::Sin(angle + FMath::FPiDiv2) , FMath::Sin(angle));
        }

        // Initialization of associated_Legendre_function_Array:
        legendreFunction( inSphere.cosTheta, inSphere.sinTheta, this->legendre);

        // r^l
        double r_l = 1.0;
        for (int l = 0 ; l <= FMB_Info_P ; ++l, r_l *= inSphere.r){
            double magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * r_l * (*ptr_associated_Legendre_function_Array);
            // m<l:
            int m = 0;
            for(; m < l ; ++m, ++p_term, ++p_theta_derivated_term, ++p_spherical_harmonic_Inner_coefficients_array, ++ptr_associated_Legendre_function_Array){
                // Computation of Inner_l^m(r, theta, phi):
                p_term->setReal( magnitude * (cosSin + m)->getReal());
                p_term->setImag( magnitude * (cosSin + m)->getImag());
                // Computation of {\partial Inner_l^m(r, theta, phi)}/{\partial theta}:
                magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * r_l * ((l*inSphere.cosTheta*(*ptr_associated_Legendre_function_Array)
                                        - (l+m)*(*(start_ptr_associated_Legendre_function_Array + expansion_Redirection_array_for_j[l-1] + m))) / inSphere.sinTheta);
                p_theta_derivated_term->setReal(magnitude * (cosSin + m)->getReal());
                p_theta_derivated_term->setImag(magnitude * (cosSin + m)->getImag());

                //printf("magnitude=%f r_l=%f p_spherical_harmonic_Inner_coefficients_array=%f real=%f imag=%f\n",
                //       magnitude,r_l,*p_spherical_harmonic_Inner_coefficients_array,p_theta_derivated_term->getReal(),p_theta_derivated_term->getImag());
            }

            // m=l:
            // Computation of Inner_m^m(r, theta, phi):
            magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * r_l * (*ptr_associated_Legendre_function_Array);
            p_term->setReal(magnitude * (cosSin + m)->getReal());
            p_term->setImag(magnitude * (cosSin + m)->getImag());
            // Computation of {\partial Inner_m^m(r, theta, phi)}/{\partial theta}:
            magnitude = (*p_spherical_harmonic_Inner_coefficients_array) * r_l * (m * inSphere.cosTheta * (*ptr_associated_Legendre_function_Array) / inSphere.sinTheta);
            p_theta_derivated_term->setReal(magnitude * (cosSin + m)->getReal());
            p_theta_derivated_term->setImag(magnitude * (cosSin + m)->getImag());

            //printf("magnitude=%f r_l=%f p_spherical_harmonic_Inner_coefficients_array=%f real=%f imag=%f\n",
            //       magnitude,r_l,*p_spherical_harmonic_Inner_coefficients_array,p_theta_derivated_term->getReal(),p_theta_derivated_term->getImag());

            ++p_term;
            ++p_theta_derivated_term;
            ++p_spherical_harmonic_Inner_coefficients_array;
            ++ptr_associated_Legendre_function_Array;
        }
    }

    //[fmb] ff_matrix_Convert_exp_2_transfer_M2L_matrix
    /*[Blas] void expffmatrix(const FComplexe* transfer_exp, FComplexe* const ff_transfer_matrix){
        FComplexe *p_ff_transfer_matrix = ff_transfer_matrix;

        for( int M = 0; M <= this->FMB_Info_P ; ++M ){
            for(int  m = 0; m <= M ; ++m ){
                for(int  N = 0; N <= FMB_Info_P ; ++N ){
                    for(int n = 0; n <= 2 * N ; ++n , ++p_ff_transfer_matrix ){
                        const int k = N - n - m;
                        if(k < 0){
                            const int pow_of_minus_1_k = ((k%2) ? -1 : 1);
                            p_ff_transfer_matrix->setReal(pow_of_minus_1_k  * transfer_exp[ expansion_Redirection_array_for_j[M + N] - k].getReal());
                            p_ff_transfer_matrix->setImag((-pow_of_minus_1_k) * transfer_exp[ expansion_Redirection_array_for_j[M + N] - k].getImag());
                        }
                        else{
                            int temp = expansion_Redirection_array_for_j[M + N];
                            *p_ff_transfer_matrix = transfer_exp[ expansion_Redirection_array_for_j[M + N] + k];
                        }
                    }
                }
            }
        }
    }*/

    //////////////////////////////////////////////////////////////////
    // Precompute
    //////////////////////////////////////////////////////////////////

    // transfer_M2M_Precompute_all_levels
    // transfer_L2L_Precompute_all_levels
    void precomputeM2M(){
        double treeWidthAtLevel = this->treeWidthAtRoot/2;

        for(int idxLevel = 0 ; idxLevel < this->treeHeight ; ++idxLevel ){
            const F3DPosition father(treeWidthAtLevel,treeWidthAtLevel,treeWidthAtLevel);
            treeWidthAtLevel /= 2;

            //std::cout << "[precomputeM2M]treeWidthAtLevel=" << treeWidthAtLevel << "\n";
            //printf("\tidxLevel=%d\tFather.x=%f\tFather.y=%f\tFather.z=%f\n",idxLevel,father.getX(),father.getY(),father.getZ());

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild ){
                FTreeCoordinate childBox;
                childBox.setPositionFromMorton(idxChild,1);

                const F3DPosition M2MVector (
                        father.getX() - (treeWidthAtLevel * (1 + (childBox.getX() * 2))),
                        father.getY() - (treeWidthAtLevel * (1 + (childBox.getY() * 2))),
                        father.getZ() - (treeWidthAtLevel * (1 + (childBox.getZ() * 2)))
                );
                Spherical sphericalM2M;
                positionToSphere(M2MVector,&sphericalM2M);
                harmonicInner(sphericalM2M,this->transitionM2M[idxLevel][idxChild]);

                const F3DPosition L2LVector (
                        (treeWidthAtLevel * (1 + (childBox.getX() * 2))) - father.getX(),
                        (treeWidthAtLevel * (1 + (childBox.getY() * 2))) - father.getY(),
                        (treeWidthAtLevel * (1 + (childBox.getZ() * 2))) - father.getZ()
                );
                Spherical sphericalL2L;
                positionToSphere(L2LVector,&sphericalL2L);
                harmonicInner(sphericalL2L,this->transitionL2L[idxLevel][idxChild]);

                //std::cout << "[M2M_vector]" << idxLevel << "/" << idxChild << " = " << M2MVector.getX() << "/" << M2MVector.getY() << "/" << M2MVector.getZ() << "\n";
                //std::cout << "[M2M_vectorSpherical]" << idxLevel << "/" << idxChild << " = " << sphericalM2M.r << "/" << sphericalM2M.cosTheta << "/" << sphericalM2M.sinTheta << "/" << sphericalM2M.phi << "\n";
                //for(int idxExpSize = 0 ; idxExpSize < FMB_Info_exp_size ; ++idxExpSize){
                    //std::cout << "transitionL2L[" << idxLevel << "][" << idxChild << "][" << idxExpSize << "]=" << this->transitionL2L[idxLevel][idxChild][idxExpSize].getReal()<<"/"<<this->transitionL2L[idxLevel][idxChild][idxExpSize].getImag()<< "\n";
                    //std::cout << "transitionM2M[" << idxLevel << "][" << idxChild << "][" << idxExpSize << "]=" << this->transitionM2M[idxLevel][idxChild][idxExpSize].getReal()<<"/"<<this->transitionM2M[idxLevel][idxChild][idxExpSize].getImag()<< "\n";
                //}
            }
        }

    }


    // transfer_M2L_Precompute_all_levels
    void precomputeM2L(){
        //[Blas] FComplexe tempComplexe[FMB_Info_M2L_exp_size];
        //printf("FMB_Info.M2L_exp_size = %d\n",FMB_Info_M2L_exp_size);

        double treeWidthAtLevel = this->treeWidthAtRoot;
        for(int idxLevel = 0 ; idxLevel < this->treeHeight ; ++idxLevel ){
            //printf("level = %d \t width = %lf\n",idxLevel,treeWidthAtLevel);
            for( int idxd1 = 0; idxd1 < this->size1Dim ; ++idxd1 ){

                for( int idxd2 = 0; idxd2 < this->size1Dim ; ++idxd2 ){

                    for( int idxd3 = 0; idxd3 < this->size1Dim ; ++idxd3 ){
                        const long x = idxd1 - this->halphSize1Dim;
                        const long y = idxd2 - this->halphSize1Dim;
                        const long z = idxd3 - this->halphSize1Dim;

                        //printf("x=%ld \t y=%ld \t z=%ld\n",x,y,z);

                        if( ( x*x + y*y + z*z ) >= ( 3*FMB_Info_ws*FMB_Info_ws + 0.1 ) ){
                            const F3DPosition relativePos( x*treeWidthAtLevel , y*treeWidthAtLevel , z*treeWidthAtLevel );                            

                            Spherical spherical;
                            positionToSphere(relativePos,&spherical);
                            // Not blas so
                            //printf("transferM2L[%d][%d][%d][%d]\n", idxLevel, idxd1, idxd2, idxd3);
                            harmonicOuter(spherical,this->transferM2L[idxLevel][idxd1][idxd2][idxd3]);
                            //for(int idxTemp = 0 ; idxTemp < this->FMB_Info_M2L_exp_size ; ++idxTemp){
                            //    printf("transferM2L[%d][%d][%d][%d][%d]=%f/%f\n", idxLevel, idxd1, idxd2, idxd3, idxTemp, this->transferM2L[idxLevel][idxd1][idxd2][idxd3][idxTemp].getReal(),this->transferM2L[idxLevel][idxd1][idxd2][idxd3][idxTemp].getImag());
                            //}

                            //[Blas] harmonicOuter(spherical,tempComplexe);
                            //[Blas] ff_matrix_Convert_exp_2_transfer_M2L_matrix
                            //[Blas] expffmatrix( tempComplexe , this->transferM2L[idxLevel][idxd1][idxd2][idxd3] );
                        }

                    }
                }
            }
            treeWidthAtLevel /= 2;
        }
    }

    void buildPrecompute(){
        expansion_Redirection_array_for_j_Initialize();
        sphericalHarmonicInitialize();
        transferAllocate();

        precomputeM2M();
        precomputeM2L();
    }

public:
    FAbstractFmbKernels(const int inTreeHeight, const double inTreeWidth)
        : treeHeight(inTreeHeight), treeWidthAtRoot(inTreeWidth),
          transitionM2M(0), transitionL2L(0), potential_sum(0), FMB_Info_up_to_P_in_M2L(true) {
        buildPrecompute();
    }

    FAbstractFmbKernels(const FAbstractFmbKernels& other)
        : treeHeight(other.treeHeight), treeWidthAtRoot(other.treeWidthAtRoot),
          transitionM2M(0), transitionL2L(0), potential_sum(0), FMB_Info_up_to_P_in_M2L(true) {
        buildPrecompute();
    }

    virtual void init(){
        force_sum = F3DPosition();
        potential_sum = 0;
    }

    /** Default destructor */
    virtual ~FAbstractFmbKernels(){
        transferDeallocate();
    }

    F3DPosition getForcesSum() const {
        return force_sum;
    }

    double getPotential() const{
        return potential_sum;
    }

    /////////////////////////////////////////////////////////////////////////////////
    //    Upward
    /////////////////////////////////////////////////////////////////////////////////

    /** OK!
    * expansion_P2M_add
    * Multipole expansion with m charges q_i in Q_i=(rho_i, alpha_i, beta_i)
    *whose relative coordinates according to *p_center are:
    *Q_i - *p_center = (rho'_i, alpha'_i, beta'_i);
    *
    *For j=0..P, k=-j..j, we have:
    *
    *M_j^k = (-1)^j { sum{i=1..m} q_i Inner_j^k(rho'_i, alpha'_i, beta'_i) }
    *
    *However the extern loop is over the bodies (i=1..m) in our code and as an
    *intern loop we have: j=0..P, k=-j..j
    *
    *and the potential is then given by:
    *
    * Phi(x) = sum_{n=0}^{+} sum_{m=-n}^{n} M_n^m O_n^{-m} (x - *p_center)
    *
    */
    void P2M(CellClass* const inPole, const FList<ParticuleClass*>* const inParticules) {

        for(typename FList<ParticuleClass*>::ConstBasicIterator iterParticule(*inParticules);
                                iterParticule.isValide() ; iterParticule.progress()){

            Spherical spherical;
            positionToSphere(iterParticule.value()->getPosition() - inPole->getPosition(), &spherical);

            //std::cout << "Working on part " << iterParticule.value()->getValue() << "\n";
            F3DPosition tempPos = iterParticule.value()->getPosition() - inPole->getPosition();
            //ok printf("\tpos_rel.x=%f\tpos_rel.y=%f\tpos_rel.z=%f\n",tempPos.getX(),tempPos.getY(),tempPos.getZ());
            //ok printf("\tp_center.x=%f\tp_center.y=%f\tp_center.z=%f\n",inPole->getPosition().getX(),inPole->getPosition().getY(),inPole->getPosition().getZ());
            //ok printf("\tbody.x=%f\tbody.y=%f\tbody.z=%f\n",iterParticule.value()->getPosition().getX(),iterParticule.value()->getPosition().getY(),iterParticule.value()->getPosition().getZ());

            harmonicInner(spherical,current_thread_Y);

            //ok printf("\tr=%f\tcos_theta=%f\tsin_theta=%f\tphi=%f\n",spherical.r,spherical.cosTheta,spherical.sinTheta,spherical.phi);

            FComplexe* p_exp_term = inPole->getMultipole();
            FComplexe* p_Y_term = current_thread_Y;
            double pow_of_minus_1_j = 1.0;//(-1)^j
            const double valueParticule = iterParticule.value()->getValue();

            for(int j = 0 ; j <= FMB_Info_P ; ++j, pow_of_minus_1_j = -pow_of_minus_1_j ){
                for(int k = 0 ; k <= j ; ++k, ++p_Y_term, ++p_exp_term){
                    p_Y_term->mulRealAndImag( valueParticule * pow_of_minus_1_j );
                    (*p_exp_term) += (*p_Y_term);
                    //printf("\tj=%d\tk=%d\tp_exp_term.real=%f\tp_exp_term.imag=%f\tp_Y_term.real=%f\tp_Y_term.imag=%f\tpow_of_minus_1_j=%f\n", j,k,(*p_exp_term).getReal(),(*p_exp_term).getImag(),(*p_Y_term).getReal(),(*p_Y_term).getImag(),pow_of_minus_1_j);
                }
            }
        }
    }

    /**
    *-----------------------------------
    *octree_Upward_pass_internal_cell
    *expansion_M2M_add
    *-----------------------------------
    *We compute the translation of multipole_exp_src from *p_center_of_exp_src to
    *p_center_of_exp_target, and add the result to multipole_exp_target.
    *
    * O_n^l (with n=0..P, l=-n..n) being the former multipole expansion terms
    * (whose center is *p_center_of_multipole_exp_src) we have for the new multipole
    * expansion terms (whose center is *p_center_of_multipole_exp_target):

    * M_j^k = sum{n=0..j}
    * sum{l=-n..n, |k-l|<=j-n}
    * O_n^l Inner_{j-n}^{k-l}(rho, alpha, beta)
    *
    * where (rho, alpha, beta) are the spherical coordinates of the vector :
    * p_center_of_multipole_exp_target - *p_center_of_multipole_exp_src
    *
    * Warning: if j-n < |k-l| we do nothing.
     */
    void M2M(CellClass* const inPole, const CellClass*const* const inChild, const int inLevel) {

        // We do NOT have: for(l=n-j+k; l<=j-n+k ;++l){} <=> for(l=-n; l<=n ;++l){if (j-n >= abs(k-l)){}}
        //     But we have:  for(k=MAX(0,n-j+l); k<=j-n+l; ++k){} <=> for(k=0; k<=j; ++k){if (j-n >= abs(k-l)){}}
        //     (This is not the same as in L2L since the range of values of k is not the same, compared to "n-j" or "j-n".)
        //     Therefore the loops over n and l are the outmost ones and
        //     we invert the loop over j with the summation with n:
        //     for{j=0..P} sum_{n=0}^j <-> sum_{n=0}^P for{j=n..P}
        FComplexe* const multipole_exp_target = inPole->getMultipole();

        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            if(!inChild[idxChild]) continue;
            //printf("\tChild %d\n",idxChild);

            const FComplexe* const multipole_exp_src = inChild[idxChild]->getMultipole();

            const FComplexe* const M2M_transfer = transitionM2M[inLevel][idxChild];

            for(int n = 0 ; n <= FMB_Info_P ; ++n ){
                // l<0 // (-1)^l
                double pow_of_minus_1_for_l = ( n % 2 ? -1 : 1);

                // O_n^l : here points on the source multipole expansion term of degree n and order |l|
                const FComplexe* p_src_exp_term = multipole_exp_src + expansion_Redirection_array_for_j[n]+n;
                //printf("\t[p_src_exp_term] expansion_Redirection_array_for_j[n]=%d\tn=%d\n",expansion_Redirection_array_for_j[n],n);

                int l = -n;
                for(; l<0 ; ++l, --p_src_exp_term, pow_of_minus_1_for_l = -pow_of_minus_1_for_l){

                    for(int j = n ; j<= FMB_Info_P ; ++j ){
                        // M_j^k
                        FComplexe *p_target_exp_term = multipole_exp_target + expansion_Redirection_array_for_j[j];
                        //printf("\t[p_target_exp_term] expansion_Redirection_array_for_j[j]=%d\n",expansion_Redirection_array_for_j[j]);
                        // Inner_{j-n}^{k-l} : here points on the M2M transfer function/expansion term of degree n-j and order |k-l|
                        const FComplexe *p_Inner_term= M2M_transfer + expansion_Redirection_array_for_j[j-n]-l /* k==0 */;
                        //printf("\t[p_Inner_term] expansion_Redirection_array_for_j[j-n]=%d\tl=%d\n",expansion_Redirection_array_for_j[j-n],-l);

                        // since n-j+l<0
                        for(int k=0 ; k <= (j-n+l) ; ++k, ++p_target_exp_term, ++p_Inner_term){ // l<0 && k>=0 => k-l>0
                            p_target_exp_term->incReal( pow_of_minus_1_for_l *
                                                        ((p_src_exp_term->getReal() * p_Inner_term->getReal()) +
                                                         (p_src_exp_term->getImag() * p_Inner_term->getImag())));
                            p_target_exp_term->incImag( pow_of_minus_1_for_l *
                                                        ((p_src_exp_term->getReal() * p_Inner_term->getImag()) -
                                                         (p_src_exp_term->getImag() * p_Inner_term->getReal())));

                            //printf("\tp_src_exp_term->getReal()=%f\tp_src_exp_term->getImag()=%f\n", p_src_exp_term->getReal(),p_src_exp_term->getImag());
                            //printf("\tp_Inner_term->getReal()=%f\tp_Inner_term->getImag()=%f\n", p_Inner_term->getReal(),p_Inner_term->getImag());
                            //printf("\tn[%d]l[%d]j[%d]k[%d] = %f / %f\n",n,l,j,k,p_target_exp_term->getReal(),p_target_exp_term->getImag());

                        } // for k
                    } // for j
                } // for l

                // l>=0
                for(; l <= n ; ++l, ++p_src_exp_term, pow_of_minus_1_for_l = -pow_of_minus_1_for_l){

                    for( int j=n ; j <= FMB_Info_P ; ++j ){
                        // (-1)^k
                        double pow_of_minus_1_for_k = ( FMath::Max(0,n-j+l) %2 ? -1 : 1 );
                        // M_j^k
                        FComplexe *p_target_exp_term = multipole_exp_target + expansion_Redirection_array_for_j[j] + FMath::Max(0,n-j+l);
                        // Inner_{j-n}^{k-l} : here points on the M2M transfer function/expansion term of degree n-j and order |k-l|
                        const FComplexe *p_Inner_term = M2M_transfer + expansion_Redirection_array_for_j[j-n] + l - FMath::Max(0,n-j+l);// -(k-l)

                        int k = FMath::Max(0,n-j+l);
                        for(; k <= (j-n+l) && (k-l) < 0 ; ++k, ++p_target_exp_term, --p_Inner_term, pow_of_minus_1_for_k = -pow_of_minus_1_for_k){ /* l>=0 && k-l<0 */
                            p_target_exp_term->incReal( pow_of_minus_1_for_k * pow_of_minus_1_for_l *
                                                        ((p_src_exp_term->getReal() * p_Inner_term->getReal()) +
                                                         (p_src_exp_term->getImag() * p_Inner_term->getImag())));
                            p_target_exp_term->incImag(pow_of_minus_1_for_k * pow_of_minus_1_for_l *
                                                       ((p_src_exp_term->getImag() * p_Inner_term->getReal()) -
                                                        (p_src_exp_term->getReal() * p_Inner_term->getImag())));
                            //printf("\tn[%d]l[%d]j[%d]k[%d] = %f / %f\n",n,l,j,k,p_target_exp_term->getReal(),p_target_exp_term->getImag());

                        } // for k

                        for(; k <= (j - n + l) ; ++k, ++p_target_exp_term, ++p_Inner_term){ // l>=0 && k-l>=0
                            p_target_exp_term->incReal(
                                    (p_src_exp_term->getReal() * p_Inner_term->getReal()) -
                                    (p_src_exp_term->getImag() * p_Inner_term->getImag()));
                            p_target_exp_term->incImag(
                                    (p_src_exp_term->getImag() * p_Inner_term->getReal()) +
                                    (p_src_exp_term->getReal() * p_Inner_term->getImag()));
                            //printf("\tn[%d]l[%d]j[%d]k[%d] = %f / %f\n",n,l,j,k,p_target_exp_term->getReal(),p_target_exp_term->getImag());

                        } // for k
                    } // for j
                } // for l
            } // for n
        }
    }

    /////////////////////////////////////////////////////////////////////////////////
    //    M2L
    /////////////////////////////////////////////////////////////////////////////////

    virtual void changeProgression(int*const start_for_j , FComplexe** const p_target_exp_term) = 0;

    /**
    *------------------
    * expansion_M2L_add
    *-------------------
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
    void M2L(CellClass* const pole, const CellClass*const* const distantNeighbors, const int size, const int inLevel) {
        FTreeCoordinate coordCenter;
        coordCenter.setPositionFromMorton(pole->getMortonIndex(),inLevel);

        for(int idxSize = 0 ; idxSize < size ; ++idxSize){
            FTreeCoordinate coordNeighbors;
            coordNeighbors.setPositionFromMorton(distantNeighbors[idxSize]->getMortonIndex(),inLevel);

           // printf("Morton = %lld\n",pole->getMortonIndex());
            //printf("\tMorton Neighbors = %lld\n",distantNeighbors[idxSize]->getMortonIndex());
            //printf("\tidxSize = %d\tleve = %d\tMorton = %lld\n",idxSize,inLevel,distantNeighbors[idxSize]->getMortonIndex());

            const FComplexe* const M2L_transfer = transferM2L[inLevel]
                                      [(coordCenter.getX()+halphSize1Dim-coordNeighbors.getX())]
                                      [(coordCenter.getY()+halphSize1Dim-coordNeighbors.getY())]
                                      [(coordCenter.getZ()+halphSize1Dim-coordNeighbors.getZ())];
            /*printf("level = %d\tx=%ld\ty=%ld\tz=%ld\n", inLevel,
                   (coordCenter.getX()-coordNeighbors.getX()),
                   (coordCenter.getY()-coordNeighbors.getY()),
                   (coordCenter.getZ()-coordNeighbors.getZ()));*/
            /*printf("M2L_transfer[0]= %f/%f\n",M2L_transfer->getReal(),M2L_transfer->getImag());
            printf("M2L_transfer[1]= %f/%f\n",M2L_transfer[1].getReal(),M2L_transfer[1].getImag());
            printf("M2L_transfer[2]= %f/%f\n",M2L_transfer[2].getReal(),M2L_transfer[2].getImag());*/
            const FComplexe* const multipole_exp_src = distantNeighbors[idxSize]->getMultipole();
            // L_j^k
            FComplexe* p_target_exp_term = pole->getLocal();
            int start_for_j = 0;

            changeProgression(&start_for_j,&p_target_exp_term);

            //    HPMSTART(51, "M2L computation (loops)");
            for (int j = start_for_j ; j <= FMB_Info_P ; ++j){

                int stop_for_n = FMB_Info_P;
                if (FMB_Info_up_to_P_in_M2L){
                    stop_for_n = FMB_Info_P - j;
                }

                // (-1)^k
                double pow_of_minus_1_for_k = 1.0;
                for (int k = 0 ; k <= j ; ++k, pow_of_minus_1_for_k = -pow_of_minus_1_for_k, ++p_target_exp_term){

                    // (-1)^n
                    double pow_of_minus_1_for_n = 1.0;
                    for (int n = 0 ; n <= stop_for_n ; ++n, pow_of_minus_1_for_n = -pow_of_minus_1_for_n){

                        // O_n^l : here points on the source multipole expansion term of degree n and order |l|
                        const FComplexe *p_src_exp_term = multipole_exp_src + expansion_Redirection_array_for_j[n] + n;
                        // Outer_{j+n}^{-k-l} : here points on the M2L transfer function/expansion term of degree j+n and order |-k-l|
                        const FComplexe *p_Outer_term = M2L_transfer + expansion_Redirection_array_for_j[n+j] + k+n;
                        //printf("expansion_Get_p_term(M2L_transfer, n+j, k+n)-M2L_transfer = %d \t(n=%d)\n",
                        //       p_Outer_term-M2L_transfer , n );
                        //printf("TRUC = %d\n",expansion_Redirection_array_for_j[n+j] + k+n);
                        double pow_of_minus_1_for_l = pow_of_minus_1_for_n; // (-1)^l
                        // We start with l=n (and not l=-n) so that we always set p_Outer_term to a correct value in the first loop.
                        int l=n;
                        for ( ; l>0 ; --l, pow_of_minus_1_for_l = -pow_of_minus_1_for_l, --p_src_exp_term, --p_Outer_term){ // we have -k-l<0 and l>0
                            p_target_exp_term->incReal( pow_of_minus_1_for_l * pow_of_minus_1_for_k *
                                                        ((p_src_exp_term->getReal() * p_Outer_term->getReal()) +
                                                         (p_src_exp_term->getImag() * p_Outer_term->getImag())));
                            p_target_exp_term->incImag( pow_of_minus_1_for_l * pow_of_minus_1_for_k *
                                                        ((p_src_exp_term->getImag() * p_Outer_term->getReal()) -
                                                         (p_src_exp_term->getReal() * p_Outer_term->getImag())));
//                            printf("\t p_target_exp_term->real = %lf \t p_target_exp_term->imag = %lf \n",
//                                                           p_target_exp_term->getReal(),p_target_exp_term->getImag());
//                            printf("\t p_src_exp_term->real = %lf \t p_src_exp_term->imag = %lf \n",
//                                                           p_src_exp_term->getReal(),p_src_exp_term->getImag());
//                            printf("\t p_Outer_term->real = %f \t p_Outer_term->imag = %f \n",
//                                                           p_Outer_term->getReal(),p_Outer_term->getImag());
                            //printf("p_Outer_term-M2L_transfer = %d\n",
                            //                               p_Outer_term-M2L_transfer);
                        }

                        for (; l>=-n && -k-l<0 ; --l, pow_of_minus_1_for_l = -pow_of_minus_1_for_l, ++p_src_exp_term, --p_Outer_term){ // we have -k-l<0 and l<=0
                            p_target_exp_term->incReal( pow_of_minus_1_for_k *
                                                        ((p_src_exp_term->getReal() * p_Outer_term->getReal()) -
                                                         (p_src_exp_term->getImag() * p_Outer_term->getImag())));
                            p_target_exp_term->decImag(  pow_of_minus_1_for_k *
                                                         ((p_src_exp_term->getImag() * p_Outer_term->getReal()) +
                                                          (p_src_exp_term->getReal() * p_Outer_term->getImag())));
//                            printf("\t\t p_target_exp_term->real = %lf \t p_target_exp_term->imag = %lf \n",
//                                                           p_target_exp_term->getReal(),p_target_exp_term->getImag());
//                            printf("\t\t p_src_exp_term->real = %lf \t p_src_exp_term->imag = %lf \n",
//                                                           p_src_exp_term->getReal(),p_src_exp_term->getImag());
//                            printf("\t\t p_Outer_term->real = %f \t p_Outer_term->imag = %f \n",
//                                                           p_Outer_term->getReal(),p_Outer_term->getImag());
                        }

                        for (; l>=-n; --l, pow_of_minus_1_for_l = -pow_of_minus_1_for_l, ++p_src_exp_term, ++p_Outer_term){ // we have -k-l>=0 and l<=0
                            p_target_exp_term->incReal( pow_of_minus_1_for_l *
                                                        ((p_src_exp_term->getReal() * p_Outer_term->getReal()) +
                                                         (p_src_exp_term->getImag() * p_Outer_term->getImag())));
                            p_target_exp_term->incImag( pow_of_minus_1_for_l *
                                                        ((p_src_exp_term->getReal() * p_Outer_term->getImag()) -
                                                         (p_src_exp_term->getImag() * p_Outer_term->getReal())));
//                            printf("\t\t p_target_exp_term->real = %lf \t p_target_exp_term->imag = %lf \n",
//                                                           p_target_exp_term->getReal(),p_target_exp_term->getImag());
//                            printf("\t\t p_src_exp_term->real = %lf \t p_src_exp_term->imag = %lf \n",
//                                                           p_src_exp_term->getReal(),p_src_exp_term->getImag());
//                            printf("\t\t p_Outer_term->real = %f \t p_Outer_term->imag = %f \n",
//                                                           p_Outer_term->getReal(),p_Outer_term->getImag());
                        }
                        //printf("\tj=%d\tk=%d\tn=%d\tl=%d\n",j,k,n,l);
                        //printf("\t p_target_exp_term->real = %lf \t p_target_exp_term->imag = %lf \n",
                        //       p_target_exp_term->getReal(),p_target_exp_term->getImag());
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////
    //    Downard
    /////////////////////////////////////////////////////////////////////////////////

    /**
      *We compute the shift of local_exp_src from *p_center_of_exp_src to
      *p_center_of_exp_target, and set the result to local_exp_target.
      *
      *O_n^l (with n=0..P, l=-n..n) being the former local expansion terms
      *(whose center is *p_center_of_exp_src) we have for the new local
      *expansion terms (whose center is *p_center_of_exp_target):
      *
      *L_j^k = sum{n=j..P}
      *sum{l=-n..n}
      *O_n^l Inner_{n-j}^{l-k}(rho, alpha, beta)
      *
      *where (rho, alpha, beta) are the spherical coordinates of the vector :
      *p_center_of_exp_target - *p_center_of_exp_src
      *
      *Warning: if |l-k| > n-j, we do nothing.
      */
    void L2L(const CellClass* const pole, CellClass** const child, const int inLevel) {

        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            // if no child at this position
            if(!child[idxChild]) continue;

            const FComplexe* const L2L_tranfer = transitionL2L[inLevel][idxChild];
            const FComplexe* const local_exp_src = pole->getLocal();
            FComplexe* const local_exp_target = child[idxChild]->getLocal();

            //printf("Level %d\n", inLevel);
            //printf("Father morton %lld\n", pole->getMortonIndex());
            //printf("Child morton %lld\n", child[idxChild]->getMortonIndex());

            // L_j^k
            FComplexe* p_target_exp_term = local_exp_target;
            for (int j=0 ; j<= FMB_Info_P ; ++j){
                // (-1)^k
                double pow_of_minus_1_for_k = 1.0;
                for (int k=0 ; k <= j ; ++k, pow_of_minus_1_for_k = -pow_of_minus_1_for_k, ++p_target_exp_term){
                    for (int n=j; n<=FMB_Info_P;++n){
                        // O_n^l : here points on the source multipole expansion term of degree n and order |l|
                        const FComplexe* p_src_exp_term = local_exp_src + expansion_Redirection_array_for_j[n] + n-j+k;
                        //printf("expansion_Redirection_array_for_j[n] + n-j+k %d\n", expansion_Redirection_array_for_j[n] + n-j+k);
                        int l = n-j+k;
                        // Inner_{n-j}^{l-k} : here points on the L2L transfer function/expansion term of degree n-j and order |l-k|
                        const FComplexe* p_Inner_term = L2L_tranfer + expansion_Redirection_array_for_j[n-j] + l-k;

                        //printf("1\n");
                        for ( ; l-k>0;  --l, --p_src_exp_term, --p_Inner_term){ /* l>0 && l-k>0 */
                                p_target_exp_term->incReal( (p_src_exp_term->getReal() * p_Inner_term->getReal()) -
                                                            (p_src_exp_term->getImag() * p_Inner_term->getImag()));
                                p_target_exp_term->incImag( (p_src_exp_term->getImag() * p_Inner_term->getReal()) +
                                                            (p_src_exp_term->getReal() * p_Inner_term->getImag()));
                                /*printf("\t p_src_exp_term->real = %lf \t p_src_exp_term->imag = %lf \n",
                                       p_src_exp_term->getReal(),p_src_exp_term->getImag());
                                printf("\t p_Inner_term->real = %lf \t p_Inner_term->imag = %lf \n",
                                       p_Inner_term->getReal(),p_Inner_term->getImag());
                                printf("\t\t p_target_exp_term->real = %lf \t p_target_exp_term->imag = %lf \n",
                                       p_target_exp_term->getReal(),p_target_exp_term->getImag());
                                printf("\tp_target_exp_term = %d\n",p_target_exp_term-local_exp_target);*/
                        }

                        //printf("2\n");
                        // (-1)^l
                        double pow_of_minus_1_for_l = ((l%2) ? -1 : 1);
                        for (; l>0 && l>=j-n+k; --l, pow_of_minus_1_for_l = -pow_of_minus_1_for_l, --p_src_exp_term, ++p_Inner_term){ /* l>0 && l-k<=0 */
                                p_target_exp_term->incReal( pow_of_minus_1_for_l * pow_of_minus_1_for_k *
                                                ((p_src_exp_term->getReal() * p_Inner_term->getReal()) +
                                                (p_src_exp_term->getImag() * p_Inner_term->getImag())));
                                p_target_exp_term->incImag( pow_of_minus_1_for_l * pow_of_minus_1_for_k *
                                                ((p_src_exp_term->getImag() * p_Inner_term->getReal()) -
                                                (p_src_exp_term->getReal() * p_Inner_term->getImag())));
                                /*printf("\t p_src_exp_term->real = %lf \t p_src_exp_term->imag = %lf \n",
                                       p_src_exp_term->getReal(),p_src_exp_term->getImag());
                                printf("\t p_Inner_term->real = %lf \t p_Inner_term->imag = %lf \n",
                                       p_Inner_term->getReal(),p_Inner_term->getImag());
                                printf("\t\t p_target_exp_term->real = %lf \t p_target_exp_term->imag = %lf \n",
                                       p_target_exp_term->getReal(),p_target_exp_term->getImag());
                                printf("\tp_target_exp_term = %d\n",p_target_exp_term-local_exp_target);*/
                        }

                        //printf("3\n");
                        // l<=0 && l-k<=0
                        for (; l>=j-n+k; --l, ++p_src_exp_term, ++p_Inner_term){
                                p_target_exp_term->incReal( pow_of_minus_1_for_k *
                                                ((p_src_exp_term->getReal() * p_Inner_term->getReal()) -
                                                (p_src_exp_term->getImag() * p_Inner_term->getImag())));
                                p_target_exp_term->decImag( pow_of_minus_1_for_k *
                                                ((p_src_exp_term->getImag() * p_Inner_term->getReal()) +
                                                (p_src_exp_term->getReal() * p_Inner_term->getImag())));
                                /*printf("\t p_src_exp_term->real = %lf \t p_src_exp_term->imag = %lf \n",
                                       p_src_exp_term->getReal(),p_src_exp_term->getImag());
                                printf("\t p_Inner_term->real = %lf \t p_Inner_term->imag = %lf \n",
                                       p_Inner_term->getReal(),p_Inner_term->getImag());
                                printf("\t\t p_target_exp_term->real = %lf \t p_target_exp_term->imag = %lf \n",
                                       p_target_exp_term->getReal(),p_target_exp_term->getImag());
                                printf("\tj=%d\tk=%d\tn=%d\tl=%d\tpow_of_minus_1_for_k=%f\n",j,k,n,l,pow_of_minus_1_for_k);
                                printf("\tp_target_exp_term = %d\n",p_target_exp_term-local_exp_target);*/
                        }
                        //printf("\tj=%d\tk=%d\tn=%d\tl=%d\n",j,k,n,l);
                        //printf("\t\t p_target_exp_term->real = %lf \t p_target_exp_term->imag = %lf \n",
                        //       p_target_exp_term->getReal(),p_target_exp_term->getImag());
                    }
                }
            }
        }
    }

    virtual void L2P(const CellClass* const local, FList<ParticuleClass*>* const particules) = 0;

    virtual void P2P(FList<ParticuleClass*>* const currentBox, const FList<ParticuleClass*>*const* directNeighbors, const int size) = 0;
};


template< class ParticuleClass, class CellClass>
const double FAbstractFmbKernels<ParticuleClass,CellClass>::PiArrayInner[4] = {0, FMath::FPiDiv2, FMath::FPi, -FMath::FPiDiv2};


template< class ParticuleClass, class CellClass>
const double FAbstractFmbKernels<ParticuleClass,CellClass>::PiArrayOuter[4] = {0, -FMath::FPiDiv2, FMath::FPi, FMath::FPiDiv2};



#endif //FABSTRACTFMBKERNELS_HPP

// [--LICENSE--]
