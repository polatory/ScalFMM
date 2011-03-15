#ifndef FFMBKERNELS_HPP
#define FFMBKERNELS_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FAbstractKernels.hpp"

#include "../Containers/FTreeCoordinate.hpp"

#include "../Utils/F3DPosition.hpp"
#include "../Utils/FComplexe.hpp"

#include <iostream>


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class AbstractKernels
* @brief
* Please read the license
*
* This kernels simply shows the details of the information
* it receives
*/
template< class ParticuleClass, class CellClass>
class FFmbKernels : public FAbstractKernels<ParticuleClass,CellClass> {
    // P is a input parameter
    static const int FMB_Info_P = 2;

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
    static const int halphSize1Dim =  (2*(FMB_Info_ws)+1) +1;

    // EXPANSION_SIZE(FMB_Info.P)
    static const int FMB_Info_exp_size = ((FMB_Info_P)+1) * ((FMB_Info_P)+2) * 0.5;
    // NEXP_SIZE(FMB_Info.P)
    static const int FMB_Info_nexp_size = (FMB_Info_P + 1) * (FMB_Info_P + 1);

    // Level of the tree
    const int treeLevel;
    // Width of the box at the root level
    const double treeWidth;

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

    // [OK?] p_precomputed_cos_and_sin_array
    FComplexe cosSin[FMB_Info_M2L_P + 1];
    // [OK?] p_associated_Legendre_function_Array
    double legendre[FMB_Info_M2L_exp_size];

    // pow_of_I_array
    static const double PiArray[4];

    // To store spherical position
    struct Spherical {
        double r, cosTheta, sinTheta, phi;
    };

    //[Blas] int expansion_Redirection_array_for_j[FMB_Info_M2L_P + 1 ];


    //////////////////////////////////////////////////////////////////
    // Allocation
    //////////////////////////////////////////////////////////////////

    /*[Blas] void expansion_Redirection_array_for_j_Initialize() {
            for( int h = 0; h <= FMB_Info_M2L_P ; ++h ){
                expansion_Redirection_array_for_j[h] = static_cast<int>( h * ( h + 1 ) * 0.5 );
            }
    }*/

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
        this->transitionM2M = new FComplexe**[this->treeLevel];
        this->transitionL2L = new FComplexe**[this->treeLevel];

        for(int idxLevel = 0 ; idxLevel < this->treeLevel ; ++idxLevel){

            this->transitionM2M[idxLevel] = new FComplexe*[8];
            this->transitionL2L[idxLevel] = new FComplexe*[8];

            for(long idxChild = 0; idxChild < 8; ++idxChild){
                this->transitionM2M[idxLevel][idxChild] = new FComplexe[FMB_Info_exp_size];
                this->transitionL2L[idxLevel][idxChild] = new FComplexe[FMB_Info_exp_size];
            }
        }
        // M2L
        this->transferM2L = new FComplexe****[this->treeLevel+1];

        for(int idxLevel = 0; idxLevel <= this->treeLevel; ++idxLevel){
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
        for(long idxLevel = 0 ; idxLevel < this->treeLevel ; ++idxLevel){
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
        for(int idxLevel = 0 ; idxLevel <= this->treeLevel; ++idxLevel){
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
    void positionToSphere(const F3DPosition& vector, Spherical* const sphere ){
        const double x2y2 = (vector.getX() * vector.getX()) + (vector.getY() * vector.getY());

        sphere->r = FMath::Sqrt( x2y2 + (vector.getZ() * vector.getZ()));
        sphere->phi = FMath::Atan2(vector.getY(),vector.getX());
        sphere->cosTheta = vector.getZ() / sphere->r; // cos_th = z/r
        sphere->sinTheta = FMath::Sqrt(x2y2) / sphere->r; // sin_th = sqrt(x^2 + y^2)/r
    }

    // associated_Legendre_function_Fill_complete_array_of_values_for_cos
    void legendreFunction( const double cosTheta, const double sinTheta, double* const results ){
        // l=0:         results[current++] = 1.0; // P_0^0(cosTheta) = 1
        int idxCurrent = 0;
        results[idxCurrent++] = 1.0;

        // l=1:
        // Compute P_1^{0} using (3) and store it into results_array:
        results[idxCurrent++] = cosTheta;

        // Compute P_1^1 using (2 bis) and store it into results_array
        const double invSinTheta = -sinTheta; // somx2 => -sinTheta
        results[idxCurrent++] = invSinTheta;

        // l>1:
        int idxCurrent1m = 1; //pointer on P_{l-1}^m P_1^0
        int idxCurrent2m = 0; //pointer on P_{l-2}^m P_0^0
        double fact = 3.0;

        // Remark: p_results_array_l_minus_1_m and p_results_array_l_minus_2_m
        // just need to be incremented at each iteration.
        for(int idxl = 2; idxl <= this->LMax ; ++idxl ){
                for( int idxm = 0; idxm <= idxl - 2 ; ++idxm , ++idxCurrent , ++idxCurrent1m , ++idxCurrent2m ){
                        // Compute P_l^m, l >= m+2, using (1) and store it into results_array:
                        results[idxCurrent] = (cosTheta * ( 2 * idxl - 1 ) * results[idxCurrent1m] - ( idxl + idxm - 1 )
                                        * results[idxCurrent2m] ) / ( idxl - idxm );
                }
                // p_results_array_l_minus_1_m now points on P_{l-1}^{l-1}

                // Compute P_l^{l-1} using (3) and store it into ptrResults:
                results[idxCurrent++] = cosTheta * ( 2 * idxl - 1 ) * results[idxCurrent1m];

                // Compute P_l^l using (2 bis) and store it into results_array:
                results[idxCurrent++] = fact * invSinTheta * results[idxCurrent1m];

                fact += 2.0;
                ++idxCurrent1m;
        }
    }

    // spherical_harmonic_Inner
    void harmonicInner(const Spherical& sphere, FComplexe* const results){

        FComplexe* ptrCosSin = this->cosSin;
        for(int idxl = 0 , idxlMod4 = 0; idxl <= LMax ; ++idxl, ++idxlMod4, ++ptrCosSin){
            if(idxlMod4 == 4) idxlMod4 = 0;
            const double angle = idxl * sphere.phi + this->PiArray[idxlMod4];

            ptrCosSin->setReal( FMath::Sin(angle + FMath::FPiDiv2) );
            ptrCosSin->setImag( FMath::Sin(angle) );
            //std::cout<< idxl << "=" << ptrCosSin->getReal() << "/" << ptrCosSin->getImag() << " (" << idxl << "/" << sphere.phi << "/" << PiArray[lMod4] << ")\n";
        }

        legendreFunction(sphere.cosTheta, sphere.sinTheta, this->legendre);
        //printf("FMB_Info_M2L_exp_size=%d\n",FMB_Info_M2L_exp_size);
        //for(int temp = 0 ; temp < FMB_Info_M2L_exp_size ; ++temp){
        //    printf("%f\n",this->legendre[temp]);
        //}

        FComplexe* currentResult = results;
        int idxLegendre = 0;//ptr_associated_Legendre_function_Array
        int idxSphereHarmoCoef = 0;
        double idxRl = 1.0 ;

        for(int idxl = 0; idxl <= this->LMax ; ++idxl, idxRl *= sphere.r){
            for(int idxm = 0 ; idxm <= idxl ; ++idxm, ++currentResult, ++idxSphereHarmoCoef, ++idxLegendre){
                const double magnitude = this->sphereHarmoInnerCoef[idxSphereHarmoCoef] * idxRl * legendre[idxLegendre];
                currentResult->setReal( magnitude * this->cosSin[idxm].getReal() );
                currentResult->setImag( magnitude * this->cosSin[idxm].getImag() );
                //printf("magnitude=%f idxRl=%f sphereHarmoInnerCoef=%f real=%f imag=%f\n",magnitude,idxRl,this->sphereHarmoInnerCoef[idxSphereHarmoCoef],currentResult->getReal(),currentResult->getImag());
            }
        }

    }
    // spherical_harmonic_Outer
    void harmonicOuter(const Spherical& sphere, FComplexe* const results){

        FComplexe* ptrCosSin = this->cosSin;
        for(int idxl = 0, idxlMod4 = 0; idxl < LMax ; ++idxl, ++idxlMod4, ++ptrCosSin){
            if(idxlMod4 == 4) idxlMod4 = 0;
            const double angle = idxl * sphere.phi + this->PiArray[idxlMod4];

            ptrCosSin->setReal( FMath::Sin(angle + FMath::FPiDiv2) );
            ptrCosSin->setImag( FMath::Sin(angle) );
        }

        legendreFunction(sphere.cosTheta, sphere.sinTheta, this->legendre);

        int idxLegendre = 0;
        FComplexe* currentResult = results;

        const double invR = 1/sphere.r;
        double idxRl1 = invR;
        for(int idxl = 0 ; idxl <= this->LMax ; ++idxl, idxRl1 *= invR){
            for(int idxm = 0 ; idxm <= idxl ; ++idxm, ++currentResult, ++idxLegendre){
                const double magnitude = this->sphereHarmoOuterCoef[idxl-idxm] * idxRl1 * this->legendre[idxLegendre];
                currentResult->setReal( magnitude * this->cosSin[idxm].getReal() );
                currentResult->setImag( magnitude * this->cosSin[idxm].getImag() );
            }
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
        double treeWidthAtLevel = this->treeWidth/2;

        for(int idxLevel = 0 ; idxLevel < this->treeLevel ; ++idxLevel ){
            const F3DPosition father(treeWidthAtLevel,treeWidthAtLevel,treeWidthAtLevel);
            treeWidthAtLevel /= 2;

            //std::cout << "[precomputeM2M]treeWidthAtLevel=" << treeWidthAtLevel << "\n";

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

        double treeWidthAtLevel = this->treeWidth;
        for(int idxLevel = 0 ; idxLevel < this->treeLevel ; ++idxLevel ){

            for( int idxd1 = 0; idxd1 < this->size1Dim ; ++idxd1 ){

                for( int idxd2 = 0; idxd2 < this->size1Dim ; ++idxd2 ){

                    for( int idxd3 = 0; idxd3 < this->size1Dim ; ++idxd3 ){
                        const int x = idxd1 - this->halphSize1Dim;
                        const int y = idxd2 - this->halphSize1Dim;
                        const int z = idxd3 - this->halphSize1Dim;

                        if( ( x*x + y*y + z*z ) >= ( 3*FMB_Info_ws*FMB_Info_ws + 0.1 ) ){
                            const F3DPosition relativePos( x*treeWidthAtLevel , y*treeWidthAtLevel , z*treeWidthAtLevel );

                            Spherical spherical;
                            positionToSphere(relativePos,&spherical);
                            // Not blas so
                            harmonicOuter(spherical,this->transferM2L[idxLevel][idxd1][idxd2][idxd3]);
                            for(int idxTemp = 0 ; idxTemp < this->FMB_Info_M2L_exp_size ; ++idxTemp){
                                std::cout << "transferM2L[" << idxLevel << "][" << idxd1 << "][" << idxd2 << "][" << idxd3 << "][" << idxTemp << "]=" << this->transferM2L[idxLevel][idxd1][idxd2][idxd3][idxTemp].getReal()<<"/"<<this->transferM2L[idxLevel][idxd1][idxd2][idxd3][idxTemp].getImag()<< "\n";
                            }

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
        //[Blas] expansion_Redirection_array_for_j_Initialize();
        sphericalHarmonicInitialize();
        transferAllocate();

        precomputeM2M();
        precomputeM2L();
    }

public:
    FFmbKernels(const int inTreeLevel, const double inTreeWidth)
        : treeLevel(inTreeLevel), treeWidth(inTreeWidth),
          transitionM2M(0), transitionL2L(0) {
        buildPrecompute();
    }

    /** Default destructor */
    virtual ~FFmbKernels(){
        transferDeallocate();
    }

    /////////////////////////////////////////////////////////////////////////////////
    //    Init
    /////////////////////////////////////////////////////////////////////////////////

    void init(){
        // Nothing todo
    }

    /////////////////////////////////////////////////////////////////////////////////
    //    Upward
    /////////////////////////////////////////////////////////////////////////////////

    void P2M(CellClass* const pole, FList<ParticuleClass*>* const particules) {

    }

    void M2M(CellClass* const pole, CellClass** const child) {

    }

    void M2L(CellClass* const pole, CellClass** const distantNeighbors, const int size) {

    }

    /////////////////////////////////////////////////////////////////////////////////
    //    Downard
    /////////////////////////////////////////////////////////////////////////////////

    void L2L(CellClass* const pole, CellClass** const child) {

    }

    void L2P(CellClass* const pole, FList<ParticuleClass*>* const particules){

    }

    void P2P(FList<ParticuleClass*>* const currentBox, FList<ParticuleClass*>** directNeighbors, const int size) {

    }
};

template< class ParticuleClass, class CellClass>
const double FFmbKernels<ParticuleClass,CellClass>::PiArray[4] = {0, FMath::FPiDiv2, FMath::FPi, -FMath::FPiDiv2};



#endif //FFMBKERNELS_HPP

// [--LICENSE--]
