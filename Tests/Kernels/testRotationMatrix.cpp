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

#include <limits>
#include <iostream>

#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"
#include "../../Src/Kernels/Spherical/FSphericalRotationKernel.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"


int atLm(const int l, const int m){
        return (l*(l+1))/2 + m;
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


FReal analytical(const FReal cosTheta, const FReal sinTheta, const int l , const int m, const int k){
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

FReal analyticalDachsel(const FReal cosTheta, const FReal sinTheta, const int l , const int m, const int k){
    /*if(l >= 0 && -l <= m && m <= l && FMath::Abs(m) <= k && k <= l) {
        return FMath::pow(-1.0, k+m) * analyticalDachsel(cosTheta,sinTheta, l, k ,m);
    }
    else if( l > 0 && -l <= m && m <= l-1 && -l <= k && k <= -(m+1) ){
        return FMath::pow(-1.0, m+k) * analyticalDachsel(cosTheta,sinTheta, l, -m, -k);
    }*/
    if( m < k ){
        return FMath::pow(-1.0, k+m) * analyticalDachsel(cosTheta,sinTheta, l, k ,m);
    }
    else if( m < -k ){
        return FMath::pow(-1.0, m+k) * analyticalDachsel(cosTheta,sinTheta, l, -m, -k);
    }
    else {
        FReal sum = 0;
        for(int n = 0 ; n <= l-m ; ++n){
            sum += FMath::pow(-1.0, l-m-n) * combin(l-k, n) * combin(l+k, l-m-n) * FMath::pow(1.0+cosTheta,n) * FMath::pow(1.0-cosTheta, l-m-n);
        }
        return (1.0/FMath::pow(2.0,l)) * FMath::Sqrt(FReal(fact(l-m)*fact(l+m))/FReal(fact(l-k)*fact(l+k)))
                * FMath::pow(1.0 + Sgn(k)*cosTheta, FMath::Abs(k)) * FMath::pow(sinTheta, m - FMath::Abs(k)) * sum;
    }
}



int main(){
    static const int P = 2;
    static const int SizeArray = (P+1)*(P+1);

    //const FReal cosTheta = FMath::Cos(FMath::FPiDiv2/FReal(2.0));
    //const FReal sinTheta = FMath::Sin(FMath::FPiDiv2/FReal(2.0));
    const FPoint relativPos(FReal(0.1),FReal(0.1),FReal(0.1));
    const FSpherical relativPosSphere(relativPos);
    //Theta 0.955317 Phi 2.356194
    const FReal Betha = 0.955317;//FMath::FPi/2; //2*FMath::FPi-FMath::FPi/2
    const FReal cosTheta = FMath::Cos(Betha); // relativPosSphere.getCosTheta();
    const FReal sinTheta = FMath::Sin(Betha); // relativPosSphere.getSinTheta();

    std::cout << "cosTheta = " << cosTheta << "\n";
    std::cout << "sinTheta = " << sinTheta << "\n";

    // Rotation matrix d[l][m][k]
    FReal d[P+1][2*P+1][2*P+1];
    FMemUtils::setall((FReal*)d, FReal(std::numeric_limits<FReal>::quiet_NaN()), (P+1)*(2*P+1)*(2*P+1));

    /////////////////////////////////////////////////////////////
    // First method, Page 3,4
    // Fast and accurate determination of Wigner rotation matrices
    // in FMM
    /////////////////////////////////////////////////////////////

    FReal legendre_f[SizeArray];
    { // init legendre-f values
        // Equ 24
        // P~{0,0} = 1
        legendre_f[0] = 1.0;

        { // P~{l,l} = sqrt(2l-1 / 2l) sin(theta) P~{l-1,l-1} , For l > 0
            for(int l = 1; l <= P ; ++l ){
                legendre_f[atLm(l,l)] = FMath::Sqrt(FReal((2.0*l-1) / (2.0*l))) * sinTheta * legendre_f[atLm(l-1,l-1)];
            }
        }
        { // P~{l,l-1} = sqrt(2l-1) cos(theta) P~{l-1,l-1} , For l > 0
            for(int l = 1; l <= P ; ++l ){
                legendre_f[atLm(l,l-1)] = FMath::Sqrt(FReal(2*l-1)) * cosTheta * legendre_f[atLm(l-1,l-1)];
            }
        }
        {
            // For l > 1, 0 <= k < l-1
            // P~{l,k} = (2l-1) cos(theta) P{l-1,k} - sqrt((l-k-1)(l+k-1) P~{l-2,k}
            //          / sqrt((l-k)(l+k))
            for(int l = 2; l <= P ; ++l ){
                for( int k = 0 ; k < l - 1 ; ++k ){
                    legendre_f[atLm(l,k)] =
                            (FReal(2*l - 1 ) * cosTheta * legendre_f[atLm(l-1,k)] - FMath::Sqrt(FReal((l-k-1)*(l+k-1))) * legendre_f[atLm(l-2,k)])
                                                          / FMath::Sqrt(FReal((l-k)*(l+k)));
                }
            }
        }
    }

    { // Initial condition
        // Equ 21
        // For l > 0,
        // & -l <= k < 0, d{l,0,k} = sqrt((l+k)!/(l-k)!) P{l,-k} = P~{l,-k}
        // & k = 0, d{l,0,k} = P{l,0}
        // & 0 < k <= l, d{l,0,k} = -1^k sqrt((l-k)!/(l+k)!) P{l,k} = -1^k P~{l,k}
        for(int l = 0 ; l <= P ; ++l){
            // k < 0
            for(int k = -l ; k < 0 ; ++k){
                d[l][P+0][P+k] = legendre_f[atLm(l,-k)];
            }
            // k == 0
            d[l][P+0][P+0] = legendre_f[atLm(l,0)];
            // 0 < k
            for(int k = 1 ; k <= l ; ++k){
                d[l][P+0][P+k] = FMath::pow(FReal(-1),k) * legendre_f[atLm(l,k)];
            }
        }
    }
    {
        for(int l = 1 ; l <= P ; ++l){
            for(int m = 0 ; m < l ; ++m){
                // Equ 18
                // For l > 0, 0 <= m < l, -l < k <= l, cos(theta) >= 0
                // d{l,m+1,k} = sqrt( l(l+1) - k(k-1) / l(l+1) - m(m+1)) d{l,m,k-1}
                //            - (m+k) sin(theta) d{l,m,k} / sqrt(l(l+1) - m(m+1)) (1+cos(theta))
                for(int k = -l+1 ; k <= l ; ++k){
                    d[l][P+m+1][P+k] =
                            FMath::Sqrt(FReal(l*(l+1)-k*(k-1))/FReal(l*(l+1)-m*(m+1))) * d[l][P+m][P+k-1]
                            - FReal(m+k)*sinTheta*d[l][P+m][P+k]/(FMath::Sqrt(FReal(l*(l+1)-m*(m+1)))*(1+cosTheta));
                }
                // Equ 19
                // For l > 0, 0 <= m < l, cos(theta) >= 0
                // d{l,m+1,-l} = (l-m) sin(theta) d{l,m,-l}
                //             / sqrt(l(l+1)-m(m+1)) (1+cos(theta))
                d[l][P+m+1][P-l] = FReal(l-m)*sinTheta*d[l][P+m][P-l]/(FMath::Sqrt(FReal(l*(l+1)-m*(m+1)))*(1+cosTheta));
            }
            // Equ 20
            // d{l,m,k} = -1^(m+k) d{l,-m,-k}  , For l > 0, -l <= m < 0, -l <= k <= l
            for(int m = -l ; m < 0 ; ++m){
                for(int k = -l ; k <= l ; ++k){
                    d[l][P+m][P+k] = FMath::pow( FReal(-1.0), m+k) * d[l][P-m][P-k];
                }
            }
        }
    }

    // Print
    std::cout << "Print result\n";
    for(int l = 0 ; l <= P ; ++l){
        for(int m = -l ; m <= l ; ++m ){
            for(int k = -l ; k <= l ; ++k ){
                std::cout << d[l][P+m][P+k] << "\t";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    FMemUtils::setall((FReal*)d, FReal(std::numeric_limits<FReal>::quiet_NaN()), (P+1)*(2*P+1)*(2*P+1));

    /////////////////////////////////////////////////////////////
    // Second method, Page 4
    // Fast and accurate determination of Wigner rotation matrices
    // in FMM
    /////////////////////////////////////////////////////////////

    FReal g[SizeArray];
    {// Equ 29
        // g{0,0} = 1
        g[0] = 1;

        // g{l,0} = sqrt( (2l - 1) / 2l) g{l-1,0}  for l > 0
        for(int l = 1; l <= P ; ++l ){
            g[atLm(l,0)] = FMath::Sqrt(FReal((l*2.0-1.0)/(l*2.0))) * g[atLm(l-1,0)];
        }

        // g{l,m} = sqrt( (l - m + 1) / (l+m)) g{l,m-1}  for l > 0, 0 < m <= l
        for(int l = 1; l <= P ; ++l ){
            for(int m = 1; m <= l ; ++m ){
                g[atLm(l,m)] = FMath::Sqrt(FReal((l-m+1))/FReal((l+m))) * g[atLm(l,m-1)];
            }
        }
    }
    { // initial
        // Equ 28
        // d{l,m,l} = -1^(l+m) g{l,m} (1+cos(theta))^m sin(theta)^(l-m) , For l > 0, 0 <= m <= l
        for(int l = 0 ; l <= P ; ++l){
            for(int m = 0 ; m <= l ; ++m){
                d[l][P+m][P+l] = FMath::pow( FReal(-1.0), l+m) * g[atLm(l,m)] * FMath::pow( FReal(1) + cosTheta, m) * FMath::pow(sinTheta, l-m);
            }
        }
    }
    {
        for(int l = 1 ; l <= P ; ++l){
            for(int k = l ; k > -l ; --k){
                // Equ 25
                // For l > 0, 0 <= m < l, -l < k <= l, cos(theta) >= 0
                // d{l,m,k-1} = sqrt( l(l+1) - m(m+1) / l(l+1) - k(k-1)) d{l,m+1,k}
                //            + (m+k) sin(theta) d{l,m,k} / sqrt(l(l+1) - k(k-1)) (1+cos(theta))
                for(int m = 0 ; m < l ; ++m){
                    d[l][P+m][P+k-1] =
                            (FMath::Sqrt(FReal(l*(l+1)-m*(m+1))/FReal(l*(l+1)-k*(k-1))) * d[l][P+m+1][P+k])
                            + (FReal(m+k)*sinTheta*d[l][P+m][P+k]/(FMath::Sqrt(FReal(l*(l+1)-k*(k-1)))*(1+cosTheta)));
                }
                // Equ 26
                // For l > 0, -l < k <= l, cos(theta) >= 0
                // d{l,l,k-1} = (l+k) sin(theta) d{l,l,k}
                //             / sqrt(l(l+1)-k(k-1)) (1+cos(theta))
                d[l][P+l][P+k-1] = FReal(l+k)*sinTheta*d[l][P+l][P+k]/(FMath::Sqrt(FReal(l*(l+1)-k*(k-1)))*(1+cosTheta));
            }
            // Equ 27
            // d{l,m,k} = -1^(m+k) d{l,-m,-k}  , For l > 0, -l <= m < 0, -l <= k <= l
            for(int m = -l ; m < 0 ; ++m){
                for(int k = -l ; k <= l ; ++k){
                    d[l][P+m][P+k] = FMath::pow( FReal(-1), m+k) * d[l][P-m][P-k];
                }
            }
        }
    }

    // Print
    std::cout << "Print result\n";
    for(int l = 0 ; l <= P ; ++l){
        for(int m = -l ; m <= l ; ++m ){
            for(int k = -l ; k <= l ; ++k ){
                std::cout << d[l][P+m][P+k] << "\t";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    FMemUtils::setall((FReal*)d, FReal(std::numeric_limits<FReal>::quiet_NaN()), (P+1)*(2*P+1)*(2*P+1));


    /////////////////////////////////////////////////////////////
    // Analatycall values
    /////////////////////////////////////////////////////////////
    std::cout << "Print result analytical\n";
    for(int l = 0 ; l <= P ; ++l){
        for(int m = -l ; m <= l ; ++m ){
            for(int k = -l ; k <= l ; ++k ){
                std::cout << analytical(cosTheta,sinTheta,l,m,k) << "\t";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    /////////////////////////////////////////////////////////////
    // Analatycall values V2
    /////////////////////////////////////////////////////////////
    std::cout << "Print result analytical V2\n";
    for(int l = 0 ; l <= P ; ++l){
        for(int m = -l ; m <= l ; ++m ){
            for(int k = -l ; k <= l ; ++k ){
                std::cout << analyticalDachsel(cosTheta,sinTheta,l,m,k) << "\t";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    /////////////////////////////////////////////////////////////
    // correct result for l=1
    /////////////////////////////////////////////////////////////
    std::cout << "Print L=1\n";
    {
        FReal resL[9];
        resL[0] = FMath::pow(FMath::Cos(Betha/2.0),2);
        resL[1] = -FMath::Sin(Betha)/FMath::Sqrt(2.0);
        resL[2] = FMath::pow(FMath::Sin(Betha/2.0),2);
        resL[3] = FMath::Sin(Betha)/FMath::Sqrt(2.0);;
        resL[4] = FMath::Cos(Betha);
        resL[5] = -FMath::Sin(Betha)/FMath::Sqrt(2.0);;
        resL[6] = FMath::pow(FMath::Sin(Betha/2.0),2);
        resL[7] = FMath::Sin(Betha)/FMath::Sqrt(2.0);
        resL[8] = FMath::pow(FMath::Cos(Betha/2.0),2);

        for(int m = 0 ; m < 3 ; ++m ){
            for(int k = 0 ; k < 3 ; ++k ){
                std::cout << resL[m*3 + k] << "\t";
            }
            std::cout << "\n";
        }
    }
    /////////////////////////////////////////////////////////////
    // correct result for l=2
    /////////////////////////////////////////////////////////////
    std::cout << "Print L=2\n";
    {
        FReal resL[25];
        resL[0] = FMath::pow(FMath::Cos(Betha/2.0),4);
        resL[1] = (-FMath::pow(FMath::Cos(Betha/2.0),2))*FMath::Sin(Betha);
        resL[2] = (1.0/2.0)*FMath::Sqrt(3.0/2.0)*FMath::pow(FMath::Sin(Betha),2);
        resL[3] = (-FMath::pow(FMath::Sin(Betha/2.0),2))*FMath::Sin(Betha);
        resL[4] = FMath::pow(FMath::Sin(Betha/2.0),4);

        resL[5] = FMath::pow(FMath::Cos(Betha/2.0),2)*FMath::Sin(Betha);
        resL[6] = FMath::pow(FMath::Cos(Betha),2)- FMath::pow(FMath::Sin(Betha/2),2);
        resL[7] = (-1.0/2.0)*FMath::Sqrt(3.0/2.0)*FMath::Sin(2*Betha);
        resL[8] = FMath::pow(FMath::Cos(Betha/2),2)- FMath::pow(FMath::Cos(Betha),2);
        resL[9] = (-FMath::pow(FMath::Sin(Betha/2.0),2))*FMath::Sin(Betha);

        resL[10] = (1.0/2.0)*FMath::Sqrt(3.0/2.0)*FMath::pow(FMath::Sin(Betha),2);
        resL[11] = (1.0/2.0)*FMath::Sqrt(3.0/2.0)*FMath::Sin(Betha*2);
        resL[12] = (1.0/2.0)*(-1+3*FMath::pow(FMath::Cos(Betha),2));
        resL[13] = (-1.0/2.0)*FMath::Sqrt(3.0/2.0)*FMath::Sin(Betha*2);
        resL[14] = (1.0/2.0)*FMath::Sqrt(3.0/2.0)*FMath::pow(FMath::Sin(Betha),2);

        resL[15] = FMath::pow(FMath::Sin(Betha/2.0),2)*FMath::Sin(Betha);
        resL[16] = FMath::pow(FMath::Cos(Betha/2),2)- FMath::pow(FMath::Cos(Betha),2);
        resL[17] = (1.0/2.0)*FMath::Sqrt(3.0/2.0)*FMath::Sin(2*Betha);
        resL[18] = FMath::pow(FMath::Cos(Betha),2)- FMath::pow(FMath::Sin(Betha/2),2);
        resL[19] = (-FMath::pow(FMath::Cos(Betha/2.0),2))*FMath::Sin(Betha);

        resL[20] = FMath::pow(FMath::Sin(Betha/2.0),4);
        resL[21] = FMath::pow(FMath::Sin(Betha/2.0),2)*FMath::Sin(Betha);
        resL[22] = (1.0/2.0)*FMath::Sqrt(3.0/2.0)*FMath::pow(FMath::Sin(Betha),2);
        resL[23] = FMath::pow(FMath::Cos(Betha/2.0),2)*FMath::Sin(Betha);
        resL[24] = FMath::pow(FMath::Cos(Betha/2.0),4);


        for(int m = 0 ; m < 5 ; ++m ){
            for(int k = 0 ; k < 5 ; ++k ){
                std::cout << resL[m*5 + k] << "\t";
            }
            std::cout << "\n";
        }
    }
    /////////////////////////////////////////////////////////////
    // Third method, Page 5 (+105)
    // Rotating around the quartic angular momentum barrier in fast multipole
    /////////////////////////////////////////////////////////////

//    FReal legendre[SizeArray];
//    {
//        const FReal x = cosTheta;
//        const FReal minus_sqrt_1_minus_x_pow_2 = -sinTheta;

//        legendre[0] = 1.0;        // P_0,0(cosTheta) = 1
//        legendre[1] = x;          // P_1,0(cosTheta) = cosTheta
//        legendre[2] = minus_sqrt_1_minus_x_pow_2;// P_1,1(cosTheta) = -sinTheta = -sqrt(1-x^2)

//        FReal fact = 3.0;

//        for(int l = 2; l <= P ; ++l ){
//            const FReal x_l2_minus_1 = x * FReal( 2 * l - 1 );
//            // m from 0 to l - 2
//            for( int m = 0; m <= l - 2 ; ++m ){
//                //         cosTheta x (2 x l - 1) x P_l-1_m - (l + m - 1) x P_l-2_m
//                // P_l_m = --------------------------------------------------------
//                //                               l - m
//                legendre[atLm(l,m)] = (x_l2_minus_1 * legendre[atLm(l-1,m)] - FReal( l + m - 1 ) * legendre[atLm(l-2,m)] )
//                                                                    / FReal( l - m );
//            }
//            // Compute P_l,{l-1}
//            legendre[atLm(l,l-1)] = x_l2_minus_1 * legendre[atLm(l,l-1)];

//            // Compute P_l,l
//            legendre[atLm(l,l)] = fact * minus_sqrt_1_minus_x_pow_2 * legendre[atLm(l,l-1)];

//            fact += FReal(2.0);
//        }
//    }
//    { // initial
//        // Equ A8
//        // d{l,0,m} = -1^(m) P{l,m} , For l >= 0, -l <= m <= l
//        for(int l = 0 ; l <= P ; ++l){
//            for(int m = -l ; m <= l ; ++m){
//                d[l][P+0][P+m] = FMath::pow( FReal(-1.0), m) * legendre[atLm(l,m)];
//            }
//        }
//    }
//    {
//        // Equ A4
//        // For l > 0, l < m <= l, -l <= k < l, cos(theta) >= 0
//        // d{l,k+1,m} = -(m+k) sin(theta) d{l,k,m} / sqrt(l(l+1) - k(k+1)) (1+cos(theta))
//        //              + sqrt( l(l+1) - m(m+1) / l(l+1) - k(k+1)) d{l,k,m-1}
//        //
//        for(int l = 1 ; l <= P ; ++l){
//            for(int k = 0 ; k < l ; ++k){
//                for(int m = -l+1 ; m <= l ; ++m){
//                    d[l][P+k+1][P+m] =
//                            -FReal(m+k)*sinTheta*d[l][P+k][P+m]/(FMath::Sqrt(FReal(l*(l+1)-k*(k+1)))*(1+cosTheta))
//                            + FMath::Sqrt(FReal(l*(l+1)-m*(m+1))/FReal(l*(l+1)-k*(k+1))) * d[l][P+k][P+m-1];
//                }
//                d[l][P+k+1][P-l] = -FReal(k-l)*sinTheta*d[l][P+k][P-l]/(FMath::Sqrt(FReal(l*(l+1)-k*(k+1)))*(1+cosTheta));
//            }
//            // Sym
//            for(int k = -1 ; k >= -l ; --k){
//                for(int m = -l ; m <= l ; ++m){
//                    d[l][P+k][P+m] = FMath::pow( FReal(-1), m+k) * d[l][P-k][P-m];
//                }
//            }
//        }
//    }

//    // Print
//    std::cout << "Print result\n";
//    for(int l = 0 ; l <= P ; ++l){
//        for(int m = -l ; m <= l ; ++m ){
//            for(int k = -l ; k <= l ; ++k ){
//                std::cout << d[l][P+m][P+k] << "\t";
//            }
//            std::cout << "\n";
//        }
//        std::cout << "\n";
//    }
//    FMemUtils::setall((FReal*)d, FReal(std::numeric_limits<FReal>::quiet_NaN()), (P+1)*(2*P+1)*(2*P+1));

    return 0;
}

