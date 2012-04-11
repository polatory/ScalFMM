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
#ifndef FSPHERICAL_HPP
#define FSPHERICAL_HPP

#include "FGlobal.hpp"
#include "FMath.hpp"
#include "FPoint.hpp"

/** This class is a Spherical position
* This class is currently in its minimum version
*/
//! The spherical coordinate system is the following
//! x = r sin(theta) cos(phi)
//! y = r sin(theta) sin(phi)
//! z = r cos(theta)
//! This system is defined in p 872 of the paper of Epton and Dembart, SIAM J Sci Comput 1995.
//!
//! \bug Error of naming : phi is in fact theta and theta is phi
//!  z = r cos Phi --> CosPhi = Z/r
//!  getPhi returns in fact the theta angle !!!!
//! To change
class FSpherical {
    // The attributes of a sphere
    FReal r;
    FReal cosTheta;
    FReal sinTheta;
    FReal phi;
//    FReal theta;

public:
    /** From now, we just need a constructor based on a 3D position */
    explicit FSpherical(const FPoint& inVector){
        //
//        // The good code for me but expansions are wrong
//        this->r = FMath::Sqrt( inVector.getX() * inVector.getX() + inVector.getY() * inVector.getY()
//                               + inVector.getZ() * inVector.getZ());
//        this->theta = FMath::Atan2(inVector.getY(),inVector.getX());
//        cosTheta = FMath::Cos(this->theta);
//        sinTheta = FMath::Sin(this->theta);
//        if( r < FMath::Epsilon ){ // if r == 0 we cannot divide!
//            this->phi = 0;
//        }
//        else {
//            this->phi = FMath::ACos(inVector.getZ() / r );
//        }
        const FReal x2y2 = (inVector.getX() * inVector.getX()) + (inVector.getY() * inVector.getY());
        this->r = FMath::Sqrt( x2y2 + (inVector.getZ() * inVector.getZ()));
        this->phi = FMath::Atan2(inVector.getY(),inVector.getX());
        if( r < FMath::Epsilon ){ // if r == 0 we cannot divide!
            this->cosTheta = FReal(1);
            this->sinTheta = FReal(1);
        }
        else {
            this->cosTheta = inVector.getZ() / r;
            this->sinTheta = FMath::Sqrt(x2y2) / r;
        }
    }

    /** Get the rayon */
    FReal getR() const{
        return r;
    }

    /** Get the cos theta = z / r */
    FReal getCosTheta() const{
        return cosTheta;
    }

    /** Get the sin theta = sqrt(x2y2) / r */
    FReal getSinTheta() const{
        return sinTheta;
    }

    /** Get the phi = atan2(y,x) */
    FReal getPhi() const{
        return phi;
    }
};

#endif // FSPHERICAL_HPP
