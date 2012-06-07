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
#include "FDebug.hpp"

/** This class is a Spherical position
* This class is currently in its minimum version
*/

//! Brief Spherical coordinate system

//! The spherical coordinate system (r radius, theta inclination, phi azimuth) is the following <br>
//! x = r sin(theta) cos(phi) <br>
//! y = r sin(theta) sin(phi) <br>
//! z = r cos(theta) <br>
//! This system is defined in p 872 of the paper of Epton and Dembart, SIAM J Sci Comput 1995.
//!
class FSpherical {
    // The attributes of a sphere
    FReal r;         //!< the radius
    FReal cosTheta;
    FReal sinTheta;
    FReal theta;     //!< the inclination angle
    FReal phi;       //!< the azimuth angle
public:
    /** From now, we just need a constructor based on a 3D position */
    explicit FSpherical(const FPoint& inVector){
        const FReal x2y2 = (inVector.getX() * inVector.getX()) + (inVector.getY() * inVector.getY());
        this->r          = FMath::Sqrt( x2y2 + (inVector.getZ() * inVector.getZ()));
        this->phi        = FMath::Atan2(inVector.getY(),inVector.getX());
        if( r < FMath::Epsilon ){ // if r == 0 we cannot divide!
            FDEBUG( FDebug::Controller << "!!! In FSpherical, r == 0!\n"; )
            this->cosTheta = FReal(1);
            this->sinTheta = FReal(1);
        }
        else {
            this->cosTheta = inVector.getZ() / r;
            this->sinTheta = FMath::Sqrt(x2y2) / r;
        }
        this->theta    = FMath::ACos(this->cosTheta);
    }

    /** Get the radius */
    FReal getR() const{
        return r;
    }

    /** Get the cosine of theta = z / r */
    FReal getCosTheta() const{
        return cosTheta;
    }

    /** Get the sin theta = sqrt(x2y2) / r */
    FReal getSinTheta() const{
        return sinTheta;
    }
    /** Get the inclination angle theta = acos(z/r) */
    FReal getTheta() const{
        return theta;
    }
    /** Get the azimuth angle phi = atan2(y,x) */
    FReal getPhi() const{
        return phi;
    }
};

#endif // FSPHERICAL_HPP
