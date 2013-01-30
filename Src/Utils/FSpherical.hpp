// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FSPHERICAL_HPP
#define FSPHERICAL_HPP

#include "FGlobal.hpp"
#include "FMath.hpp"
#include "FPoint.hpp"
#include "FDebug.hpp"

/**
* This class is a Spherical position
*
* @brief Spherical coordinate system
*
* The spherical coordinate system (r radius, theta inclination, phi azimuth) is the following <br>
* x = r sin(theta) cos(phi) <br>
* y = r sin(theta) sin(phi) <br>
* z = r cos(theta) <br>
* This system is defined in p 872 of the paper of Epton and Dembart, SIAM J Sci Comput 1995.<br>
*
* Even if it can look different from usual expression (where theta and phi are inversed),
* such expression is used to match the SH expression.
*/
class FSpherical {
    // The attributes of a sphere
    FReal r;         //!< the radius
    FReal theta;     //!< the inclination angle [0, pi] - colatitude
    FReal phi;       //!< the azimuth angle [-pi,pi] - longitude - around z axis
    FReal cosTheta;
    FReal sinTheta;
public:
    /** Default Constructor, set attributes to 0 */
    FSpherical()
        : r(0), theta(0), phi(0), cosTheta(0), sinTheta(0) {
    }

    /** From now, we just need a constructor based on a 3D position */
    explicit FSpherical(const FPoint& inVector){
        const FReal x2y2 = (inVector.getX() * inVector.getX()) + (inVector.getY() * inVector.getY());
        this->r          = FMath::Sqrt( x2y2 + (inVector.getZ() * inVector.getZ()));

        this->phi        = FMath::Atan2(inVector.getY(),inVector.getX());

        this->cosTheta = inVector.getZ() / r;
        this->sinTheta = FMath::Sqrt(x2y2) / r;
        this->theta    = FMath::ACos(this->cosTheta);
        // if r == 0 we cannot divide!
        FDEBUG(if( r < FMath::Epsilon ) FDebug::Controller << "!!! In FSpherical, r == 0!\n"; )
    }

    /** Get the radius */
    FReal getR() const{
        return r;
    }

    /** Get the inclination angle theta = acos(z/r) [0, pi] */
    FReal getTheta() const{
        return theta;
    }
    /** Get the azimuth angle phi = atan2(y,x) [-pi,pi] */
    FReal getPhi() const{
        return phi;
    }

    /** Get the inclination angle [0, pi] */
    FReal getInclination() const{
        return theta;
    }
    /** Get the azimuth angle [0,2pi] */
    FReal getAzimuth() const{
        return (phi < 0 ? FMath::FPi*2 + phi : phi);
    }

    /** Get the cos of theta = z / r */
    FReal getCosTheta() const{
        return cosTheta;
    }

    /** Get the sin of theta = sqrt(x2y2) / r */
    FReal getSinTheta() const{
        return sinTheta;
    }
};

#endif // FSPHERICAL_HPP
