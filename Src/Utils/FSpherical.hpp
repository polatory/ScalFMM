#ifndef FSPHERICAL_HPP
#define FSPHERICAL_HPP

#include "FGlobal.hpp"
#include "FMath.hpp"
#include "F3DPosition.hpp"

/** This class is a Spherical position
* This class is currently in its minimum version
*/
class FSpherical {
    // The attributes of a sphere
    FReal r;
    FReal cosTheta;
    FReal sinTheta;
    FReal phi;

public:
    /** From now, we just need a constructor based on a 3D position */
    explicit FSpherical(const F3DPosition& inVector){
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
