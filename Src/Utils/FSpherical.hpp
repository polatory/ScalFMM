#ifndef FSPHERICAL_HPP
#define FSPHERICAL_HPP

#include "FGlobal.hpp"
#include "F3DPosition.hpp"

/** This class is a Spherical position
* This class is currently in its minimum version
*/
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

#endif // FSPHERICAL_HPP
