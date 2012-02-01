#ifndef FEXTENDVELOCITY_HPP
#define FEXTENDVELOCITY_HPP

// [--License--]

#include "../Utils/FGlobal.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendVelocity
* Please read the license
*
* This class is an extenssion.
* It proposes a 3d array as a velocity vector.
*/
class FExtendVelocity {
protected:
    F3DPosition velocity; //< 3D vector stored in a position object

public:
    /** Default constructor */
    FExtendVelocity() {
    }

    /** Copy constructor */
    FExtendVelocity(const FExtendVelocity& other) : velocity(other.velocity) {
    }

    /** Copy operator */
    FExtendVelocity& operator=(const FExtendVelocity& other) {
        this->velocity = other.velocity;
        return *this;
    }

    /** Return the velocity */
    const F3DPosition& getVelocity() const {
        return this->velocity;
    }

    /** Set Velocity */
    void incVelocity(const F3DPosition& inVelocity) {
        this->velocity += inVelocity;
    }

    /** Set Velocity with 3 FReals */
    void incVelocity(const FReal inVx, const FReal inVy, const FReal inVz) {
        this->velocity.incX(inVx);
        this->velocity.incY(inVy);
        this->velocity.incZ(inVz);
    }

    /** set the velocity from 3 variables */
    void setVelocity(const FReal inVx, const FReal inVy, const FReal inVz) {
        this->velocity.setPosition(inVx , inVy, inVz);
    }
};


#endif //FEXTENDVELOCITY_HPP

// [--END--]
