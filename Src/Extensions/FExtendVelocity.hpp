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
#ifndef FEXTENDVELOCITY_HPP
#define FEXTENDVELOCITY_HPP



#include "../Utils/FGlobal.hpp"
#include "../Utils/FPoint.hpp"
#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendVelocity
* Please read the license
*
* This class is an extenssion.
* It proposes a FPoint as a velocity vector.
*/
class FExtendVelocity {
protected:
    FPoint velocity; //< 3D vector stored in a position object

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
    const FPoint& getVelocity() const {
        return this->velocity;
    }

    /** Set Velocity */
    void incVelocity(const FPoint& inVelocity) {
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

    /** Save current object */
    void save(FBufferWriter& buffer) const {
        buffer << velocity;
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        buffer >> velocity;
    }
};


#endif //FEXTENDVELOCITY_HPP


