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
#ifndef FEXTENDFORCES_HPP
#define FEXTENDFORCES_HPP


#include "../Utils/FGlobal.hpp"
#include "../Utils/FPoint.hpp"
#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendForces
* Please read the license
*
* This class is an extenssion.
* It proposes a 3d array as a forces vector.
*/
class FExtendForces {
protected:
    FPoint forces; //< 3D vector stored in a position object

public:
    /** Default constructor */
    FExtendForces() {
    }

    /** Copy constructor */
    FExtendForces(const FExtendForces& other) : forces(other.forces) {
    }

    /** Copy operator */
    FExtendForces& operator=(const FExtendForces& other) {
        this->forces = other.forces;
        return *this;
    }

    /** Return the forces */
    const FPoint& getForces() const {
        return this->forces;
    }

    /** Set Forces */
    void incForces(const FPoint& inForces) {
        this->forces += inForces;
    }

    /** Set Forces with 3 FReals */
    void incForces(const FReal inFx, const FReal inFy, const FReal inFz) {
        this->forces.incX(inFx);
        this->forces.incY(inFy);
        this->forces.incZ(inFz);
    }

    /** set the forces from 3 variables */
    void setForces(const FReal inFx, const FReal inFy, const FReal inFz) {
        this->forces.setPosition(inFx , inFy, inFz);
    }

    /** Save current object */
    void save(FBufferWriter& buffer) const {
        forces.save(buffer);
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        forces.restore(buffer);
    }
};


#endif //FEXTENDFORCES_HPP


