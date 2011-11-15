#ifndef FEXTENDFORCES_HPP
#define FEXTENDFORCES_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FGlobal.hpp"
#include "../Utils/F3DPosition.hpp"

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
    F3DPosition forces; //< 3D vector stored in a position object

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
    const F3DPosition& getForces() const {
        return this->forces;
    }

    /** Set Forces */
    void incForces(const F3DPosition& inForces) {
        this->forces += inForces;
    }

    /** Set Forces with 3 FReals */
    void incForces(const FReal inFx, const FReal inFy, const FReal inFz) {
        this->forces.incX(inFx);
        this->forces.incY(inFy);
        this->forces.incZ(inFz);
    }

};


#endif //FEXTENDFORCES_HPP

// [--LICENSE--]
