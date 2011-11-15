#ifndef FEXTENDPOTENTIAL_HPP
#define FEXTENDPOTENTIAL_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FGlobal.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendPotential
* Please read the license
* This class is an extenssion.
* It proposes a Potential (FReal).
*/
class FExtendPotential {
protected:
    FReal potential;   //< The potential extended

public:
    /** Default constructor */
    FExtendPotential() : potential(0) {
    }

    /** Copy constructor */
    FExtendPotential(const FExtendPotential& other) : potential(other.potential) {
    }

    /** Destructor */
    virtual ~FExtendPotential(){
    }

    /** Copy operator */
    FExtendPotential& operator=(const FExtendPotential& other) {
        this->potential = other.potential;
        return *this;
    }

    /** To get the potential */
    FReal getPotential() const {
        return this->potential;
    }

    /** To set the potential */
    void setPotential(const FReal inPotential) {
        this->potential = inPotential;
    }

    /** To inc the potential */
    void incPotential(const FReal inPotential) {
        this->potential += inPotential;
    }

};


#endif //FEXTENDPOTENTIAL_HPP

// [--LICENSE--]
