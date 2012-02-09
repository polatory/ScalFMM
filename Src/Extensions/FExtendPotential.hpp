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
#ifndef FEXTENDPOTENTIAL_HPP
#define FEXTENDPOTENTIAL_HPP


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


