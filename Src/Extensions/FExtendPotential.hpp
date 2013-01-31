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
#ifndef FEXTENDPOTENTIAL_HPP
#define FEXTENDPOTENTIAL_HPP


#include "../Utils/FGlobal.hpp"
#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"

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

    /** Save current object */
    void save(FBufferWriter& buffer) const {
        buffer << potential;
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        buffer >> potential;
    }
};


#endif //FEXTENDPOTENTIAL_HPP


