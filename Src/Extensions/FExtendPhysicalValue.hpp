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
#ifndef FExtendPhysicalValue_HPP
#define FExtendPhysicalValue_HPP


#include "../Utils/FGlobal.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendPhysicalValue
* Please read the license
* This class is an extenssion.
* It proposes a physicalValue (FReal).
*/
class FExtendPhysicalValue {
protected:
    FReal physicalValue;   //< A simple physicalValue

public:
    /** Default constructor */
    FExtendPhysicalValue() : physicalValue(0) {
    }

    /** Copy constructor */
    FExtendPhysicalValue(const FExtendPhysicalValue& other) : physicalValue(other.physicalValue) {
    }

    /** Copy Constructor */
    FExtendPhysicalValue& operator=(const FExtendPhysicalValue& other) {
        this->physicalValue = other.physicalValue;
        return *this;
    }

    /** To get the physicalValue */
    FReal getPhysicalValue() const {
        return this->physicalValue;
    }

    /** To set the physicalValue */
    void setPhysicalValue(const FReal inphysicalValue) {
        this->physicalValue = inphysicalValue;
    }

};


#endif //FExtendPhysicalValue_HPP


