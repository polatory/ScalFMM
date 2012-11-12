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
#ifndef FExtendPhysicalValue_HPP
#define FExtendPhysicalValue_HPP


#include "../Utils/FGlobal.hpp"
#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"

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

    /** Save current object */
    void save(FBufferWriter& buffer) const {
        buffer << physicalValue;
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        buffer >> physicalValue;
    }

};


#endif //FExtendPhysicalValue_HPP


