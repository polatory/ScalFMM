#ifndef FExtendPhysicalValue_HPP
#define FExtendPhysicalValue_HPP
// /!\ Please, you must read the license at the bottom of this page

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

    /** Destructor */
    virtual ~FExtendPhysicalValue(){
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

// [--LICENSE--]
