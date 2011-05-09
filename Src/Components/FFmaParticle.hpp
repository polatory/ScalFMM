#ifndef FFmaPARTICLE_HPP
#define FFmaPARTICLE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FBasicParticle.hpp"
#include "../Extensions/FExtendPhysicalValue.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmaParticle
* Please read the license
*
* This class defines a particle for FMA loader.
* As defined in FFmaLoader it needs {FBasicParticle,FExtendPhysicalValue}
*/
class FFmaParticle : public FBasicParticle, public FExtendPhysicalValue {
public:
    /** Default destructor */
    virtual ~FFmaParticle(){
    }
};


#endif //FFmaPARTICLE_HPP

// [--LICENSE--]
