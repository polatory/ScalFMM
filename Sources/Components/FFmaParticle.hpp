#ifndef FFMAPARTICLE_HPP
#define FFMAPARTICLE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FBasicParticle.hpp"
#include "../Extenssions/FExtendPhysicalValue.hpp"

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


#endif //FFMAPARTICLE_HPP

// [--LICENSE--]
