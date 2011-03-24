#ifndef FFMAPARTICULE_HPP
#define FFMAPARTICULE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FBasicParticule.hpp"
#include "FExtendValue.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmaParticule
* Please read the license
*
* This class defines a particule for FMA loader.
* As defined in FFmaLoader it needs {FBasicParticule,FExtendValue}
*/
class FFmaParticule : public FBasicParticule, public FExtendValue {
public:
    /** Default destructor */
    virtual ~FFmaParticule(){
    }
};


#endif //FFMAPARTICULE_HPP

// [--LICENSE--]
