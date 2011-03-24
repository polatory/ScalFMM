#ifndef FBASICPARTICULE_HPP
#define FBASICPARTICULE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FExtendPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicParticule
* Please read the license
*
* This class defines a basic particule used for examples. It extends
* the mininum, only what is needed by FOctree and FFMMAlgorithm
* to make the things working.
* By using this extension it will implement the FAbstractParticule without
* inheriting from it.
*/
class FBasicParticule : public FExtendPosition{
public:
    /** Default destructor */
    virtual ~FBasicParticule(){
    }
};


#endif //FBASICPARTICULE_HPP

// [--LICENSE--]
