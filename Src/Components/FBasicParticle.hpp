#ifndef FBASICPARTICLE_HPP
#define FBASICPARTICLE_HPP
// [--License--]

#include "../Extensions/FExtendPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicParticle
* Please read the license
*
* This class defines a basic particle used for examples. It extends
* the mininum, only what is needed by FOctree and FFmmAlgorithm
* to make the things working.
* By using this extension it will implement the FAbstractParticle without
* inheriting from it.
*/
class FBasicParticle : public FExtendPosition{
public:
    /** Default destructor */
    virtual ~FBasicParticle(){
    }
};


#endif //FBASICPARTICLE_HPP

// [--END--]
