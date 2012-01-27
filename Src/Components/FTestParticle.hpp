#ifndef FTESTPARTICLE_HPP
#define FTESTPARTICLE_HPP
// [--License--]

#include "FBasicParticle.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTestParticle
* Please read the license
*
* This class is used in the FTestKernels, please
* look at this class to know whit it is.
*
* Particles just need the data down.
*/
class FTestParticle : public FBasicParticle {
protected:
    // To store data during downard pass
    long long int dataDown;
public:
    FTestParticle(): dataDown(0){
    }

    /** Default destructor */
    virtual ~FTestParticle(){
    }

    /** Get the down data */
    long long int getDataDown() const {
        return this->dataDown;
    }

    /** Set down data */
    void setDataDown(const long long int inData){
        this->dataDown = inData;
    }
};


#endif //FTESTPARTICLE_HPP

// [--END--]
