#ifndef FTESTPARTICULE_HPP
#define FTESTPARTICULE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FBasicParticule.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTestParticule
* Please read the license
*
* This class is used in the FTestKernels, please
* look at this class to know whit it is.
*/
class FTestParticule : public FBasicParticule {
protected:
    // To store data during downard pass
    long dataDown;
public:
    FTestParticule(): dataDown(0){
    }

    /** Default destructor */
    virtual ~FTestParticule(){
    }

    long getDataDown() const {
        return this->dataDown;
    }

    void setDataDown(const long inData){
        this->dataDown = inData;
    }
};


#endif //FTESTPARTICULE_HPP

// [--LICENSE--]
