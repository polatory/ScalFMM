#ifndef FFMAPARTICULE_HPP
#define FFMAPARTICULE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FBasicParticule.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmaParticule
* Please read the license
*
* This class defines a basic particule used for examples.
*/
class FFmaParticule : public FBasicParticule{
protected:
    double value;    //< Particule's fma data

public:
    /**
    * Constructor with a position
    * @param inX x position
    * @param inY y position
    * @param inZ z position
    * @param inValue particule value
    */
    FFmaParticule(const double inX, const double inY, const double inZ, const double inValue)
	: FBasicParticule(inX,inY,inZ), value(inValue) {
    }

    /**
    * Constructor with a position object
    * @param inPos particule position
    */
    FFmaParticule(const F3DPosition& inPos) : FBasicParticule(inPos), value(0){
    }

    /** Default constructor */
    FFmaParticule() : value(0) {
    }

    /** Default destructor */
    virtual ~FFmaParticule(){
    }

    /**
    * From the Fmb needed
    * @return the value of the current cell
    */
    double getValue() const {
        return this->value;
    }

    /**
      * This function is needed by the fma loader to fill the current particule
      * @param inValue the data given by the fma loader
      */
    void setValue(const double inValue) {
        this->value = inValue;
    }
};


#endif //FFMAPARTICULE_HPP

// [--LICENSE--]
