#ifndef FBASICPARTICULE_HPP
#define FBASICPARTICULE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FAbstractParticule.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicParticule
* Please read the license
*
* This class defines a basic particule used for examples.
*/
class FBasicParticule : public FAbstractParticule{
protected:
    F3DPosition pos;    //< Particule's position

public:
    /**
    * Constructor with a position
    * @param inX x position
    * @param inY y position
    * @param inZ z position
    */
    FBasicParticule(const double inX, const double inY, const double inZ) : pos(inX,inY,inZ) {
    }

    /**
    * Constructor with a position object
    * @param inPos particule position
    */
    FBasicParticule(const F3DPosition& inPos) : pos(inPos) {
    }

    /** Default constructor */
    FBasicParticule(){
    }

    /** Default destructor */
    virtual ~FBasicParticule(){
    }

    /**
    * From the FAbstractParticule definition
    * @return the position of the current cell
    */
    F3DPosition getPosition() const {
        return pos;
    }

    /**
      * This function is needed by the basic loader to fill the current particule
      * @param inPos the position given by the basic loader
      */
    void setPosition(const F3DPosition& inPos) {
        pos = inPos;
    }
};


#endif //FBASICPARTICULE_HPP

// [--LICENSE--]
