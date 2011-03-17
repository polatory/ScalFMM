#ifndef FBASICCELL_HPP
#define FBASICCELL_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FAbstractCell.hpp"
#include "../Containers/FTreeCoordinate.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicCell
* Please read the license
*
* This class defines a basic cell used for examples.
*/
class FBasicCell : public FAbstractCell {
protected:
    MortonIndex index;    //< Cell's position in the tree
    F3DPosition position; //< Cell's spacial position

public:
    /** Default constructor */
    FBasicCell() : index(0), position(0,0,0){
    }

    /** Default destructor */
    virtual ~FBasicCell(){
    }

    /**
    * From the FAbstractParticule definition
    * @return the position of the current cell
    */
    MortonIndex getMortonIndex() const {
        return this->index;
    }

    /**
      * This function is needed by the basic loader to fill the current particule
      * @param inPos the position given by the basic loader
      */
    void setMortonIndex(const MortonIndex inIndex) {
        this->index = inIndex;
    }

    /**
    * @return the position of the current cell
    */
    F3DPosition getPosition() const{
        return this->position;
    }

    /**
    * @param inPosition the position of the current cell
    */
    void setPosition(const F3DPosition& inPosition){
        this->position = inPosition;
    }
};


#endif //FBASICCELL_HPP

// [--LICENSE--]
