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

public:
    /**
    * Constructor with a position
    * @param inIndex the morton index
    */
    FBasicCell(const MortonIndex inIndex) : index(inIndex) {
    }


    /** Default constructor */
    FBasicCell(){
    }

    /** Default destructor */
    virtual ~FBasicCell(){
    }

    /**
    * From the FAbstractParticule definition
    * @return the position of the current cell
    */
    MortonIndex getMortonIndex() const {
        return index;
    }

    /**
      * This function is needed by the basic loader to fill the current particule
      * @param inPos the position given by the basic loader
      */
    void setMortonIndex(const MortonIndex inIndex) {
        index = inIndex;
    }
};


#endif //FBASICCELL_HPP

// [--LICENSE--]
