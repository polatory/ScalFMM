#ifndef FBASICCELL_HPP
#define FBASICCELL_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Extensions/FExtendPosition.hpp"
#include "../Extensions/FExtendMortonIndex.hpp"
#include "../Extensions/FExtendCoordinate.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicCell
* Please read the license
*
* This class defines a basic cell used for examples. It extends
* the mininum, only what is needed by FOctree and FFmmAlgorithm
* to make the things working.
* By using this extension it will implement the FAbstractCell without
* inheriting from it.
*/
class FBasicCell : public FExtendPosition, public FExtendMortonIndex, public FExtendCoordinate {
public:
    /** Default destructor */
    virtual ~FBasicCell(){
    }
};


#endif //FBASICCELL_HPP

// [--LICENSE--]
