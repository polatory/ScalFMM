#ifndef FCONVERT_HPP
#define FCONVERT_HPP

#include <sstream>
#include <iostream>
#include <string.h>

#include "../Containers/FTreeCoordinate.hpp"
#include "../Utils/F3DPosition.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class
* Please read the license
* This class proposes some convertion functions
* FConvert is used here as a namespace.
*/
class FConvert {
public :
    /** To get spatial position (F3DPosition) from morton data
      * @param inIndex the morton index of the object
      * @param inLevel the level of the current object
      * @param inWidthAtLevel the width of the box at this level
      * return outPosition the result
      */
    static F3DPosition MortonToPosition(const MortonIndex inIndex, const int inLevel, const FReal inWidthAtLevel){
        FTreeCoordinate treePosition;
        treePosition.setPositionFromMorton(inIndex, inLevel);

        F3DPosition outPosition(
            treePosition.getX() * inWidthAtLevel + inWidthAtLevel/2,
            treePosition.getY() * inWidthAtLevel + inWidthAtLevel/2,
            treePosition.getZ() * inWidthAtLevel + inWidthAtLevel/2
        );

        return outPosition;
    }

    /**
    * This function gives a parameter in a standart type
    * @parameter inArg parameter position has to be strictly less than argc/userParemetersCount
    * @return argv[inArg] in the template VariableType form
    * @warning VariableType need to work with istream >> operator
    * <code> const int argInt = userParemetersAt<int>(1,-1); </code>
    */
    template <class VariableType>
    static const VariableType charToOther(const char* const str, const VariableType& defaultValue = VariableType()){
            std::istringstream iss(str,std::istringstream::in);
            VariableType value;
            iss >> value;
            if( /*iss.tellg()*/ iss.eof() ) return value;
            return defaultValue;
    }
};

#endif

