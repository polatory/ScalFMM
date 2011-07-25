#include "F3DPosition.hpp"

// Binary operator for F3DPosition Class


/**
* Operator F3Position minus FReal
* This substract inValue to all dimensions of the inPosition
* @param inPosition the position to compute
* @param inValue the value to decrease/substract position
* @return the resulting position
*/
F3DPosition operator-(const F3DPosition& inPosition, const FReal inValue){
	F3DPosition position(inPosition);
	position -= inValue;
	return position;
}

/**
* Operator F3Position plus FReal
* This affect from inValue all dimensions of the inPosition
* @param inPosition the position to compute
* @param inValue the value to increase/affect position
* @return the resulting position
*/
F3DPosition operator+(const F3DPosition& inPosition, const FReal inValue){
	F3DPosition position(inPosition);
	position += inValue;
	return position;
}

/**
* Operator F3Position minus F3Position
* This substract one from anther
* @param inPosition the position to reduce
* @param inOther the position to decrease/substract inPosition
* @return the resulting position
*/
F3DPosition operator-(const F3DPosition& inPosition, const F3DPosition& inOther){
        F3DPosition position(inPosition);
        position -= inOther;
        return position;
}

/**
* Operator F3Position plus F3Position
* This substract one from anther
* @param inPosition the position to reduce
* @param inOther the position to increase inPosition
* @return the resulting position
*/
F3DPosition operator+(const F3DPosition& inPosition, const F3DPosition& inOther){
        F3DPosition position(inPosition);
        position += inOther;
        return position;
}


