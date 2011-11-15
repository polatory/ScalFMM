#ifndef F3DPOSITION_HPP
#define F3DPOSITION_HPP
// /!\ Please, you must read the license at the bottom of this page

// To get memcpy
#include <cstring>
#include "FGlobal.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class F3DPosition
* Please read the license
*
* This class is a 3D vector. It can be used as a position
* or as a 3d forces vector etc.
*/
class F3DPosition{
private:
    FReal x; //< x position
    FReal y; //< y position
    FReal z; //< z position

public:	
    /** Default constructor (sets position to 0/0/0) */
    F3DPosition() : x(0.0), y(0.0), z(0.0){
    }

    /** Constructor from values */
    explicit F3DPosition(const FReal inX,const FReal inY,const FReal inZ)
        : x(inX), y(inY), z(inZ){
    }

    /**
 * Copy constructor
 * @param other the source class to copy
 */
    F3DPosition(const F3DPosition& other): x(other.x), y(other.y), z(other.z){
    }

    /**
        * Assignement operator
        * @param other the source class to copy
        */
    F3DPosition(const F3DPosition& other, const FReal addset)
        : x(other.x + addset), y(other.y + addset), z(other.z + addset){
    }

    /**
 * Copy constructor
 * @param other the source class to copy
 * @return this a reference to the current class
 */
    F3DPosition& operator=(const F3DPosition& other){
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
        return *this;
    }

    /**
 * Position setter
 * @param other the source class to copy
 * @return this a reference to the current class
 */
    void setPosition(const FReal inX,const FReal inY,const FReal inZ){
        this->x = inX;
        this->y = inY;
        this->z = inZ;
    }

    /**
 * Get x
 * @return this->x
 */
    FReal getX() const{
        return this->x;
    }

    /**
 * Get y
 * @return this->y
 */
    FReal getY() const{
        return this->y;
    }

    /**
 * Get z
 * @return this->z
 */
    FReal getZ() const{
        return this->z;
    }

    /**
 * Set x
 * @param the new x
 */
    void setX(const FReal inX){
        this->x = inX;
    }

    /**
 * Set y
 * @param the new y
 */
    void setY(const FReal inY){
        this->y = inY;
    }

    /**
 * Set z
 * @param the new z
 */
    void setZ(const FReal inZ){
        this->z = inZ;
    }

    /**
 * Set x
 * @param the new x
 */
    void incX(const FReal inX){
        this->x += inX;
    }

    /**
 * Set y
 * @param the new y
 */
    void incY(const FReal inY){
        this->y += inY;
    }

    /**
 * Set z
 * @param the new z
 */
    void incZ(const FReal inZ){
        this->z += inZ;
    }

    /**
 * Subtract to all dim the inValue
 * @param inValue the value to substract
 * @return the current object after being subtracted
 */
    F3DPosition& operator-=(const FReal inValue){
        this->x -= inValue;
        this->y -= inValue;
        this->z -= inValue;
        return *this;
    }

    /**
 * Affect to all dim the inValue
 * @param inValue the value to afect
 * @return the current object after being affected
 */
    F3DPosition& operator+=(const FReal inValue){
        this->x += inValue;
        this->y += inValue;
        this->z += inValue;
        return *this;
    }

    /**
        * Subtract to all dim the other position
        * @param other the value to substract
        * @return the current object after being subtracted
        */
    F3DPosition& operator-=(const F3DPosition& other){
        this->x -= other.x;
        this->y -= other.y;
        this->z -= other.z;
        return *this;
    }

    /**
        * Affect to all dim the other position
        * @param other the value to afect
        * @return the current object after being affected
        */
    F3DPosition& operator+=(const F3DPosition& other){
        this->x += other.x;
        this->y += other.y;
        this->z += other.z;
        return *this;
    }

    /**
        * Affect to all dim the other position
        * @param other the value to afect
        * @return the current object after being affected
        */
    F3DPosition& operator*=(const FReal value){
        this->x *= value;
        this->y *= value;
        this->z *= value;
        return *this;
    }

    /**
        * Operator F3Position minus FReal
        * This substract inValue to all dimensions of the inPosition
        * @param inPosition the position to compute
        * @param inValue the value to decrease/substract position
        * @return the resulting position
        */
    friend inline F3DPosition operator-(const F3DPosition& inPosition, const FReal inValue){
        return F3DPosition(inPosition, -inValue);
    }

    /**
        * Operator F3Position plus FReal
        * This affect from inValue all dimensions of the inPosition
        * @param inPosition the position to compute
        * @param inValue the value to increase/affect position
        * @return the resulting position
        */
    friend inline F3DPosition operator+(const F3DPosition& inPosition, const FReal inValue){
        return F3DPosition(inPosition, inValue);
    }

    /**
        * Operator F3Position minus F3Position
        * This substract one from anther
        * @param inPosition the position to reduce
        * @param inOther the position to decrease/substract inPosition
        * @return the resulting position
        */
    friend inline F3DPosition operator-(const F3DPosition& inPosition, const F3DPosition& inOther){
        return F3DPosition(inPosition.x - inOther.x, inPosition.y - inOther.y, inPosition.z - inOther.z);
    }

    /**
        * Operator F3Position plus F3Position
        * This substract one from anther
        * @param inPosition the position to reduce
        * @param inOther the position to increase inPosition
        * @return the resulting position
        */
    friend inline F3DPosition operator+(const F3DPosition& inPosition, const F3DPosition& inOther){
        return F3DPosition(inPosition.x + inOther.x, inPosition.y + inOther.y, inPosition.z + inOther.z);
    }
};





#endif //F3DPOSITION_HPP

// [--LICENSE--]
