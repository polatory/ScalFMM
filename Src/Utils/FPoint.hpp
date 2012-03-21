// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietardata[1] license.
// ===================================================================================
#ifndef FPOINT_HPP
#define FPOINT_HPP


// To get memcpy
#include <cstring>
#include <iostream>

#include "FGlobal.hpp"
#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FPoint
* Please read the license
*
* This class is a 3D vector. It can be used as a position
* or as a 3d forces vector etc.
*/
class FPoint {
private:
    FReal data[3]; //< all positions x y z

public:	
    /** Default constructor (sets position to 0/0/0) */
    FPoint(){
        data[0] = data[1] = data[2] = FReal(0);
    }

    /** Constructor from values */
    explicit FPoint(const FReal inX,const FReal inY,const FReal inZ){
        data[0] = inX;
        data[1] = inY;
        data[2] = inZ;
    }

    /**
 * Copdata[1] constructor
 * @param other the source class to copy
 */
    FPoint(const FPoint& other) {
        data[0] = other.data[0];
        data[1] = other.data[1];
        data[2] = other.data[2];
    }

    /**
        * Assignement operator
        * @param other the source class to copy
        */
    FPoint(const FPoint& other, const FReal addset) {
        data[0] = other.data[0] + addset;
        data[1] = other.data[1] + addset;
        data[2] = other.data[2] + addset;
    }

    /**
 * Copdata[1] constructor
 * @param other the source class to copy
 * @return this a reference to the current class
 */
    FPoint& operator=(const FPoint& other){
        this->data[0] = other.data[0];
        this->data[1] = other.data[1];
        this->data[2] = other.data[2];
        return *this;
    }

    /**
 * Position setter
 * @param other the source class to copy
 * @return this a reference to the current class
 */
    void setPosition(const FReal inX,const FReal inY,const FReal inZ){
        this->data[0] = inX;
        this->data[1] = inY;
        this->data[2] = inZ;
    }

    /**
 * Get x
 * @return this->data[0]
 */
    FReal getX() const{
        return this->data[0];
    }

    /**
 * Get y
 * @return this->data[1]
 */
    FReal getY() const{
        return this->data[1];
    }

    /**
 * Get z
 * @return this->data[2]
 */
    FReal getZ() const{
        return this->data[2];
    }

    /**
 * Set x
 * @param the new x
 */
    void setX(const FReal inX){
        this->data[0] = inX;
    }

    /**
 * Set y
 * @param the new y
 */
    void setY(const FReal inY){
        this->data[1] = inY;
    }

    /**
 * Set z
 * @param the new z
 */
    void setZ(const FReal inZ){
        this->data[2] = inZ;
    }

    /**
 * Set x
 * @param the new x
 */
    void incX(const FReal inX){
        this->data[0] += inX;
    }

    /**
 * Set y
 * @param the new y
 */
    void incY(const FReal inY){
        this->data[1] += inY;
    }

    /**
 * Set z
 * @param the new z
 */
    void incZ(const FReal inZ){
        this->data[2] += inZ;
    }

    /**
 * Subtract to all dim the inValue
 * @param inValue the value to substract
 * @return the current object after being subtracted
 */
    FPoint& operator-=(const FReal inValue){
        this->data[0] -= inValue;
        this->data[1] -= inValue;
        this->data[2] -= inValue;
        return *this;
    }

    /**
 * Affect to all dim the inValue
 * @param inValue the value to afect
 * @return the current object after being affected
 */
    FPoint& operator+=(const FReal inValue){
        this->data[0] += inValue;
        this->data[1] += inValue;
        this->data[2] += inValue;
        return *this;
    }

    /**
        * Subtract to all dim the other position
        * @param other the value to substract
        * @return the current object after being subtracted
        */
    FPoint& operator-=(const FPoint& other){
        this->data[0] -= other.data[0];
        this->data[1] -= other.data[1];
        this->data[2] -= other.data[2];
        return *this;
    }

    /**
        * Affect to all dim the other position
        * @param other the value to afect
        * @return the current object after being affected
        */
    FPoint& operator+=(const FPoint& other){
        this->data[0] += other.data[0];
        this->data[1] += other.data[1];
        this->data[2] += other.data[2];
        return *this;
    }

    /**
        * Affect to all dim the other position
        * @param other the value to afect
        * @return the current object after being affected
        */
    FPoint& operator*=(const FReal value){
        this->data[0] *= value;
        this->data[1] *= value;
        this->data[2] *= value;
        return *this;
    }

    /**
        * Operator F3Position minus FReal
        * This substract inValue to all dimensions of the inPosition
        * @param inPosition the position to compute
        * @param inValue the value to decrease/substract position
        * @return the resulting position
        */
    friend inline FPoint operator-(const FPoint& inPosition, const FReal inValue){
        return FPoint(inPosition, -inValue);
    }

    /**
        * Operator F3Position plus FReal
        * This affect from inValue all dimensions of the inPosition
        * @param inPosition the position to compute
        * @param inValue the value to increase/affect position
        * @return the resulting position
        */
    friend inline FPoint operator+(const FPoint& inPosition, const FReal inValue){
        return FPoint(inPosition, inValue);
    }

    /**
        * Operator F3Position minus F3Position
        * This substract one from anther
        * @param inPosition the position to reduce
        * @param inOther the position to decrease/substract inPosition
        * @return the resulting position
        */
    friend inline FPoint operator-(const FPoint& inPosition, const FPoint& inOther){
        return FPoint(inPosition.data[0] - inOther.data[0], inPosition.data[1] - inOther.data[1], inPosition.data[2] - inOther.data[2]);
    }

    /**
        * Operator F3Position plus F3Position
        * This substract one from anther
        * @param inPosition the position to reduce
        * @param inOther the position to increase inPosition
        * @return the resulting position
        */
    friend inline FPoint operator+(const FPoint& inPosition, const FPoint& inOther){
        return FPoint(inPosition.data[0] + inOther.data[0], inPosition.data[1] + inOther.data[1], inPosition.data[2] + inOther.data[2]);
    }

    /**
     * Operator stream FPoint to std::ostream
     * This can be used to simpldata[1] write out a position
     * @param[in,out] output where to write the position
     * @param[in] inPosition the position to write out
     * @return the output for multiple << operators
     */
    friend std::ostream& operator<<(std::ostream& output, const FPoint& inPosition){
        output << "(" <<  inPosition.getX() << ", " << inPosition.getY() << ", " << inPosition.getZ() <<")";
        return output;  // for multiple << operators.
    }

    /** Save current object */
    void save(FBufferWriter& buffer) const {
        buffer << data[0] << data[1] << data[2];
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        buffer >> data[0] >> data[1] >> data[2];
    }
};




#endif //FPOINT_HPP

