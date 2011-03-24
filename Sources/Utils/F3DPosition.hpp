#ifndef F3DPOSITION_HPP
#define F3DPOSITION_HPP
// /!\ Please, you must read the license at the bottom of this page

// To get memcpy
#include <cstring>


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class F3DPosition
* Please read the license
*
*/
class F3DPosition{
private:
	double x; //< x position
	double y; //< y position
	double z; //< z position

public:	
	/** Default constructor */
	F3DPosition() : x(0.0), y(0.0), z(0.0){
	}

	/** Default constructor */
	F3DPosition(const double inX,const double inY,const double inZ)
		: x(inX), y(inY), z(inZ){
	}

	/** Default destructor */
	virtual ~F3DPosition(){
	}

	/**
	* Copy constructor
	* @param other the source class to copy
	*/
        F3DPosition(const F3DPosition& other): x(other.x), y(other.y), z(other.z){
	}

	/**
	* Copy constructor
	* @param other the source class to copy
	* @return this a reference to the current class
	*/
	F3DPosition& operator=(const F3DPosition& other){
                memcpy(this, &other, sizeof(F3DPosition));
		return *this;
	}

	/**
	* Position setter
	* @param other the source class to copy
	* @return this a reference to the current class
	*/
	void setPosition(const double inX,const double inY,const double inZ){
		this->x = inX;
		this->y = inY;
		this->z = inZ;
	}

	/**
	* Get x
	* @return this->x
	*/
	double getX() const{
		return this->x;
	}

	/**
	* Get y
	* @return this->y
	*/
	double getY() const{
		return this->y;
	}

	/**
	* Get z
	* @return this->z
	*/
	double getZ() const{
		return this->z;
	}

	/**
	* Set x
	* @param the new x
	*/
	void setX(const double inX){
		this->x = inX;
	}

	/**
	* Set y
	* @param the new y
	*/
	void setY(const double inY){
		this->y = inY;
	}

	/**
	* Set z
	* @param the new z
	*/
	void setZ(const double inZ){
		this->z = inZ;
	}

	/**
	* Subtract to all dim the inValue
	* @param inValue the value to substract
	* @return the current object after being subtracted
	*/
	F3DPosition& operator-=(const double inValue){
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
	F3DPosition& operator+=(const double inValue){
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
        F3DPosition& operator*=(const double value){
                this->x *= value;
                this->y *= value;
                this->z *= value;
                return *this;
        }

};

/**
* Operator F3Position minus double
* This substract inValue to all dimensions of the inPosition
* @param inPosition the position to compute
* @param inValue the value to decrease/substract position
* @return the resulting position
*/
F3DPosition operator-(const F3DPosition& inPosition, const double inValue){
	F3DPosition position(inPosition);
	position -= inValue;
	return position;
}

/**
* Operator F3Position plus double
* This affect from inValue all dimensions of the inPosition
* @param inPosition the position to compute
* @param inValue the value to increase/affect position
* @return the resulting position
*/
F3DPosition operator+(const F3DPosition& inPosition, const double inValue){
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

#endif //F3DPOSITION_HPP

// [--LICENSE--]
