#ifndef FABSTRACTPARTICULE_HPP
#define FABSTRACTPARTICULE_HPP
// /!\ Please, you must read the license at the bottom of this page

/* forward declaration to avoid include */
class F3DPosition;

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractBody
* @brief
* Please read the license
*
* This class define the method that every particule class
* has to implement.
* @warning Inherite from this class when implement a specific particule type
*/
class FAbstractParticule{
public:	
	/** Default destructor */
	virtual ~FAbstractParticule(){
	}

	/**
	* Must be implemented by each user Particule class
	* @return the position of the current cell
	*/
	virtual F3DPosition getPosition() const = 0;
};


#endif //FABSTRACTPARTICULE_HPP

// [--LICENSE--]
