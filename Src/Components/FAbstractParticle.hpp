#ifndef FABSTRACTPARTICLE_HPP
#define FABSTRACTPARTICLE_HPP
// /!\ Please, you must read the license at the bottom of this page

/* forward declaration to avoid include */
class F3DPosition;

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractBody
* @brief
* Please read the license
*
* This class define the method that every particle class
* has to implement.
*
* In fact FOctree & FFmmAlgorithm need this function to be implemented.
* But you cannot use this interface with the extension (as an example :
* because the compiler will faill to know if getPosition is coming
* from this interface or from the extension)
*
*
* @warning Inherite from this class when implement a specific particle type
*/
class FAbstractParticle{
public:	
	/** Default destructor */
	virtual ~FAbstractParticle(){
	}

	/**
	* Must be implemented by each user Particle class
	* @return the position of the current cell
	*/
        virtual const F3DPosition& getPosition() const = 0;
};


#endif //FABSTRACTPARTICLE_HPP

// [--LICENSE--]
