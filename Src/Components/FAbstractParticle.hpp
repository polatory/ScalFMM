// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FABSTRACTPARTICLE_HPP
#define FABSTRACTPARTICLE_HPP


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


