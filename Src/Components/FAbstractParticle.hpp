// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FABSTRACTPARTICLE_HPP
#define FABSTRACTPARTICLE_HPP


/* forward declaration to avoid include */
class FPoint;

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
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
        virtual const FPoint& getPosition() const = 0;
};


#endif //FABSTRACTPARTICLE_HPP


