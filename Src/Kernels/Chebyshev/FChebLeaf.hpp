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
#ifndef FCHEBLEAF_HPP
#define FCHEBLEAF_HPP


#include "../../Utils/FNoCopyable.hpp"
#include "../../Containers/FVector.hpp"

class FChebParticle;

/**
 * @author Matthias Messner (matthias.messner@inria.fr)
 *
 * @class FChebLeaf
 *
 * @brief Please read the license
 *
 * This class is used as a leaf in the Chebyshev FMM approach.
 */
template<class ParticleClass, class ContainerClass>
class FChebLeaf : FNoCopyable
{
private:
	ContainerClass particles; //!< Stores all particles contained by this leaf
	
public:
	~FChebLeaf() {}
	
	/**
	 * @param particle The new particle to be added to the leaf
	 */
        void push(const ParticleClass& inParticle)
	{
		particles.push(inParticle);
	}
	
	/**
	 * @return Array containing source particles
	 */
	ContainerClass* getSrc()
	{
		return &particles;
	}

	/**
	 * @return Array containing target particles
	 */
	ContainerClass* getTargets()
	{
		return &particles;
	}

};


#endif //FSIMPLELEAF_HPP


