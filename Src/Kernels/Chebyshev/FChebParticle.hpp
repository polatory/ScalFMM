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
#ifndef FCHEBPARTICLE_HPP
#define FCHEBPARTICLE_HPP

#include <stdexcept>
#include <cassert>

#include "../../Extensions/FExtendPosition.hpp"
#include "../../Extensions/FExtendPhysicalValue.hpp"
#include "../../Extensions/FExtendPotential.hpp"
#include "../../Extensions/FExtendForces.hpp"

/**
 * @author Matthias Messner (matthias.matthias@inria.fr)
 * @class FChebParticle
 * Please read the license
 *
 * The class @p FChebParticle defines the particle used for the Chebyshev FMM
 * approach.
 */
class FChebParticle : public FExtendPosition,
											public FExtendPhysicalValue,
											public FExtendPotential,
											public FExtendForces
{
public:
	~FChebParticle() {}
};


#endif
