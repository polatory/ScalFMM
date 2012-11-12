// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Matthias Messner.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Matthias Messner.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
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
