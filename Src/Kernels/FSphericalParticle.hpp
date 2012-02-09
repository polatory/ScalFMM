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
#ifndef FSPHERICALPARTICLE_HPP
#define FSPHERICALPARTICLE_HPP

#include "../Extensions/FExtendForces.hpp"
#include "../Extensions/FExtendPotential.hpp"
#include "../Extensions/FExtendParticleType.hpp"
#include "../Components/FFmaParticle.hpp"

class FSphericalParticle : public FExtendForces, public FFmaParticle, public FExtendPotential {
public:
};

class FTypedSphericalParticle : public FSphericalParticle, public FExtendParticleType {
public:
};

#endif // FSPHERICALPARTICLE_HPP
