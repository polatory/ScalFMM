#ifndef FCHEBPARTICLE_HPP
#define FCHEBPARTICLE_HPP

#include <stdexcept>
#include <cassert>

#include "../Extensions/FExtendPosition.hpp"
#include "../Extensions/FExtendPhysicalValue.hpp"

/**
 * @author Matthias Messner (matthias.matthias@inria.fr)
 * @class FChebParticle
 * Please read the license
 *
 * The class @p FChebParticle defines the particle used for the Chebyshev FMM
 * approach.
 */
class FChebParticle : public FExtendPosition,
											public FExtendPhysicalValue
{
public:
	~FChebParticle() {}
};


#endif
