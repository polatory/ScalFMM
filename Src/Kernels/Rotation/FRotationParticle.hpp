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
#ifndef FROTATIONLPARTICLE_HPP
#define FROTATIONLPARTICLE_HPP

#include "../../Extensions/FExtendForces.hpp"
#include "../../Extensions/FExtendPotential.hpp"
#include "../../Extensions/FExtendParticleType.hpp"
#include "../../Components/FFmaParticle.hpp"

#include "../../Components/FAbstractSerializable.hpp"

class FRotationParticle : public FExtendForces, public FFmaParticle, public FExtendPotential {
public:

    /** Save current object */
    void save(FBufferWriter& buffer) const {
        FExtendForces::save(buffer);
        FFmaParticle::save(buffer);
        FExtendPotential::save(buffer);
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        FExtendForces::restore(buffer);
        FFmaParticle::restore(buffer);
        FExtendPotential::restore(buffer);
    }
};


class FTypedRotationParticle : public FRotationParticle, public FExtendParticleType {
public:
    /** Save current object */
    void save(FBufferWriter& buffer) const {
        FRotationParticle::save(buffer);
        FExtendParticleType::save(buffer);
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        FRotationParticle::restore(buffer);
        FExtendParticleType::restore(buffer);
    }
};

#endif // FROTATIONLPARTICLE_HPP
