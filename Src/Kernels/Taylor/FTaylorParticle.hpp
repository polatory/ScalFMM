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
#ifndef FTAYLORPARTICLE_HPP
#define FTAYLORPARTICLE_HPP

#include "../../Extensions/FExtendForces.hpp"
#include "../../Extensions/FExtendPotential.hpp"
#include "../../Extensions/FExtendPosition.hpp"

/** 
 * @author Cyrille Piacibello
 * @class FTaylorParticle
 *
 * This class is a particle used for the Taylor Kernel
 */

class FTaylorParticle : public FExtendPotential, public FFmaParticle, public FExtendPosition {
  
public:
  
  /** Save current object */
  void save(FBufferWriter& buffer) const {
    FExtendPosition::save(buffer);
    FFmaParticle::save(buffer);
    FExtendPotential::save(buffer);
  }
  /** Retrieve current object */
  void restore(FBufferReader& buffer) {
    FExtendPosition::restore(buffer);
    FFmaParticle::restore(buffer);
    FExtendPotential::restore(buffer);
  }
  
  
};

#endif
