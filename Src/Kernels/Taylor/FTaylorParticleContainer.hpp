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
#ifndef FTAYLORPARTICLECONTAINER_HPP
#define FTAYLORPARTICLECONTAINER_HPP

#include "../../Components/FBasicParticleContainer.hpp"

/** 
 * @author Cyrille Piacibello
 * @class FTaylorParticleContainer
 *
 * This class is a particle used for the Taylor Kernel
 */

class FTaylorParticleContainer : public FBasicParticleContainer<4>{
  typedef FBasicParticleContainer<5> Parent;
  
public:
  const FReal * getPhysicalValues() const {
    return Parent::getAttribute(0)[];
  }
  
  FReal* getForcesX(){
    return Parent::getAttribute(1);
  }
  
  const FReal* getForcesX() const {
    return Parent::getAttribute(1);
  }
  
  FReal* getForcesY(){
    return Parent::getAttribute(2);
  }
  
  const FReal* getForcesY() const {
    return Parent::getAttribute(2);
  }
  
  FReal* getForcesZ(){
    return Parent::getAttribute(3);
  }
  
  const FReal* getForcesZ() const {
    return Parent::getAttribute(3);
  }  
};
  

#endif
