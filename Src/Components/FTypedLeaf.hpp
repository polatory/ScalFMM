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
#ifndef FTYPEDLEAF_HPP
#define FTYPEDLEAF_HPP


#include "../Utils/FAssertable.hpp"
#include "FAbstractLeaf.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTypedLeaf
* @brief
* Please read the license
* This class is used to enable the use of typed particles
* (source XOR target) or simple system (source AND target).
*
* Particles should be typed to enable targets/sources difference.
*/
template< class ParticleClass , class ContainerClass>
class FTypedLeaf  : public FAbstractLeaf<ParticleClass,ContainerClass>, public FAssertable {
    ContainerClass sources; //< The sources containers
    ContainerClass targets; //< The targets containers

public:
    /** Default destructor */
    virtual ~FTypedLeaf(){
    }

    /**
        * To add a new particle in the leaf
        * @param particle the new particle
        */
    void push(const ParticleClass& particle){
        // Test if is source or target particle
        if(particle.isTarget()) this->targets.push(particle);
        else if(particle.isSource()) this->sources.push(particle);
        else fassert(false, "Error particle has undefined type.", __LINE__, __FILE__);
    }

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particles that are sources
        */
    ContainerClass* getSrc() {
        return &this->sources;
    }

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particles that are targets
        */
    ContainerClass* getTargets() {
        return &this->targets;
    }

};


#endif //FTYPEDLEAF_HPP


