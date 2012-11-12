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
#ifndef FSIMPLELEAF_HPP
#define FSIMPLELEAF_HPP


#include "FAbstractLeaf.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FSimpleLeaf
* @brief
* Please read the license
* This class is used as a leaf in simple system (source AND target)
* here there is only one list that store all particles.
*/
template< class ParticleClass , class ContainerClass>
class FSimpleLeaf : public FAbstractLeaf<ParticleClass,ContainerClass> {
    ContainerClass particles;

public:
    /** Default destructor */
    virtual ~FSimpleLeaf(){
    }

    /**
        * To add a new particle in the leaf
        * @param particle the new particle to store in the current leaf
        */
    void push(const ParticleClass& particle){
        this->particles.push(particle);
    }

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particles that are sources
        */
    ContainerClass* getSrc() {
        return &this->particles;
    }

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particles that are targets
        */
    ContainerClass* getTargets() {
        return &this->particles;
    }

};


#endif //FSIMPLELEAF_HPP


