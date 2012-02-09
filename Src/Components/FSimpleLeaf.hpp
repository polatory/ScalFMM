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
        * @param particle the new particle
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


