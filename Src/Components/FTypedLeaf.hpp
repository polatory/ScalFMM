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
* (source XOR target) or simple system (source AND target)
*/
template< class ParticleClass , class ContainerClass>
class FTypedLeaf  : public FAbstractLeaf<ParticleClass,ContainerClass>, public FAssertable {
    ContainerClass sources;
    ContainerClass targets;

public:
    /** Default destructor */
    virtual ~FTypedLeaf(){
    }

    /**
        * To add a new particle in the leaf
        * @param particle the new particle
        */
    void push(const ParticleClass& particle){
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


