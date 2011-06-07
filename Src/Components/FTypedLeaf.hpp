#ifndef FTYPEDLEAF_HPP
#define FTYPEDLEAF_HPP
// /!\ Please, you must read the license at the bottom of this page

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
template< class ParticleClass ,template <class ParticleClass> class ContainerClass>
class FTypedLeaf  : public FAbstractLeaf<ParticleClass,ContainerClass>, public FAssertable {
    ContainerClass<ParticleClass> sources;
    ContainerClass<ParticleClass> targets;

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
        else assert(false, "Error particle has undefined type.", __LINE__, __FILE__);
    }

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particles that are sources
        */
    ContainerClass<ParticleClass>* getSrc() {
        return &this->sources;
    }

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particles that are targets
        */
    ContainerClass<ParticleClass>* getTargets() {
        return &this->targets;
    }

};


#endif //FTYPEDLEAF_HPP

// [--LICENSE--]
