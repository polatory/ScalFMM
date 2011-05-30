#ifndef FTYPEDLEAF_HPP
#define FTYPEDLEAF_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Containers/FList.hpp"
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
template< class ParticleClass >
class FTypedLeaf  : public FAbstractLeaf<ParticleClass>, public FAssertable {
    FList<ParticleClass> sources;
    FList<ParticleClass> targets;

public:
    /** Default destructor */
    virtual ~FTypedLeaf(){
    }

    /**
        * To add a new particle in the leaf
        * @param particle the new particle
        */
    void push(const ParticleClass& particle){
        if(particle.isTarget()) this->targets.pushFront(particle);
        else if(particle.isSource()) this->sources.pushFront(particle);
        else assert(false, "Error particle has undefined type.", __LINE__, __FILE__);
    }

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particles that are sources
        */
    FList<ParticleClass>* getSrc() {
        return &this->sources;
    }

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particles that are targets
        */
    FList<ParticleClass>* getTargets() {
        return &this->targets;
    }

};


#endif //FTYPEDLEAF_HPP

// [--LICENSE--]
