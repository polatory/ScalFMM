#ifndef FSIMPLELEAF_HPP
#define FSIMPLELEAF_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FAbstractLeaf.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FSimpleLeaf
* @brief
* Please read the license
* This class is used as a leaf in simple system (source AND target)
* here there is only one list that store all particles.
*/
template< class ParticleClass ,template <class ParticleClass> class ContainerClass>
class FSimpleLeaf : public FAbstractLeaf<ParticleClass,ContainerClass> {
    ContainerClass<ParticleClass> particles;

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
    ContainerClass<ParticleClass>* getSrc() {
        return &this->particles;
    }

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particles that are targets
        */
    ContainerClass<ParticleClass>* getTargets() {
        return &this->particles;
    }

};


#endif //FSIMPLELEAF_HPP

// [--LICENSE--]
