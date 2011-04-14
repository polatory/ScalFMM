#ifndef FSIMPLELEAF_HPP
#define FSIMPLELEAF_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Containers/FList.hpp"
#include "FAbstractLeaf.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FSimpleLeaf
* @brief
* Please read the license
* This class is used as a leaf in simple system (source AND target)
* here there is only one list that store all particles.
*/
template< class ParticleClass >
class FSimpleLeaf : public FAbstractLeaf<ParticleClass> {
    FList<ParticleClass*> particles;

public:
    /** Default destructor */
    virtual ~FSimpleLeaf(){
    }

    /**
        * To add a new particle in the leaf
        * @param particle the new particle
        */
    void push(ParticleClass* const particle){
        this->particles.pushFront(particle);
    }

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particles that are sources
        */
    FList<ParticleClass*>* getSources() {
        return &this->particles;
    }

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particles that are targets
        */
    FList<ParticleClass*>* getTargets() {
        return &this->particles;
    }

};


#endif //FSIMPLELEAF_HPP

// [--LICENSE--]
