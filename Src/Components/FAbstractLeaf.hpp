#ifndef FABSTRACTLEAF_HPP
#define FABSTRACTLEAF_HPP
// /!\ Please, you must read the license at the bottom of this page


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractLeaf
* @brief
* Please read the license
* This class is used to enable the use of typed particles
* (source XOR target) or simple system (source AND target)
*/
template< class ParticleClass >
class FAbstractLeaf {
public:
    /** Default destructor */
    virtual ~FAbstractLeaf(){
    }

    /**
        * To add a new particle in the leaf
        * @param particle the new particle
        * Depending on the system to use the class that inherit
        * this interface can sort the particle as they like.
        */
    virtual void push(const ParticleClass& particle) = 0;

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particles that are sources
        * Depending on the system to use the class that inherit
        * this interface can sort the particle as they like.
        */
    virtual FList<ParticleClass>* getSrc() = 0;

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particles that are targets
        * Depending on the system to use the class that inherit
        * this interface can sort the particle as they like.
        */
    virtual FList<ParticleClass>* getTargets() = 0;

};


#endif //FABSTRACTLEAF_HPP

// [--LICENSE--]
