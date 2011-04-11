#ifndef FABSTRACTLEAF_HPP
#define FABSTRACTLEAF_HPP
// /!\ Please, you must read the license at the bottom of this page


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractLeaf
* @brief
* Please read the license
* This class is used to enable the use of typed particules
* (source XOR target) or simple system (source AND target)
*/
template< class ParticuleClass >
class FAbstractLeaf {
public:
    /** Default destructor */
    virtual ~FAbstractLeaf(){
    }

    /**
        * To add a new particule in the leaf
        * @param particule the new particule
        * Depending on the system to use the class that inherit
        * this interface can sort the particule as they like.
        */
    virtual void push(ParticuleClass* const particule) = 0;

    /**
        * To get all the sources in a leaf
        * @return a pointer to the list of particules that are sources
        * Depending on the system to use the class that inherit
        * this interface can sort the particule as they like.
        */
    virtual FList<ParticuleClass*>* getSources() = 0;

    /**
        * To get all the target in a leaf
        * @return a pointer to the list of particules that are targets
        * Depending on the system to use the class that inherit
        * this interface can sort the particule as they like.
        */
    virtual FList<ParticuleClass*>* getTargets() = 0;

};


#endif //FABSTRACTLEAF_HPP

// [--LICENSE--]
