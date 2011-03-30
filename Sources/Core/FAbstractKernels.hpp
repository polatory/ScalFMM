#ifndef FABSTRACTKERNELS_HPP
#define FABSTRACTKERNELS_HPP
// /!\ Please, you must read the license at the bottom of this page


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractKernels
* @brief
* Please read the license
*
* If you want to create you own kernels you have to inherit from this class.
*/
template< class ParticuleClass, class CellClass>
class FAbstractKernels{
public:
    /** Default destructor */
    virtual ~FAbstractKernels(){
    }

    /**
        * Init the Kernels
        * This function is called just before to start the computation
        */
    virtual void init(){}

    /**
        * P2M
        * particules to multipole
        * @param pole the multipole to fill using the particules
        * @param particules the particules from the same spacial boxe
        */
    virtual void P2M(CellClass* const pole, const FList<ParticuleClass*>* const particules) = 0;

    /**
        * M2M
        * Multipole to multipole
        * @param pole the father (the boxe that contains other ones)
        * @param child the boxe to take values from
        * @param the current computation level
        * the child array has a size of 8 elements (address if exists or 0 otherwise).
        * You must test if a pointer is 0 to know if an element exists inside this array
        */
    virtual void M2M(CellClass* const pole, const CellClass*const* const child, const int inLevel) = 0;

    /**
        * M2L
        * Multipole to local
        * @param local the element to fill using distant neighbors
        * @param distantNeighbors is an array containing fathers's direct neighbors's child - direct neigbors
        * @param size the number of neighbors
        * @param inLevel the current level of the computation
        */
    virtual void M2L(CellClass* const local, const CellClass*const* const distantNeighbors, const int size, const int inLevel) = 0;

    /**
        * L2L
        * Local to local
        * @param the father to take value from
        * @param the child to downward values (child may have already been impacted by M2L)
        * @param level the current level of computation
        * the child array has a size of 8 elements (address if exists or 0 otherwise).
        * You must test if a pointer is 0 to know if an element exists inside this array
        */
    virtual void L2L(const CellClass* const local, CellClass** const child, const int inLevel) = 0;

    /**
        * L2P
        * Local to particules
        * @param local the leaf element (smaller boxe local element)
        * @param particules the list of particules inside this boxe
        */
    virtual void L2P(const CellClass* const local, FList<ParticuleClass*>* const particules) = 0;

    /**
        * P2P
        * Particules to particules
        * @param particules current boxe particules
        * @param directNeighborsParticules the particules from direct neighbors (this is an array of list)
        * @param size the number of direct neighbors (the size of the array directNeighborsParticules)
        */
    virtual void P2P(FList<ParticuleClass*>* const particules, const FList<ParticuleClass*>*const* const directNeighborsParticules, const int size) = 0;
};


#endif //FABSTRACTKERNELS_HPP

// [--LICENSE--]
