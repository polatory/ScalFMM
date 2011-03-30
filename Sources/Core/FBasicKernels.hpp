#ifndef FBASICKERNELS_HPP
#define FBASICKERNELS_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FAbstractKernels.hpp"

#include "../Utils/FDebug.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class AbstractKernels
* @brief
* Please read the license
*
* This kernels simply shows the details of the information
* it receives (in debug)
*/
template< class ParticuleClass, class CellClass>
class FBasicKernels : public FAbstractKernels<ParticuleClass,CellClass> {
public:
    /** Default destructor */
    virtual ~FBasicKernels(){
    }

    /** When init the kernel */
    virtual void init(){}

    /** Print the number of particules */
    virtual void P2M(CellClass* const pole, const FList<ParticuleClass*>* const particules) {
        //FDEBUG( FDebug::Controller << "P2M : " << particules->getSize() << "\n" );
    }

    /** Print the morton index */
    virtual void M2M(CellClass* const FRestrict pole, const CellClass*const FRestrict *const FRestrict child, const int inLevel) {
        //FDEBUG( FDebug::Controller << "M2M : " << pole->getMortonIndex() << "\n" );
    }

    /** Print the morton index */
    virtual void M2L(CellClass* const FRestrict pole, const CellClass*const FRestrict *const FRestrict distantNeighbors, const int size, const int inLevel) {
        //FDEBUG( FDebug::Controller << "M2L : " << pole->getMortonIndex() << " (" << size << ")\n" );
    }

    /** Print the morton index */
    virtual void L2L(const CellClass* const FRestrict local, CellClass* FRestrict *const FRestrict  child, const int inLevel) {
        //FDEBUG( FDebug::Controller << "L2L : " << local->getMortonIndex() << "\n" );
    }

    /** Print the number of particules */
    virtual void L2P(const CellClass* const pole, FList<ParticuleClass*>* const particules){
        //FDEBUG( FDebug::Controller << "L2P : " << particules->getSize() << "\n" );
    }

    /** Print the number of particules */
    virtual void P2P(FList<ParticuleClass*>* const FRestrict currentBox, const FList<ParticuleClass*>* FRestrict const* FRestrict directNeighbors, const int size) {
        //FDEBUG( FDebug::Controller << "P2P : " << currentBox->getSize() << " (" << size << ")\n" );
    }
};


#endif //FBASICKERNELS_HPP

// [--LICENSE--]
