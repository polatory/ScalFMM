#ifndef FTESTKERNELS_HPP
#define FTESTKERNELS_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FAbstractKernels.hpp"

#include "../Containers/FList.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class AbstractKernels
* @brief
* Please read the license
*
* This kernels is a virtual kernels to validate that the fmm algorithm is
* correctly done on particules.
* It used FTestCell and FTestParticule
*/
template< class ParticuleClass, class CellClass>
class FTestKernels : public FAbstractKernels<ParticuleClass,CellClass> {
public:
    /** Default destructor */
    virtual ~FTestKernels(){
    }

    // Before upward
    void P2M(CellClass* const pole, const FList<ParticuleClass*>* const particules) {
        // the pole represents all particules under
        pole->setDataUp(particules->getSize());
    }
    // During upward
    void M2M(CellClass* const FRestrict pole, const CellClass *const FRestrict *const FRestrict child, const int inLevel) {
        // A parent represents the sum of the child
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child[idx]){
                pole->setDataUp(pole->getDataUp() + child[idx]->getDataUp());
            }
        }
    }
    // Before Downward
    void M2L(CellClass* const FRestrict pole, const CellClass*const FRestrict *const distantNeighbors, const int size, const int inLevel) {
        // The pole is impacted by what represent other poles
        for(int idx = 0 ; idx < size ; ++idx){
            pole->setDataDown(pole->getDataDown() + distantNeighbors[idx]->getDataUp());
        }
    }
    // During Downward
    void L2L(const CellClass*const FRestrict local, CellClass* FRestrict *const FRestrict child, const int inLevel) {
        // Each child is impacted by the father
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child[idx]){
                child[idx]->setDataDown(local->getDataDown() + child[idx]->getDataDown());
            }
        }
    }
    // After Downward
    void L2P(const CellClass* const  local, FList<ParticuleClass*>*const particules){
        // The particules is impacted by the parent cell
        typename FList<ParticuleClass*>::BasicIterator iter(*particules);
        while( iter.isValide() ){
            iter.value()->setDataDown(iter.value()->getDataDown() + local->getDataDown());
            iter.progress();
        }
    }
    // After Downward
    void P2P(FList<ParticuleClass*>*const FRestrict currentBox, const FList<ParticuleClass*>* FRestrict const* FRestrict directNeighbors, const int size) {
        // Each particules targeted is impacted by the particules sources
        long inc = currentBox->getSize() - 1;
        for(int idx = 0 ; idx < size ; ++idx){
            inc += directNeighbors[idx]->getSize();
        }

        typename FList<ParticuleClass*>::BasicIterator iter(*currentBox);
        while( iter.isValide() ){
            iter.value()->setDataDown(iter.value()->getDataDown() + inc);
            iter.progress();
        }
    }
};


#endif //FTESTKERNELS_HPP

// [--LICENSE--]
