#ifndef FABSTRACTKERNELS_HPP
#define FABSTRACTKERNELS_HPP
// /!\ Please, you must read the license at the bottom of this page


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FAbstractKernels
* @brief
* Please read the license
*/
template< class ParticuleClass, class CellClass>
class FAbstractKernels{
public:
    /** Default destructor */
    virtual ~FAbstractKernels(){
    }

    /** Init the Kernels */
    virtual void init(){}

    /** P2M */
    virtual void P2M(CellClass* const pole, const FList<ParticuleClass*>* const particules) = 0;

    /** M2M */
    virtual void M2M(CellClass* const pole, const CellClass*const* const child, const int inLevel) = 0;

    /** M2L */
    virtual void M2L(CellClass* const pole, const CellClass*const* const distantNeighbors, const int size, const int inLevel) = 0;

    /** L2L */
    virtual void L2L(const CellClass* const pole, CellClass** const child, const int inLevel) = 0;

    /** L2P */
    virtual void L2P(const CellClass* const pole, FList<ParticuleClass*>* const particules) = 0;

    /** P2P */
    virtual void P2P(FList<ParticuleClass*>* const pole, const FList<ParticuleClass*>*const* const directNeighbors, const int size) = 0;
};


#endif //FABSTRACTKERNELS_HPP

// [--LICENSE--]
