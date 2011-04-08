#ifndef FBASICKERNELS_HPP
#define FBASICKERNELS_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FAbstractKernels.hpp"

#include "../Utils/FTrace.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class AbstractKernels
* @brief
* Please read the license
*
* This kernels simply shows the details of the information
* it receives (in debug)
*/
template< class ParticuleClass, class CellClass, int TreeHeight>
class FBasicKernels : public FAbstractKernels<ParticuleClass,CellClass,TreeHeight> {
public:
    /** Default destructor */
    virtual ~FBasicKernels(){
    }

    /** When init the kernel */
    virtual void init(){}

    /** Print the number of particules */
    virtual void P2M(CellClass* const , const FList<ParticuleClass*>* const ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    /** Print the morton index */
    virtual void M2M(CellClass* const FRestrict , const CellClass*const FRestrict *const FRestrict , const int ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    /** Print the morton index */
    virtual void M2L(CellClass* const FRestrict , const CellClass*const FRestrict *const FRestrict , const int , const int ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    /** Print the morton index */
    virtual void L2L(const CellClass* const FRestrict , CellClass* FRestrict *const FRestrict  , const int ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    /** Print the number of particules */
    virtual void L2P(const CellClass* const , FList<ParticuleClass*>* const ){
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    /** Print the number of particules */
    virtual void P2P(FList<ParticuleClass*>* const FRestrict , const FList<ParticuleClass*>* const FRestrict ,
                     const FList<ParticuleClass*>* FRestrict const* FRestrict , const int ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }
};


#endif //FBASICKERNELS_HPP

// [--LICENSE--]
