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
template< class ParticleClass, class CellClass, int TreeHeight>
class FBasicKernels : public FAbstractKernels<ParticleClass,CellClass,TreeHeight> {
public:
    /** Default destructor */
    virtual ~FBasicKernels(){
    }

    /** Print the number of particles */
    virtual void P2M(CellClass* const , const FList<ParticleClass*>* const ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    /** Print the morton index */
    virtual void M2M(CellClass* const FRestrict , const CellClass*const FRestrict *const FRestrict , const int ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    /** Print the morton index */
    virtual void M2L(CellClass* const FRestrict , const CellClass* [], const int , const int ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    /** Print the morton index */
    virtual void L2L(const CellClass* const FRestrict , CellClass* FRestrict *const FRestrict  , const int ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    /** Print the number of particles */
    virtual void L2P(const CellClass* const , FList<ParticleClass*>* const ){
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    /** Print the number of particles */
    virtual void P2P(FList<ParticleClass*>* const FRestrict , const FList<ParticleClass*>* const FRestrict ,
                     const FList<ParticleClass*>* const [26], const int ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    virtual void P2P(const MortonIndex ,
                     FList<ParticleClass*>* const FRestrict , const FList<ParticleClass*>* const FRestrict ,
                     FList<ParticleClass*>* const [26], const MortonIndex [26], const int ){
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }
};


#endif //FBASICKERNELS_HPP

// [--LICENSE--]
