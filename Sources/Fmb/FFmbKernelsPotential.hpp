#ifndef FFMBKERNELSPOTENTIAL_HPP
#define FFMBKERNELSPOTENTIAL_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FGlobal.hpp"
#include "FAbstractFmbKernels.hpp"



/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class AbstractKernels
* @brief
* Please read the license
*
*
* This class is a Fmb Kernels.
*/
template< class ParticuleClass, class CellClass, int TreeHeight>
class FFmbKernelsPotential : public FAbstractFmbKernels<ParticuleClass,CellClass, TreeHeight> {
public:
    FFmbKernelsPotential(const FReal inTreeWidth)
        : FAbstractFmbKernels<ParticuleClass,CellClass,TreeHeight>(inTreeWidth) {
    }


    void changeProgression(int*const start_for_j , FComplexe** const p_target_exp_term){}

    /** bodies_L2P
      * expansion_L2P_add_to_force_vector
      */
    void L2P(const CellClass* const local, FList<ParticuleClass*>* const particules){
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        typename FList<ParticuleClass*>::BasicIterator iterTarget(*particules);
        while( iterTarget.isValide() ){
            //printf("Morton %lld\n",local->getMortonIndex());

            // expansion_Evaluate_local
            harmonicInner( positionToSphere( iterTarget.value()->getPosition() - local->getPosition()),
                           FAbstractFmbKernels<ParticuleClass,CellClass,TreeHeight>::current_thread_Y);

            FReal potential;
            expansion_Evaluate_local_with_Y_already_computed(local->getLocal(),&potential);
            iterTarget.value()->setPotential(potential);

            //printf("\t fx = %f \t fy = %f \t fz = %f \n",iterTarget.value()->getForces().getX(),iterTarget.value()->getForces().getY(),iterTarget.value()->getForces().getZ());
            //printf("p_potential = %lf\n", potential);

            iterTarget.progress();
        }
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }


    void expansion_Evaluate_local_with_Y_already_computed(const FComplexe* local_exp,
                                                          FReal* const p_result){
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );

        FReal result = 0.0;

        FComplexe* p_Y_term = FAbstractFmbKernels<ParticuleClass,CellClass,TreeHeight>::current_thread_Y;
        for(int j = 0 ; j<= FMB_Info_P ; ++j){
            // k=0
            (*p_Y_term) *= (*local_exp);
            result += p_Y_term->getReal();
            //printf("\t\t p_Y_term->real = %f p_Y_term->imag = %f \t local_exp->real = %f local_exp->imag = %f \n",
            //       p_Y_term->getReal(), p_Y_term->getImag(), local_exp->getReal(), local_exp->getImag());
            ++p_Y_term;
            ++local_exp;

            // k>0
            for (int k=1; k<=j ;++k, ++p_Y_term, ++local_exp){
                (*p_Y_term) *= (*local_exp);
                result += 2 * p_Y_term->getReal();
                //printf("\t\t p_Y_term->real = %f p_Y_term->imag = %f \t local_exp->real = %f local_exp->imag = %f \n",
                //       p_Y_term->getReal(), p_Y_term->getImag(), local_exp->getReal(), local_exp->getImag());
            }
        }

        *p_result = result;

        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }




    /** void bodies_Compute_direct_interaction 	(
      *          bodies_t *FMB_RESTRICT  	p_b_target,
      *          bodies_t *FMB_RESTRICT  	p_b_src,
      *          bool  	mutual
      *  )
      *
      */
    void P2P(FList<ParticuleClass*>* const FRestrict targets, const FList<ParticuleClass*>* const FRestrict sources,
             const FList<ParticuleClass*>* FRestrict const* FRestrict directNeighbors, const int size) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        typename FList<ParticuleClass*>::BasicIterator iterTarget(*targets);
        while( iterTarget.isValide() ){
            for(int idxDirectNeighbors = 0 ; idxDirectNeighbors < size ; ++idxDirectNeighbors){
                typename FList<ParticuleClass*>::ConstBasicIterator iterSource(*directNeighbors[idxDirectNeighbors]);
                while( iterSource.isValide() ){
                  DIRECT_COMPUTATION_NO_MUTUAL_SOFT(&iterTarget.value(), iterSource.value());
                  iterSource.progress();
                }
             }

            typename FList<ParticuleClass*>::ConstBasicIterator iterSameBox(*sources);
            while( iterSameBox.isValide() ){
                if(iterSameBox.value() != iterTarget.value()){
                    DIRECT_COMPUTATION_NO_MUTUAL_SOFT(&iterTarget.value(), iterSameBox.value());
                }
                iterSameBox.progress();
            }


            iterTarget.progress();
        }
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }
    void DIRECT_COMPUTATION_NO_MUTUAL_SOFT(ParticuleClass** const target, const ParticuleClass* const source){
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
      const FReal dx = (*target)->getPosition().getX() - source->getPosition().getX();
      const FReal dy = (*target)->getPosition().getY() - source->getPosition().getY();
      const FReal dz = (*target)->getPosition().getZ() - source->getPosition().getZ();

      FReal inv_distance = 1.0/FMath::Sqrt(dx*dx + dy*dy + dz*dz + FAbstractFmbKernels<ParticuleClass,CellClass,TreeHeight>::FMB_Info_eps_soft_square);
      inv_distance *= (*target)->getPhysicalValue() * source->getPhysicalValue();

      (*target)->setPotential( inv_distance + (*target)->getPotential());
      //source->setPotential( inv_distance + source->getPotential());
      FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }
};


#endif //FFMBKERNELSPOTENTIAL_HPP

// [--LICENSE--]
