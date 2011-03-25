#ifndef FFMBKERNELSPOTENTIAL_HPP
#define FFMBKERNELSPOTENTIAL_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FAbstractFmbKernels.hpp"



/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class AbstractKernels
* @brief
* Please read the license
*
* This kernels simply shows the details of the information
* it receives
*/
template< class ParticuleClass, class CellClass>
class FFmbKernelsPotential : public FAbstractFmbKernels<ParticuleClass,CellClass> {
public:
    FFmbKernelsPotential(const int inTreeHeight, const double inTreeWidth)
        : FAbstractFmbKernels<ParticuleClass,CellClass>(inTreeHeight,inTreeWidth) {
    }


    void changeProgression(int*const start_for_j , FComplexe** const p_target_exp_term){}

    /** bodies_L2P
      * expansion_L2P_add_to_force_vector
      */
    void L2P(CellClass* const local, FList<ParticuleClass*>* const particules){
        typename FList<ParticuleClass*>::BasicIterator iterTarget(*particules);
        while( iterTarget.isValide() ){
            //printf("Morton %lld\n",local->getMortonIndex());

            // expansion_Evaluate_local
            typename FAbstractFmbKernels<ParticuleClass,CellClass>::Spherical spherical;
            positionToSphere( iterTarget.value()->getPosition() - local->getPosition(), &spherical );
            harmonicInner( spherical, FAbstractFmbKernels<ParticuleClass,CellClass>::current_thread_Y);

            double potential;
            expansion_Evaluate_local_with_Y_already_computed(local->getLocal(),&potential);
            iterTarget.value()->setPotential(potential);

            //printf("\t fx = %f \t fy = %f \t fz = %f \n",iterTarget.value()->getForces().getX(),iterTarget.value()->getForces().getY(),iterTarget.value()->getForces().getZ());
            //printf("p_potential = %lf\n", potential);

            iterTarget.progress();
        }
    }


    void expansion_Evaluate_local_with_Y_already_computed(const FComplexe* local_exp,
                                                          double* const p_result){

        double result = 0.0;

        FComplexe* p_Y_term = FAbstractFmbKernels<ParticuleClass,CellClass>::current_thread_Y;
        for(int j = 0 ; j<= FAbstractFmbKernels<ParticuleClass,CellClass>::FMB_Info_P ; ++j){
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

    }




    /** void bodies_Compute_direct_interaction 	(
      *          bodies_t *FMB_RESTRICT  	p_b_target,
      *          bodies_t *FMB_RESTRICT  	p_b_src,
      *          bool  	mutual
      *  )
      *
      */
    void P2P(FList<ParticuleClass*>* const currentBox, FList<ParticuleClass*>** directNeighbors, const int size) {
        typename FList<ParticuleClass*>::BasicIterator iterTarget(*currentBox);
        while( iterTarget.isValide() ){
            for(int idxDirectNeighbors = 0 ; idxDirectNeighbors < size ; ++idxDirectNeighbors){
                typename FList<ParticuleClass*>::BasicIterator iterSource(*directNeighbors[idxDirectNeighbors]);
                while( iterSource.isValide() ){
                  DIRECT_COMPUTATION_NO_MUTUAL_SOFT(&iterTarget.value(),
                                                   iterSource.value());
                  iterSource.progress();
                }
             }

            typename FList<ParticuleClass*>::BasicIterator iterSameBox(*currentBox);
            while( iterSameBox.isValide() ){
                if(iterSameBox.value() != iterTarget.value()){
                    DIRECT_COMPUTATION_NO_MUTUAL_SOFT(&iterTarget.value(),
                                                     iterSameBox.value());
                }
                iterSameBox.progress();
            }

            //printf("They contains energy (res = %f, potential = %f, value = %f)\n",
                            //potential_sum, iterTarget.value()->getPotential(), iterTarget.value()->getValue());
            FAbstractFmbKernels<ParticuleClass,CellClass>::potential_sum += iterTarget.value()->getPotential() * iterTarget.value()->getValue();

            iterTarget.progress();
        }
    }
    void DIRECT_COMPUTATION_NO_MUTUAL_SOFT(ParticuleClass** const target, const ParticuleClass* const source){
      const double dx = (*target)->getPosition().getX() - source->getPosition().getX();
      const double dy = (*target)->getPosition().getY() - source->getPosition().getY();
      const double dz = (*target)->getPosition().getZ() - source->getPosition().getZ();

      double inv_distance = 1.0/FMath::Sqrt(dx*dx + dy*dy + dz*dz + FAbstractFmbKernels<ParticuleClass,CellClass>::FMB_Info_eps_soft_square);
      inv_distance *= (*target)->getValue() * source->getValue();

      (*target)->setPotential( inv_distance + (*target)->getPotential());
      //source->setPotential( inv_distance + source->getPotential());
    }
};


#endif //FFMBKERNELSPOTENTIAL_HPP

// [--LICENSE--]
