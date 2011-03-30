#ifndef FFMBKERNELSPOTENTIALFORCES_HPP
#define FFMBKERNELSPOTENTIALFORCES_HPP
// /!\ Please, you must read the license at the bottom of this page


#include "FAbstractFmbKernels.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFmbKernelsPotentialForces
* @brief
* Please read the license
*
* This class is a Fmb Kernels.
*/
template< class ParticuleClass, class CellClass>
class FFmbKernelsPotentialForces : public FAbstractFmbKernels<ParticuleClass,CellClass> {

public:
    FFmbKernelsPotentialForces(const int inTreeHeight, const FReal inTreeWidth)
        : FAbstractFmbKernels<ParticuleClass,CellClass>(inTreeHeight,inTreeWidth) {
    }   


    void changeProgression(int*const , FComplexe** const ){}

    /** bodies_L2P
      *     expansion_L2P_add_to_force_vector_and_to_potential
      *         expansion_L2P_add_to_force_vector
      *         expansion_Evaluate_local_with_Y_already_computed
      */
    void L2P(const CellClass* const local, FList<ParticuleClass*>* const particules){
        typename FList<ParticuleClass*>::BasicIterator iterTarget(*particules);
        while( iterTarget.isValide() ){
            //printf("Morton %lld\n",local->getMortonIndex());

            F3DPosition force_vector_in_local_base;
            typename FAbstractFmbKernels<ParticuleClass,CellClass>::Spherical spherical;

            positionToSphere( iterTarget.value()->getPosition() - local->getPosition(), &spherical );

            /*printf("\t\t bodies_it_Get_p_position(&it) x = %lf \t y = %lf \t z = %lf \n",
                   (iterTarget.value()->getPosition()).getX(),
                   (iterTarget.value()->getPosition()).getY(),
                   (iterTarget.value()->getPosition()).getZ());
            printf("\t\t p_leaf_center x = %lf \t y = %lf \t z = %lf \n",
                   (local->getPosition()).getX(),
                   (local->getPosition()).getY(),
                   (local->getPosition()).getZ());*/
            /*printf("\t\t p_position_to_leaf_center x = %lf \t y = %lf \t z = %lf \n",
                    (iterTarget.value()->getPosition() - local->getPosition()).getX(),
                    (iterTarget.value()->getPosition() - local->getPosition()).getY(),
                    (iterTarget.value()->getPosition() - local->getPosition()).getZ());*/
            /*printf("\t\t phi = %lf \t cos = %lf \t sin = %lf \t r= %lf \n",
                    spherical.phi,spherical.cosTheta,spherical.sinTheta,spherical.r);*/

            harmonicInnerThetaDerivated( spherical, FAbstractFmbKernels<ParticuleClass,CellClass>::current_thread_Y, FAbstractFmbKernels<ParticuleClass,CellClass>::current_thread_Y_theta_derivated);

            // The maximum degree used here will be P.
            const FComplexe* p_Y_term = FAbstractFmbKernels<ParticuleClass,CellClass>::current_thread_Y+1;
            const FComplexe* p_Y_theta_derivated_term = FAbstractFmbKernels<ParticuleClass,CellClass>::current_thread_Y_theta_derivated+1;
            const FComplexe* p_local_exp_term = local->getLocal()+1;

            for (int j = 1 ; j <= FAbstractFmbKernels<ParticuleClass,CellClass>::FMB_Info_P ; ++j ){
                FComplexe exp_term_aux;

                // k=0:
                // F_r:
                exp_term_aux.setReal( (p_Y_term->getReal() * p_local_exp_term->getReal()) - (p_Y_term->getImag() * p_local_exp_term->getImag()) );
                exp_term_aux.setImag( (p_Y_term->getReal() * p_local_exp_term->getImag()) + (p_Y_term->getImag() * p_local_exp_term->getReal()) );

                force_vector_in_local_base.setX( force_vector_in_local_base.getX() + j * exp_term_aux.getReal());
                // F_phi: k=0 => nothing to do for F_phi
                // F_theta:
                exp_term_aux.setReal( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getReal()) - (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getImag()) );
                exp_term_aux.setImag( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getImag()) + (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getReal()) );

                force_vector_in_local_base.setY( force_vector_in_local_base.getY() + exp_term_aux.getReal());


                /*printf("\t j = %d \t exp_term_aux real = %lf imag = %lf \n",
                                                        j,
                                                        exp_term_aux.getReal(),
                                                        exp_term_aux.getImag());*/
                /*printf("\t\t\t p_Y_theta_derivated_term->getReal = %lf \t p_Y_theta_derivated_term->getImag = %lf \n",
                       p_Y_theta_derivated_term->getReal(),p_Y_theta_derivated_term->getImag());*/
                //printf("\t\t\t p_local_exp_term->getReal = %lf \t p_local_exp_term->getImag = %lf \n",
                //       p_local_exp_term->getReal(),p_local_exp_term->getImag());
                //printf("\t\t\t p_Y_term->getReal = %lf \t p_Y_term->getImag = %lf \n",p_Y_term->getReal(),p_Y_term->getImag());

                //printf("[LOOP1] force_vector_in_local_base x = %lf \t y = %lf \t z = %lf \n",
                //       force_vector_in_local_base.getX(),force_vector_in_local_base.getY(),force_vector_in_local_base.getZ());


                ++p_local_exp_term;
                ++p_Y_term;
                ++p_Y_theta_derivated_term;


                // k>0:
                for (int k=1; k<=j ;++k, ++p_local_exp_term, ++p_Y_term, ++p_Y_theta_derivated_term){
                    // F_r:

                    exp_term_aux.setReal( (p_Y_term->getReal() * p_local_exp_term->getReal()) - (p_Y_term->getImag() * p_local_exp_term->getImag()) );
                    exp_term_aux.setImag( (p_Y_term->getReal() * p_local_exp_term->getImag()) + (p_Y_term->getImag() * p_local_exp_term->getReal()) );

                    force_vector_in_local_base.setX(force_vector_in_local_base.getX() + 2 * j * exp_term_aux.getReal());
                    // F_phi:
                    force_vector_in_local_base.setZ( force_vector_in_local_base.getZ() - 2 * k * exp_term_aux.getImag());
                    // F_theta:

                    /*printf("\t\t k = %d \t j = %d \t exp_term_aux real = %e imag = %e \n",
                                                            k,j,
                                                            exp_term_aux.getReal(),
                                                            exp_term_aux.getImag());*/
                    //printf("\t\t\t p_Y_term->getReal = %e \t p_Y_term->getImag = %e \n",p_Y_term->getReal(),p_Y_term->getImag());
                    //printf("\t\t\t p_local_exp_term->getReal = %e \t p_local_exp_term->getImag = %e \n",p_local_exp_term->getReal(),p_local_exp_term->getImag());

                    exp_term_aux.setReal( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getReal()) - (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getImag()) );
                    exp_term_aux.setImag( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getImag()) + (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getReal()) );

                    force_vector_in_local_base.setY(force_vector_in_local_base.getY() + 2 * exp_term_aux.getReal());

                    /*printf("\t\t k = %d \t j = %d \t exp_term_aux real = %lf imag = %lf \n",
                                                            k,j,
                                                            exp_term_aux.getReal(),
                                                            exp_term_aux.getImag());*/
                    /*printf("\t\t\t p_Y_theta_derivated_term->getReal = %lf \t p_Y_theta_derivated_term->getImag = %lf \n",
                           p_Y_theta_derivated_term->getReal(),p_Y_theta_derivated_term->getImag());*/
                    /*printf("\t\t\t p_local_exp_term->getReal = %lf \t p_local_exp_term->getImag = %lf \n",
                           p_local_exp_term->getReal(),p_local_exp_term->getImag());*/

                    /*printf("[LOOP2] force_vector_in_local_base x = %lf \t y = %lf \t z = %lf \n",
                           force_vector_in_local_base.getX(),force_vector_in_local_base.getY(),force_vector_in_local_base.getZ());*/
                }
            }

            /*printf("[END LOOP] force_vector_in_local_base x = %lf \t y = %lf \t z = %lf \n",
                   force_vector_in_local_base.getX(),force_vector_in_local_base.getY(),force_vector_in_local_base.getZ());*/

            // We want: - gradient(POTENTIAL_SIGN potential).
            // The -(- 1.0) computing is not the most efficient programming ...
            //#define FMB_TMP_SIGN -(POTENTIAL_SIGN 1.0)
            force_vector_in_local_base.setX( force_vector_in_local_base.getX() * (-1.0) / spherical.r);
            force_vector_in_local_base.setY( force_vector_in_local_base.getY() * (-1.0) / spherical.r);
            force_vector_in_local_base.setZ( force_vector_in_local_base.getZ() * (-1.0) / (spherical.r * spherical.sinTheta));
            //#undef FMB_TMP_SIGN

            /////////////////////////////////////////////////////////////////////

            //spherical_position_Set_ph
            //FMB_INLINE COORDINATES_T angle_Convert_in_MinusPi_Pi(COORDINATES_T a){
            FReal ph = FMath::Fmod(spherical.phi, 2*FMath::FPi);
            if (ph > M_PI) ph -= 2*FMath::FPi;
            if (ph < -M_PI + FMath::Epsilon)  ph += 2 * FMath::Epsilon;

            //spherical_position_Set_th
            FReal th = FMath::Fmod(FMath::ACos(spherical.cosTheta), 2*FMath::FPi);
            if (th < 0.0) th += 2*FMath::FPi;
            if (th > FMath::FPi){
                th = 2*FMath::FPi - th;
                //spherical_position_Set_ph(p, spherical_position_Get_ph(p) + M_PI);
                    ph = FMath::Fmod(ph + FMath::FPi, 2*FMath::FPi);
                    if (ph > M_PI) ph -= 2*FMath::FPi;
                    if (ph < -M_PI + FMath::Epsilon)  ph += 2 * FMath::Epsilon;
                th = FMath::Fmod(th, 2*FMath::FPi);
                if (th > M_PI) th -= 2*FMath::FPi;
                if (th < -M_PI + FMath::Epsilon)  th += 2 * FMath::Epsilon;
            }
            //spherical_position_Set_r
            FReal rh = spherical.r;
            if (spherical.r < 0){
                rh = -spherical.r;
                //spherical_position_Set_ph(p, M_PI - spherical_position_Get_th(p));
                ph = FMath::Fmod(FMath::FPi - th, 2*FMath::FPi);
                if (ph > M_PI) ph -= 2*FMath::FPi;
                if (ph < -M_PI + FMath::Epsilon)  ph += 2 * FMath::Epsilon;
                //spherical_position_Set_th(p, spherical_position_Get_th(p) + M_PI);
                th = FMath::Fmod(th + FMath::FPi, 2*FMath::FPi);
                if (th < 0.0) th += 2*FMath::FPi;
                if (th > FMath::FPi){
                    th = 2*FMath::FPi - th;
                    //spherical_position_Set_ph(p, spherical_position_Get_ph(p) + M_PI);
                        ph = FMath::Fmod(ph + FMath::FPi, 2*FMath::FPi);
                        if (ph > M_PI) ph -= 2*FMath::FPi;
                        if (ph < -M_PI + FMath::Epsilon)  ph += 2 * FMath::Epsilon;
                    th = FMath::Fmod(th, 2*FMath::FPi);
                    if (th > M_PI) th -= 2*FMath::FPi;
                    if (th < -M_PI + FMath::Epsilon)  th += 2 * FMath::Epsilon;
                }
            }

            /*printf("[details] ph = %f , rh = %f , th = %f \n",
                   ph,rh,th);*/


            const FReal cos_theta = FMath::Cos(th);
            const FReal cos_phi = FMath::Cos(ph);
            const FReal sin_theta = FMath::Sin(th);
            const FReal sin_phi = FMath::Sin(ph);

            /*printf("[details] cos_theta = %f \t cos_phi = %f \t sin_theta = %f \t sin_phi = %f \n",
                   cos_theta, cos_phi, sin_theta, sin_phi);*/
            /*printf("[force_vector_in_local_base] x = %lf \t y = %lf \t z = %lf \n",
                   force_vector_in_local_base.getX(),force_vector_in_local_base.getY(),force_vector_in_local_base.getZ());*/

            F3DPosition force_vector_tmp;

            force_vector_tmp.setX(
                    cos_phi * sin_theta * force_vector_in_local_base.getX() +
                    cos_phi * cos_theta * force_vector_in_local_base.getY() +
                    (-sin_phi) * force_vector_in_local_base.getZ());

            force_vector_tmp.setY(
                    sin_phi * sin_theta * force_vector_in_local_base.getX() +
                    sin_phi * cos_theta * force_vector_in_local_base.getY() +
                    cos_phi * force_vector_in_local_base.getZ());

            force_vector_tmp.setZ(
                    cos_theta * force_vector_in_local_base.getX() +
                    (-sin_theta) * force_vector_in_local_base.getY());

            /*printf("[force_vector_tmp]  = %lf \t y = %lf \t z = %lf \n",
                   force_vector_tmp.getX(),force_vector_tmp.getY(),force_vector_tmp.getZ());*/

            /////////////////////////////////////////////////////////////////////

            //#ifndef _DIRECT_MATRIX_
            // When _DIRECT_MATRIX_ is defined, this multiplication is done in 'leaf_Sum_near_and_far_fields()'
            force_vector_tmp *= iterTarget.value()->getValue();
            //#endif

            /*printf("[force_vector_tmp] fx = %f \t fy = %f \t fz = %f \n",
                   force_vector_tmp.getX(),force_vector_tmp.getY(),force_vector_tmp.getZ());*/

            iterTarget.value()->setForces( iterTarget.value()->getForces() + force_vector_tmp );

            FReal potential;
            expansion_Evaluate_local_with_Y_already_computed(local->getLocal(),&potential);
            iterTarget.value()->setPotential(potential);

            /*printf("[END] fx = %f \t fy = %f \t fz = %f \n\n",
                   iterTarget.value()->getForces().getX(),iterTarget.value()->getForces().getY(),iterTarget.value()->getForces().getZ());*/
            //printf("p_potential = %lf\n", potential);

            iterTarget.progress();
        }
    }


    void expansion_Evaluate_local_with_Y_already_computed(const FComplexe* local_exp,
                                                          FReal* const p_result){

        FReal result = 0.0;

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
    void P2P(FList<ParticuleClass*>* const FRestrict currentBox, const FList<ParticuleClass*>* FRestrict const* FRestrict directNeighbors, const int size) {
        typename FList<ParticuleClass*>::BasicIterator iterTarget(*currentBox);
        while( iterTarget.isValide() ){

            for(int idxDirectNeighbors = 0 ; idxDirectNeighbors < size ; ++idxDirectNeighbors){
                typename FList<ParticuleClass*>::ConstBasicIterator iterSource(*directNeighbors[idxDirectNeighbors]);
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

            //printf("x = %f \t y = %f \t z = %f \n",iterTarget.value()->getPosition().getX(),iterTarget.value()->getPosition().getY(),iterTarget.value()->getPosition().getZ());
            //printf("\t P2P fx = %f \t fy = %f \t fz = %f \n",iterTarget.value()->getForces().getX(),iterTarget.value()->getForces().getY(),iterTarget.value()->getForces().getZ());
            //printf("\t potential = %f \n",iterTarget.value()->getPotential());

            FAbstractFmbKernels<ParticuleClass,CellClass>::force_sum += iterTarget.value()->getForces();

            FAbstractFmbKernels<ParticuleClass,CellClass>::potential_sum += iterTarget.value()->getPotential() * iterTarget.value()->getValue();
            //printf("They contain energy (res = %f, potential = %f, value = %f)\n",
            //                FAbstractFmbKernels<ParticuleClass,CellClass>::potential_sum, iterTarget.value()->getPotential(), iterTarget.value()->getValue());

            iterTarget.progress();
        }
    }
    void DIRECT_COMPUTATION_NO_MUTUAL_SOFT(ParticuleClass** const target, const ParticuleClass* const source){
      const FReal dx = (*target)->getPosition().getX() - source->getPosition().getX();
      const FReal dy = (*target)->getPosition().getY() - source->getPosition().getY();
      const FReal dz = (*target)->getPosition().getZ() - source->getPosition().getZ();

      FReal inv_square_distance = 1.0/ (dx*dx + dy*dy + dz*dz + FAbstractFmbKernels<ParticuleClass,CellClass>::FMB_Info_eps_soft_square);
      FReal inv_distance = FMath::Sqrt(inv_square_distance);
      inv_distance *= (*target)->getValue() * source->getValue();
      inv_square_distance *= inv_distance;

      (*target)->setForces(
              (*target)->getForces().getX() + dx * inv_square_distance,
              (*target)->getForces().getY() + dy * inv_square_distance,
              (*target)->getForces().getZ() + dz * inv_square_distance
      );

      (*target)->setPotential( inv_distance + (*target)->getPotential());
    }
};


#endif //FFMBKERNELSPOTENTIALFORCES_HPP

// [--LICENSE--]
