#ifndef FELECBASICKERNELS_HPP
#define FELECBASICKERNELS_HPP
// [--License--]

#include "../Components/FAbstractKernels.hpp"

#include "../Utils/FGlobal.hpp"
#include "../Utils/FTrace.hpp"

#include "FHarmonic.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @brief
* Please read the license
*
*/
template< class ParticleClass, class CellClass, class ContainerClass>
class FElecForcesKernels : public FAbstractKernels<ParticleClass,CellClass,ContainerClass> {
    int devP;
    FHarmonic harmonic;

public:
    FElecForcesKernels(const int inDevP)
        : devP(inDevP), harmonic(inDevP) {
    }

    /** Default destructor */
    ~FElecForcesKernels(){
    }

    /** P2M with a cell and all its particles */
    void P2M(CellClass* const inPole, const ContainerClass* const ) {
        FComplexe* const cellMultiPole = inPole->getMultipole();
        // Copying the position is faster than using cell position
        const F3DPosition polePosition = inPole->getPosition();
        // For all particles in the leaf box
        typename ContainerClass::ConstBasicIterator iterParticle(*inParticles);
        while( iterParticle.hasNotFinished()){
            // P2M
            particleToMultiPole(cellMultiPole, polePosition, iterParticle.data());
            iterParticle.gotoNext();
        }
    }

    /** Do nothing */
    void M2M(CellClass* const FRestrict , const CellClass*const FRestrict *const FRestrict , const int ) {

    }

    /** Do nothing */
    virtual void M2L(CellClass* const FRestrict , const CellClass* [], const int , const int ) {

    }

    /** Do nothing */
    void L2L(const CellClass* const FRestrict , CellClass* FRestrict *const FRestrict  , const int ) {

    }

    /** L2P with a cell and all its particles */
    void L2P(const CellClass* const local, ContainerClass* const particles){
        const FComplexe* const cellLocal = local->getLocal();
        // Copying the position is faster than using cell position
        const F3DPosition localPosition = local->getPosition();
        // For all particles in the leaf box
        typename ContainerClass::BasicIterator iterTarget(*particles);
        while( iterTarget.hasNotFinished() ){
            // L2P
            localToParticle(&iterTarget.data(), localPosition, cellLocal);
            iterTarget.gotoNext();
        }
    }

    /** This P2P has to be used when target != sources
      * It will proceed an direct interation no mutual
      *
      * It takes all the target particles from the current leaf,
      * then it computes the sources/targets interaction in this leaf,
      * then it computes the sources/targets inteactions between this leaf and the
      * neighbors.
      */
    void P2P(ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
             const ContainerClass* const directNeighbors[26], const int size) {

        { // Compute interaction in this leaf
            typename ContainerClass::BasicIterator iterTarget(*targets);
            while( iterTarget.hasNotFinished() ){
                // We copy the target particle to work with a particle in the heap
                ParticleClass target( iterTarget.data() );

                // For all the source particles in the same leaf
                typename ContainerClass::ConstBasicIterator iterSameBox(*sources);
                while( iterSameBox.hasNotFinished() ){
                    //(&iterSameBox.data() != &iterTarget.data())
                    directInteraction(target, iterSameBox.data());
                    iterSameBox.gotoNext();
                }
                // Set data and progress
                iterTarget.setData(target);
                iterTarget.gotoNext();
            }
        }
        { // Compute interactions with other leaves
            // For all the neigbors leaves
            for(int idxDirectNeighbors = 0 ; idxDirectNeighbors < size ; ++idxDirectNeighbors){
                // For all particles in current leaf
                typename ContainerClass::BasicIterator iterTarget(*targets);
                while( iterTarget.hasNotFinished() ){
                    // For all the particles in the other leaf
                    typename ContainerClass::ConstBasicIterator iterSource(*directNeighbors[idxDirectNeighbors]);
                    while( iterSource.hasNotFinished() ){
                        directInteraction(target, iterSource.data());
                        iterSource.gotoNext();
                    }
                    // Set data and progress
                    iterTarget.setData(target);
                    iterTarget.gotoNext();
                }
            }
        }
    }

    /** This P2P has to be used when target == sources
      * It will proceed a direct interation >> mutual
      *
      * It takes all the particles from the current leaf,
      * then it computes the interactions in this leaf,
      * then it computes the  inteactions between this leaf and the
      * neighbors.
      */
    void P2P(const MortonIndex inCurrentIndex,
             ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict /*sources*/,
             ContainerClass* const directNeighbors[26], const MortonIndex inNeighborsIndex[26], const int size){
        { // Compute interaction in this leaf
            typename ContainerClass::BasicIterator iterTarget(*targets);
            while( iterTarget.hasNotFinished() ){
                // We copy the target particle to work with a particle in the heap
                ParticleClass target( iterTarget.data() );

                // For all particles after the current one
                typename ContainerClass::ConstBasicIterator iterSameBox = iterTarget;
                iterSameBox.gotoNext();
                while( iterSameBox.hasNotFinished() ){
                    directInteractionMutual(&target, &iterSameBox.data());
                    iterSameBox.gotoNext();
                }
                // Set data and progress
                iterTarget.setData(target);
                iterTarget.gotoNext();
            }
        }
        { // Compute interactions with other leaves
            // For all the neigbors leaves
            for(int idxDirectNeighbors = 0 ; idxDirectNeighbors < size ; ++idxDirectNeighbors){
                if(inCurrentIndex < inNeighborsIndex[idxDirectNeighbors] ){
                    // For all particles in current leaf
                    typename ContainerClass::BasicIterator iterTarget(*targets);
                    while( iterTarget.hasNotFinished() ){
                        ParticleClass target( iterTarget.data() );
                        // For all the particles in the other leaf
                        typename ContainerClass::ConstBasicIterator iterSource(*directNeighbors[idxDirectNeighbors]);
                        while( iterSource.hasNotFinished() ){
                            directInteractionMutual(&target, &iterSource.data());
                            iterSource.gotoNext();
                        }
                        // Set data and progress
                        iterTarget.setData(target);
                        iterTarget.gotoNext();
                    }
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    //                                  Periodic
    ///////////////////////////////////////////////////////////////////////////////

    /** Before Downward */
    void M2L(CellClass* const FRestrict , const CellClass* [189], FTreeCoordinate [189], const int , const int ) {
    }


    /** After Downward */
    void P2P(const MortonIndex ,
             ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
             ContainerClass* const [26], const FTreeCoordinate [26], const int ) {
    }

    ///////////////////////////////////////////////////////////////////////////////
    //                                  Computation
    ///////////////////////////////////////////////////////////////////////////////


    /** P2M computation
      */
    void particleToMultiPole(FComplexe*const cellMultiPole, const F3DPosition& inPolePosition ,
                             const CellClass& particle){
        harmonic.computeInner( FSpherical(particle.getPosition() - inPolePosition) );

        FReal minus_one_pow_j = 1.0;//(-1)^j
        const FReal valueParticle = particle.getPhysicalValue();
        int p_exp_term = 0; // p_Y_term

        for(int jIdx = 0 ; jIdx <= devP ; ++jIdx){
            for(int kIdx = 0 ; kIdx <= jIdx ; ++kIdx, ++p_Y_term, ++p_exp_term){
                harmonic.result(p_exp_term).mulRealAndImag( valueParticle * minus_one_pow_j );
                cellMultiPole[p_exp_term] += harmonic.result(p_exp_term);
            }

            minus_one_pow_j = -minus_one_pow_j;
        }
    }

    void localToParticle(ParticleClass*const particle, const F3DPosition& local_position,
                         const FComplexe*const local_exp){
        FReal force_vector_in_local_base_x = 0;
        FReal force_vector_in_local_base_y = 0;
        FReal force_vector_in_local_base_z = 0;

        const FSpherical spherical(particle->getPosition() - local_position);
        harmonic.computeInnerTheta( spherical );

        // The maximum degree used here will be P.
        const FComplexe* p_Y_term = current_thread_Y+1;
        const FComplexe* p_Y_theta_derivated_term = current_thread_Y_theta_derivated+1;
        const FComplexe* p_local_exp_term = local_exp+1;

        for (int j = 1 ; j <= devP ; ++j ){
            FReal exp_term_aux_real = 0.0;
            FReal exp_term_aux_imag = 0.0;

            // k=0:
            // F_r:
            exp_term_aux_real = ( (p_Y_term->getReal() * p_local_exp_term->getReal()) - (p_Y_term->getImag() * p_local_exp_term->getImag()) );
            exp_term_aux_imag = ( (p_Y_term->getReal() * p_local_exp_term->getImag()) + (p_Y_term->getImag() * p_local_exp_term->getReal()) );

            force_vector_in_local_base_x = ( force_vector_in_local_base_x  + FReal(j) * exp_term_aux_real );
            // F_phi: k=0 => nothing to do for F_phi
            // F_theta:
            exp_term_aux_real = ( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getReal()) - (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getImag()) );
            exp_term_aux_imag = ( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getImag()) + (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getReal()) );

            force_vector_in_local_base_y = ( force_vector_in_local_base_y + exp_term_aux_real );

            ++p_local_exp_term;
            ++p_Y_term;
            ++p_Y_theta_derivated_term;


            // k>0:
            for (int k=1; k<=j ;++k, ++p_local_exp_term, ++p_Y_term, ++p_Y_theta_derivated_term){
                // F_r:

                exp_term_aux_real = ( (p_Y_term->getReal() * p_local_exp_term->getReal()) - (p_Y_term->getImag() * p_local_exp_term->getImag()) );
                exp_term_aux_imag = ( (p_Y_term->getReal() * p_local_exp_term->getImag()) + (p_Y_term->getImag() * p_local_exp_term->getReal()) );

                force_vector_in_local_base_x = (force_vector_in_local_base_x  + FReal(2 * j) * exp_term_aux_real );
                // F_phi:
                force_vector_in_local_base_z = ( force_vector_in_local_base_z - FReal(2 * k) * exp_term_aux_imag);
                // F_theta:

                exp_term_aux_real = ( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getReal()) - (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getImag()) );
                exp_term_aux_imag = ( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getImag()) + (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getReal()) );

                force_vector_in_local_base_y = (force_vector_in_local_base_y + FReal(2.0) * exp_term_aux_real );

           }

        }
        // We want: - gradient(POTENTIAL_SIGN potential).
        // The -(- 1.0) computing is not the most efficient programming ...
        force_vector_in_local_base_x = ( force_vector_in_local_base_x  * FReal(-1.0) / spherical.getR());
        force_vector_in_local_base_y = ( force_vector_in_local_base_y * FReal(-1.0) / spherical.getR());
        force_vector_in_local_base_z = ( force_vector_in_local_base_z * FReal(-1.0) / (spherical.getR() * spherical.getSinTheta()));

        /////////////////////////////////////////////////////////////////////

        //spherical_position_Set_ph
        //FMB_INLINE COORDINATES_T angle_Convert_in_MinusPi_Pi(COORDINATES_T a){
        FReal ph = FMath::Fmod(spherical.getPhi(), FReal(2)*FMath::FPi);
        if (ph > M_PI) ph -= FReal(2) * FMath::FPi;
        if (ph < -M_PI + FMath::Epsilon)  ph += FReal(2) * FMath::Epsilon;

        //spherical_position_Set_th
        FReal th = FMath::Fmod(FMath::ACos(spherical.getCosTheta()), FReal(2) * FMath::FPi);
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
        //FReal rh = spherical.r;
        if (spherical.r < 0){
            //rh = -spherical.r;
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

        const FReal cos_theta   = FMath::Cos(th);
        const FReal cos_phi     = FMath::Cos(ph);
        const FReal sin_theta   = FMath::Sin(th);
        const FReal sin_phi     = FMath::Sin(ph);

        FReal force_vector_tmp_x = (
                cos_phi * sin_theta * force_vector_in_local_base_x  +
                cos_phi * cos_theta * force_vector_in_local_base_y +
                (-sin_phi) * force_vector_in_local_base_z);

        FReal force_vector_tmp_y = (
                sin_phi * sin_theta * force_vector_in_local_base_x  +
                sin_phi * cos_theta * force_vector_in_local_base_y +
                cos_phi * force_vector_in_local_base_z);

        FReal force_vector_tmp_z = (
                cos_theta * force_vector_in_local_base_x +
                (-sin_theta) * force_vector_in_local_base_y);

        const FReal physicalValue = particle->getPhysicalValue();
        force_vector_tmp_x *= physicalValue;
        force_vector_tmp_y *= physicalValue;
        force_vector_tmp_z *= physicalValue;

        particle->incForces( force_vector_tmp_x, force_vector_tmp_y, force_vector_tmp_z );

        particle->incPotential(expansion_Evaluate_local_with_Y_already_computed(local_exp));
    }


    /** P2P mutual interaction
      * F = q * q' / r²
      */
    void directInteractionMutual(ParticleClass*const FRestrict target, ParticleClass*const FRestrict source){

        FReal dx = -(target->getPosition().getX() - source->getPosition().getX());
        FReal dy = -(target->getPosition().getY() - source->getPosition().getY());
        FReal dz = -(target->getPosition().getZ() - source->getPosition().getZ());

        FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
        FReal inv_distance = FMath::Sqrt(inv_square_distance);

        inv_square_distance *= inv_distance;
        inv_square_distance *= target->getPhysicalValue() * source->getPhysicalValue();

        dx *= inv_square_distance;
        dy *= inv_square_distance;
        dz *= inv_square_distance;

        target->incForces( dx, dy, dz);
        target->incPotential( inv_distance  * source->getPhysicalValue() );

        source.incForces( (-dx), (-dy), (-dz));
        source.incPotential( inv_distance * target->getPhysicalValue() );
    }

    /** P2P NO mutual interaction
      * F = q * q' / r²
      */
    void directInteraction(ParticleClass*const FRestrict target, const ParticleClass& source){

        FReal dx = -(target->getPosition().getX() - source.getPosition().getX());
        FReal dy = -(target->getPosition().getY() - source.getPosition().getY());
        FReal dz = -(target->getPosition().getZ() - source.getPosition().getZ());

        FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
        FReal inv_distance = FMath::Sqrt(inv_square_distance);

        inv_square_distance *= inv_distance;
        inv_square_distance *= target->getPhysicalValue() * source.getPhysicalValue();

        dx *= inv_square_distance;
        dy *= inv_square_distance;
        dz *= inv_square_distance;

        target->incForces( dx, dy, dz);
        target->incPotential( inv_distance  * source.getPhysicalValue() );
    }
};


#endif //FELECBASICKERNELS_HPP

// [--END--]
