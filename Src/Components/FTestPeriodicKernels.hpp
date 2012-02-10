// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FTESTPERIODICKERNELS_HPP
#define FTESTPERIODICKERNELS_HPP

#include "FTestKernels.hpp"

/** This kernel is the periodic version of the test kernel
  */

template< class ParticleClass, class CellClass, class ContainerClass>
class FTestPeriodicKernels : public FTestKernels<ParticleClass,CellClass,ContainerClass> {
public:


    /** After Downward */
    void P2P(const MortonIndex ,
             ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
             ContainerClass* const directNeighborsParticles[26], const FTreeCoordinate [26], const int size) {

        // Each particles targeted is impacted by the particles sources
        long long int inc = sources->getSize();
        if(targets == sources){
            inc -= 1;
        }
        for(int idx = 0 ; idx < size ; ++idx){
            inc += directNeighborsParticles[idx]->getSize();
        }

        typename ContainerClass::BasicIterator iter(*targets);
        while( iter.hasNotFinished() ){
            iter.data().setDataDown(iter.data().getDataDown() + inc);
            iter.gotoNext();
        }

    }
};

#endif // FTESTPERIODICKERNELS_HPP
