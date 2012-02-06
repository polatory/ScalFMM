// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================
#ifndef FTESTPERIODICKERNELS_HPP
#define FTESTPERIODICKERNELS_HPP

#include "FTestKernels.hpp"

/** This kernel is the periodic version of the test kernel
  */

template< class ParticleClass, class CellClass, class ContainerClass>
class FTestPeriodicKernels : public FTestKernels<ParticleClass,CellClass,ContainerClass> {
public:

    /** Before Downward */
    void M2L(CellClass* const FRestrict pole, const CellClass* distantNeighbors[189], const FTreeCoordinate [189], const int size, const int ) {
        // The pole is impacted by what represent other poles
        for(int idx = 0 ; idx < size ; ++idx){
            pole->setDataDown(pole->getDataDown() + distantNeighbors[idx]->getDataUp());
        }
    }


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
