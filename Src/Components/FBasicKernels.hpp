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
#ifndef FBASICKERNELS_HPP
#define FBASICKERNELS_HPP


#include "FAbstractKernels.hpp"

#include "../Utils/FGlobal.hpp"
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
template< class ParticleClass, class CellClass, class ContainerClass>
class FBasicKernels : public FAbstractKernels<ParticleClass,CellClass,ContainerClass> {
public:
    /** Default destructor */
    virtual ~FBasicKernels(){
    }

    /** Do nothing */
    virtual void P2M(CellClass* const , const ContainerClass* const ) {

    }

    /** Do nothing */
    virtual void M2M(CellClass* const FRestrict , const CellClass*const FRestrict *const FRestrict , const int ) {

    }

    /** Do nothing */
    virtual void M2L(CellClass* const FRestrict , const CellClass* [], const int , const int ) {

    }

    /** Do nothing */
    virtual void L2L(const CellClass* const FRestrict , CellClass* FRestrict *const FRestrict  , const int ) {

    }

    /** Do nothing */
    virtual void L2P(const CellClass* const , ContainerClass* const ){

    }

    /** Do nothing */
    virtual void P2P(ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
                     const ContainerClass* const [26], const int ) {

    }

    /** Do nothing */
    virtual void P2P(const MortonIndex ,
                     ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
                     ContainerClass* const [26], const MortonIndex [26], const int ){

    }

    // ------------------- Periodic  --------------------


    /** Before Downward */
    void M2L(CellClass* const FRestrict , const CellClass* [189], const FTreeCoordinate [189], const int , const int ) {
    }


    /** After Downward */
    void P2P(const MortonIndex ,
             ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
             ContainerClass* const [26], const FTreeCoordinate [26], const int ) {
    }
};


#endif //FBASICKERNELS_HPP


