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
#ifndef FSPHERICALCELL_HPP
#define FSPHERICALCELL_HPP


#include "../Utils/FComplexe.hpp"
#include "../Utils/FMemUtils.hpp"

#include "../Components/FBasicCell.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FSphericalCell
* Please read the license.
*
*/
class FSphericalCell : public FBasicCell {
protected:
    static int DevP;
    static int ExpP;

    FComplexe* multipole_exp; //< For multipole extenssion
    FComplexe* local_exp;     //< For local extenssion

public:
    static void Init(const int inDevP){
        DevP = inDevP;
        ExpP = int((DevP+1) * (DevP+2) * 0.5);
    }

    static int GetP(){
        return DevP;
    }

    static int GetExp(){
        return ExpP;
    }


    /** Default constructor */
    FSphericalCell()
        : multipole_exp(0), local_exp(0){
        multipole_exp = new FComplexe[ExpP];
        local_exp = new FComplexe[ExpP];
    }

    /** Constructor */
    FSphericalCell(const FSphericalCell& other)
        : multipole_exp(0), local_exp(0){
        multipole_exp = new FComplexe[ExpP];
        local_exp = new FComplexe[ExpP];
        (*this) = other;
    }

    /** Default destructor */
    virtual ~FSphericalCell(){
        delete[] multipole_exp;
        delete[] local_exp;
    }

    /** Copy constructor */
    FSphericalCell& operator=(const FSphericalCell& other) {
        FMemUtils::copyall(multipole_exp, other.multipole_exp, ExpP);
        FMemUtils::copyall(local_exp, other.local_exp, ExpP);
        return *this;
    }

    /** Get Multipole */
    const FComplexe* getMultipole() const {
        return multipole_exp;
    }
    /** Get Local */
    const FComplexe* getLocal() const {
        return local_exp;
    }

    /** Get Multipole */
    FComplexe* getMultipole() {
        return multipole_exp;
    }
    /** Get Local */
    FComplexe* getLocal() {
        return local_exp;
    }
};


int FSphericalCell::DevP(-1);
int FSphericalCell::ExpP(-1);


#endif //FSPHERICALCELL_HPP


