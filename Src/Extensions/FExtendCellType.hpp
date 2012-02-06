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
#ifndef FEXTENDCELLTYPE_HPP
#define FEXTENDCELLTYPE_HPP



/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendCellType
* Please read the license
* This class is an extenssion.
* It proposes a target/source extenssion for cell.
* Because cells may have child that contains only
* sources or targets (in Tsm system) then it is important
* to not compute for nothing.
*/
class FExtendCellType {
protected:
    /** Particle potential type */
    static const int Neither = 0;
    static const int ContainsSrc = 1;
    static const int ContainsTargets = 2;

    /** Current type */
    int type;

public:
    /** Default constructor */
    FExtendCellType() : type(Neither) {
    }

    /** Copy constructor */
    FExtendCellType(const FExtendCellType& other) : type(other.type) {
    }

    /** Copy operator */
    FExtendCellType& operator=(const FExtendCellType& other) {
        this->type = other.type;
        return *this;
    }

    /** To know if a cell has sources */
    bool hasSrcChild() const {
        return this->type & ContainsSrc;
    }

    /** To know if a cell has targets */
    bool hasTargetsChild() const {
        return this->type & ContainsTargets;
    }

    /** To set cell as sources container */
    void setSrcChildTrue() {
        this->type |= ContainsSrc;
    }

    /** To set cell as targets container */
    void setTargetsChildTrue() {
        this->type |= ContainsTargets;
    }
};


#endif //FEXTENDCELLTYPE_HPP


