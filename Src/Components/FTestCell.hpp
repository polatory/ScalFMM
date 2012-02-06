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
#ifndef FTESTCELL_HPP
#define FTESTCELL_HPP


#include "FBasicCell.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicCell
* Please read the license
*
* This class is used in the FTestKernels, please
* look at this class to know what it is.
*
* This cell simply store the data when up/down
*/
class FTestCell : public FBasicCell  {
protected:
    // To store data during upward and downward pass
    long long int dataUp, dataDown;
public:
    FTestCell(): dataUp(0) , dataDown(0){
    }
    /** Default destructor */
    virtual ~FTestCell(){
    }
    /** When doing the upward pass */
    long long int getDataUp() const {
        return this->dataUp;
    }
    /** When doing the upward pass */
    void setDataUp(const long long int inData){
        this->dataUp = inData;
    }
    /** When doing the downard pass */
    long long int getDataDown() const {
        return this->dataDown;
    }
    /** When doing the downard pass */
    void setDataDown(const long long int inData){
        this->dataDown = inData;
    }
};


#endif //FTESTCELL_HPP


