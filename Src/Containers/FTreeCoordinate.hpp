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
#ifndef FTREECOORDINATE_HPP
#define FTREECOORDINATE_HPP


#include "../Utils/FGlobal.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTreeCoordinate
* Please read the license
*
* This class represents tree coordinate. It is used to save
* the position in "box unit" (not system/space unit!).
* It is directly related to morton index
*/
class FTreeCoordinate{
private:
    int x;	//< x box-th position
    int y;	//< y box-th position
    int z;	//< z box-th position

public:	
    /** Default constructor (position = {0,0,0})*/
    FTreeCoordinate() : x(0), y(0), z(0){
    }

    /**
        * Default constructor
        * @param inX the x
        * @param inY the y
        * @param inZ the z
        */
    explicit FTreeCoordinate(const int inX,const int inY,const int inZ)
        : x(inX), y(inY), z(inZ){
    }

    /**
	* Copy constructor
	* @param other the source class to copy
	*/
    FTreeCoordinate(const FTreeCoordinate& other):x(other.x), y(other.y), z(other.z){
    }

    /**
        * Copy constructor
        * @param other the source class to copy
        */
    FTreeCoordinate(const FTreeCoordinate& other, const int inOffset)
        : x(other.x + inOffset), y(other.y + inOffset), z(other.z + inOffset){
    }

    /**
	* Copy constructor
	* @param other the source class to copy
	* @return this a reference to the current object
	*/
    FTreeCoordinate& operator=(const FTreeCoordinate& other){
        this->x = other.x;
        this->y = other.y;
        this->z = other.z;
        return *this;
    }

    /**
	* Position setter
        * @param inX the new x
        * @param inY the new y
        * @param inZ the new z
	*/
    void setPosition(const int inX,const int inY,const int inZ){
        this->x = inX;
        this->y = inY;
        this->z = inZ;
    }

    /**
	* X Getter
	* @return this->x
	*/
    int getX() const{
        return this->x;
    }

    /**
	* Y Getter
	* @return this->y
	*/
    int getY() const{
        return this->y;
    }

    /**
	* Z Getter
	* @return this->z
	*/
    int getZ() const{
        return this->z;
    }

    /**
	* X Setter, simply change x position
	* @param the new x
	*/
    void setX(const int inX){
        this->x = inX;
    }

    /**
	* Y Setter, simply change y position
	* @param the new y
	*/
    void setY(const int inY){
        this->y = inY;
    }

    /**
	* Z Setter, simply change z position
	* @param the new z
	*/
    void setZ(const int inZ){
        this->z = inZ;
    }

    /**
	* To get the morton index of the current position
	* @complexity inLevel
	* @param inLevel the level of the component
	* @return morton index
	*/
    MortonIndex getMortonIndex(const int inLevel) const{
        MortonIndex index = 0x0LL;
        MortonIndex mask = 0x1LL;
        // the ordre is xyz.xyz...
        MortonIndex mx = this->x << 2;
        MortonIndex my = this->y << 1;
        MortonIndex mz = this->z;

        for(int indexLevel = 0; indexLevel < inLevel ; ++indexLevel){
            index |= (mz & mask);
            mask <<= 1;
            index |= (my & mask);
            mask <<= 1;
            index |= (mx & mask);
            mask <<= 1;

            mz <<= 2;
            my <<= 2;
            mx <<= 2;
        }

        return index;
    }

    /** This function set the position of the current object using a morton index
          * @param inIndex the morton index to compute position
          * @param the level of the morton index
          */
    void setPositionFromMorton(MortonIndex inIndex, const int inLevel){
        MortonIndex mask = 0x1LL;

        this->x = 0;
        this->y = 0;
        this->z = 0;

        for(int indexLevel = 0; indexLevel < inLevel ; ++indexLevel){
            z |= int(inIndex & mask);
            inIndex >>= 1;
            y |= int(inIndex & mask);
            inIndex >>= 1;
            x |= int(inIndex & mask);

            mask <<= 1;
        }

    }

    /** Test equal operator
          * @param other the coordinate to compare
          * @return true if other & current object have same position
          */
    bool operator==(const FTreeCoordinate& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    /** To test difference
      *
      */
    bool operator!=(const FTreeCoordinate& other) const{
        return x != other.x || y != other.y || z != other.z;
    }
};



#endif //FTREECOORDINATE_HPP


