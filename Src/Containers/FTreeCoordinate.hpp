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
#ifndef FTREECOORDINATE_HPP
#define FTREECOORDINATE_HPP


#include "../Utils/FGlobal.hpp"
#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"


/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTreeCoordinate
* Please read the license
*
* This class represents tree coordinate. It is used to save
* the position in "box unit" (not system/space unit!).
* It is directly related to morton index, as interleaves
* bits from this coordinate make the morton index
*/
class FTreeCoordinate{
private:
    int data[3];	//< all box-th position

public:	
    /** Default constructor (position = {0,0,0})*/
    FTreeCoordinate() {
        data[0] = data[1] = data[2] = 0;
    }

    /**
        * Default constructor
        * @param inX the x
        * @param inY the y
        * @param inZ the z
        */
    explicit FTreeCoordinate(const int inX,const int inY,const int inZ) {
         data[0] = inX;
         data[1] = inY;
         data[2] = inZ;
    }

    explicit FTreeCoordinate(const int inPosition[3]) {
        data[0] = inPosition[0];
        data[1] = inPosition[1];
        data[2] = inPosition[2];
    }

    /**
	* Copy constructor
	* @param other the source class to copy
	*/
    FTreeCoordinate(const FTreeCoordinate& other) {
        data[0] = other.data[0];
        data[1] = other.data[1];
        data[2] = other.data[2];
    }

    /**
        * Copy constructor
        * @param other the source class to copy
        */
    FTreeCoordinate(const FTreeCoordinate& other, const int inOffset) {
        data[0] = other.data[0] + inOffset;
        data[1] = other.data[1] + inOffset;
        data[2] = other.data[2] + inOffset;
    }

    /**
	* Copy constructor
	* @param other the source class to copy
	* @return this a reference to the current object
	*/
    FTreeCoordinate& operator=(const FTreeCoordinate& other){
        data[0] = other.data[0];
        data[1] = other.data[1];
        data[2] = other.data[2];
        return *this;
    }

    /**
	* Position setter
        * @param inX the new x
        * @param inY the new y
        * @param inZ the new z
	*/
    void setPosition(const int inX,const int inY,const int inZ){
        data[0] = inX;
        data[1] = inY;
        data[2] = inZ;
    }

    /**
	* X Getter
        * @return data[0]
	*/
    int getX() const{
        return data[0];
    }

    /**
	* Y Getter
        * @return data[1]
	*/
    int getY() const{
        return data[1];
    }

    /**
	* Z Getter
        * @return data[2]
	*/
    int getZ() const{
        return data[2];
    }

    /**
	* X Setter, simply change x position
	* @param the new x
	*/
    void setX(const int inX){
        data[0] = inX;
    }

    /**
	* Y Setter, simply change y position
	* @param the new y
	*/
    void setY(const int inY){
        data[1] = inY;
    }

    /**
	* Z Setter, simply change z position
	* @param the new z
	*/
    void setZ(const int inZ){
        data[2] = inZ;
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
        MortonIndex mx = data[0] << 2;
        MortonIndex my = data[1] << 1;
        MortonIndex mz = data[2];

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

        data[0] = 0;
        data[1] = 0;
        data[2] = 0;

        for(int indexLevel = 0; indexLevel < inLevel ; ++indexLevel){
            data[2] |= int(inIndex & mask);
            inIndex >>= 1;
            data[1] |= int(inIndex & mask);
            inIndex >>= 1;
            data[0] |= int(inIndex & mask);

            mask <<= 1;
        }

    }

    /** Test equal operator
          * @param other the coordinate to compare
          * @return true if other & current object have same position
          */
    bool operator==(const FTreeCoordinate& other) const {
        return data[0] == other.data[0] && data[1] == other.data[1] && data[2] == other.data[2];
    }

    /** To test difference
      *
      */
    bool operator!=(const FTreeCoordinate& other) const{
        return data[0] != other.data[0] || data[1] != other.data[1] || data[2] != other.data[2];
    }

    /**
     * Operator stream FTreeCoordinate to std::ostream
     * This can be used to simply write out a tree coordinate
     * @param[in,out] output where to write the coordinate
     * @param[in] inCoordinate the coordinate to write out
     * @return the output for multiple << operators
     */    
    template <class StreamClass>
    friend StreamClass& operator<<(StreamClass& output, const FTreeCoordinate& inCoordinate){
        output << "(" <<  inCoordinate.getX() << ", " << inCoordinate.getY() << ", " << inCoordinate.getZ() <<")";
        return output;  // for multiple << operators.
    }

    /** Save current object */
    void save(FBufferWriter& buffer) const {
        buffer << data[0] << data[1] << data[2];
    }
    /** Retrieve current object */
    void restore(FBufferReader& buffer) {
        buffer >> data[0] >> data[1] >> data[2];
    }
};



#endif //FTREECOORDINATE_HPP


