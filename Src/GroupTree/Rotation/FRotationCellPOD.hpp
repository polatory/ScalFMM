//@SCALFMM_PRIVATE
#ifndef FROTATIONCELLPOD_HPP
#define FROTATIONCELLPOD_HPP


#include "../../Utils/FComplex.hpp"
#include "../../Utils/FMemUtils.hpp"
#include "../../Utils/FGlobal.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../StarPUUtils/FStarPUDefaultAlign.hpp"

struct alignas(FStarPUDefaultAlign::StructAlign) FRotationCellPODCore {
    MortonIndex mortonIndex;
    int coordinates[3];
};

template <int P>
struct alignas(FStarPUDefaultAlign::StructAlign) FRotationCellPODPole {
    //< Size of multipole vector
    static const int MultipoleSize = ((P+2)*(P+1))/2; // Artimethique suite (n+1)*n/2
    //< Multipole vector (static memory)
    FComplex multipole_exp[MultipoleSize]; //< For multipole extenssion
};

template <int P>
struct alignas(FStarPUDefaultAlign::StructAlign) FRotationCellPODLocal {
    //< Size of local vector
    static const int LocalSize = ((P+2)*(P+1))/2;     // Artimethique suite (n+1)*n/2
    //< Local vector (static memory)
    FComplex local_exp[LocalSize];         //< For local extenssion
};


template <int P>
class FRotationCellPOD
{
    FRotationCellPODCore* symb;
    FRotationCellPODPole<P>* up;
    FRotationCellPODLocal<P>* down;

public:
    FRotationCellPOD(FRotationCellPODCore* inSymb, FRotationCellPODPole<P>* inUp,
              FRotationCellPODLocal<P>* inDown): symb(inSymb), up(inUp), down(inDown){
    }

    FRotationCellPOD()
        : symb(nullptr), up(nullptr), down(nullptr){
    }

    /** To get the morton index */
    MortonIndex getMortonIndex() const {
        return symb->mortonIndex;
    }

    /** To set the morton index */
    void setMortonIndex(const MortonIndex inMortonIndex) {
        symb->mortonIndex = inMortonIndex;
    }

    /** To get the position */
    FTreeCoordinate getCoordinate() const {
        return FTreeCoordinate(symb->coordinates[0],
                symb->coordinates[1], symb->coordinates[2]);
    }

    /** To set the position */
    void setCoordinate(const FTreeCoordinate& inCoordinate) {
        symb->coordinates[0] = inCoordinate.getX();
        symb->coordinates[1] = inCoordinate.getY();
        symb->coordinates[2] = inCoordinate.getZ();
    }

    /** To set the position from 3 FReals */
    void setCoordinate(const int inX, const int inY, const int inZ) {
        symb->coordinates[0] = inX;
        symb->coordinates[1] = inY;
        symb->coordinates[2] = inZ;
    }

    /** Get Multipole */
    const FComplex* getMultipole() const
    {	return up->multipole_exp;
    }
    /** Get Local */
    const FComplex* getLocal() const{
        return down->local_exp;
    }

    /** Get Multipole */
    FComplex* getMultipole(){
        return up->multipole_exp;
    }
    /** Get Local */
    FComplex* getLocal(){
        return down->local_exp;
    }

    /** To get the leading dim of a vec */
    int getVectorSize() const{
        return down->VectorSize;
    }

    /** Make it like the begining */
    void resetToInitialState(){
        for(int idx = 0 ; idx < up->MultipoleSize ; ++idx){
            up->multipole_exp[idx].setRealImag(FReal(0.0), FReal(0.0));
        }
        for(int idx = 0 ; idx < down->LocalSize ; ++idx){
            down->local_exp[idx].setRealImag(FReal(0.0), FReal(0.0));
        }
    }
};


#endif // FROTATIONCELLPOD_HPP

