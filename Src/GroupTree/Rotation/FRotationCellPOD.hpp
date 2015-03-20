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

    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    template <class BufferWriterClass>
    void serializeUp(BufferWriterClass& buffer) const{
        buffer.write(up->multipole_exp, up->MultipoleSize);
    }
    template <class BufferReaderClass>
    void deserializeUp(BufferReaderClass& buffer){
        buffer.fillArray(up->multipole_exp, up->MultipoleSize);
    }

    template <class BufferWriterClass>
    void serializeDown(BufferWriterClass& buffer) const{
        buffer.write(down->local_exp, down->LocalSize);
    }
    template <class BufferReaderClass>
    void deserializeDown(BufferReaderClass& buffer){
        buffer.fillArray(down->local_exp, down->LocalSize);
    }

    ///////////////////////////////////////////////////////
    // to extend Serializable
    ///////////////////////////////////////////////////////
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        buffer << symb->mortonIndex << symb->coordinates[0]
               << symb->coordinates[1] << symb->coordinates[2];
        buffer.write(up->multipole_exp, up->MultipoleSize);
        buffer.write(down->local_exp, down->LocalSize);
    }
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        buffer >> symb->mortonIndex >> symb->coordinates[0]
               >> symb->coordinates[1] >> symb->coordinates[2];
        buffer.fillArray(up->multipole_exp, up->MultipoleSize);
        buffer.fillArray(down->local_exp, down->LocalSize);
    }

    int getSavedSize() const {
        return int(sizeof(FComplex)*(up->MultipoleSize + down->LocalSize) + sizeof(symb->mortonIndex) + sizeof(symb->coordinates[0]) +
                sizeof(symb->coordinates[1]) + sizeof(symb->coordinates[2]));
    }

    int getSavedSizeUp() const {
        return ((int) sizeof(FComplex)) * (up->MultipoleSize);
    }

    int getSavedSizeDown() const {
        return ((int) sizeof(FComplex)) * (down->LocalSize);
    }

};


#endif // FROTATIONCELLPOD_HPP

