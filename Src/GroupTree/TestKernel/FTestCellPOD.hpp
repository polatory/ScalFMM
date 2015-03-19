// @SCALFMM_PRIVATE
#ifndef FTESTCELLPOD_HPP
#define FTESTCELLPOD_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../StarPUUtils/FStarPUDefaultAlign.hpp"

struct alignas(FStarPUDefaultAlign::StructAlign) FTestCellPODCore {
    FTestCellPODCore(){
        mortonIndex = (0);
        coordinates[0] = 0;
        coordinates[1] = 0;
        coordinates[2] = 0;
    }

    MortonIndex mortonIndex;
    int coordinates[3];
};

typedef long long FTestCellPODData;

class FTestCellPOD {
protected:
    FTestCellPODCore* symb;
    FTestCellPODData* up;
    FTestCellPODData* down;

public:
    FTestCellPOD(FTestCellPODCore* inSymb, FTestCellPODData* inUp,
                 FTestCellPODData* inDown)
        : symb(inSymb), up(inUp), down(inDown){
    }
    FTestCellPOD()
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

    /** When doing the upward pass */
    long long int getDataUp() const {
        return (*up);
    }
    /** When doing the upward pass */
    void setDataUp(const long long int inData){
        (*up) = inData;
    }
    /** When doing the downard pass */
    long long int getDataDown() const {
        return (*down);
    }
    /** When doing the downard pass */
    void setDataDown(const long long int inData){
        (*down) = inData;
    }

    /** Make it like the begining */
    void resetToInitialState(){
        (*down) = 0;
        (*up)   = 0;
    }

    /////////////////////////////////////////////////

    /** Save the current cell in a buffer */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        buffer << symb->mortonIndex << symb->coordinates[0]
               << symb->coordinates[1] << symb->coordinates[2];
        buffer << (*down) << (*up);
    }

    /** Restore the current cell from a buffer */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        buffer >> symb->mortonIndex >> symb->coordinates[0]
               >> symb->coordinates[1] >> symb->coordinates[2];
        buffer >> (*down) >> (*up);
    }

    int getSavedSize() const {
        return int(sizeof(symb->mortonIndex) + sizeof(symb->coordinates[0]) +
                sizeof(symb->coordinates[1]) + sizeof(symb->coordinates[2]) +
                sizeof((*down)) + sizeof((*up)));
    }

    /////////////////////////////////////////////////

    /** Serialize only up data in a buffer */
    template <class BufferWriterClass>
    void serializeUp(BufferWriterClass& buffer) const {
        buffer << (*up);
    }
    /** Deserialize only up data in a buffer */
    template <class BufferReaderClass>
    void deserializeUp(BufferReaderClass& buffer){
        buffer >> (*up);
    }

    /** Serialize only down data in a buffer */
    template <class BufferWriterClass>
    void serializeDown(BufferWriterClass& buffer) const {
        buffer << (*down);
    }
    /** Deserialize only up data in a buffer */
    template <class BufferReaderClass>
    void deserializeDown(BufferReaderClass& buffer){
        buffer >> (*down);
    }

    int getSavedSizeDown() {
        return int(sizeof(FTestCellPODData));
    }

    int getSavedSizeUp() {
        return int(sizeof(FTestCellPODData));
    }
};


#endif // FTESTCELLPOD_HPP

