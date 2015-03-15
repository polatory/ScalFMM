#ifndef FTESTCELLPOD_HPP
#define FTESTCELLPOD_HPP

#include "../Utils/FGlobal.hpp"
#include "../Containers/FTreeCoordinate.hpp"
#include "FStarPUDefaultAlign.hpp"

struct alignas(FStarPUDefaultAlign::StructAlign) FTestCellPODCore {
    MortonIndex mortonIndex;
    int coordinates[3];
    long long int dataUp, dataDown;
};

class alignas(FStarPUDefaultAlign::StructAlign) FTestCellPOD {
protected:
    FTestCellPODCore data;

public:
    FTestCellPOD() {
        data.mortonIndex = (0);
        data.dataUp = (0);
        data.dataDown = (0);
        data.coordinates[0] = 0;
        data.coordinates[1] = 0;
        data.coordinates[2] = 0;
    }

    /** To get the morton index */
    MortonIndex getMortonIndex() const {
        return data.mortonIndex;
    }

    /** To set the morton index */
    void setMortonIndex(const MortonIndex inMortonIndex) {
        data.mortonIndex = inMortonIndex;
    }

    /** To get the position */
    FTreeCoordinate getCoordinate() const {
        return FTreeCoordinate(data.coordinates[0],
                data.coordinates[1], data.coordinates[2]);
    }

    /** To set the position */
    void setCoordinate(const FTreeCoordinate& inCoordinate) {
        data.coordinates[0] = inCoordinate.getX();
        data.coordinates[1] = inCoordinate.getY();
        data.coordinates[2] = inCoordinate.getZ();
    }

    /** To set the position from 3 FReals */
    void setCoordinate(const int inX, const int inY, const int inZ) {
        data.coordinates[0] = inX;
        data.coordinates[1] = inY;
        data.coordinates[2] = inZ;
    }

    /** When doing the upward pass */
    long long int getDataUp() const {
        return data.dataUp;
    }
    /** When doing the upward pass */
    void setDataUp(const long long int inData){
        data.dataUp = inData;
    }
    /** When doing the downard pass */
    long long int getDataDown() const {
        return data.dataDown;
    }
    /** When doing the downard pass */
    void setDataDown(const long long int inData){
        data.dataDown = inData;
    }

    /** Make it like the begining */
    void resetToInitialState(){
        data.dataDown = 0;
        data.dataUp   = 0;
    }

    /////////////////////////////////////////////////

    /** Save the current cell in a buffer */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        buffer << data.mortonIndex << data.coordinates[0]
               << data.coordinates[1] << data.coordinates[2];
        buffer << data.dataDown << data.dataUp;
    }

    /** Restore the current cell from a buffer */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        buffer >> data.mortonIndex >> data.coordinates[0]
               >> data.coordinates[1] >> data.coordinates[2];
        buffer >> data.dataDown >> data.dataUp;
    }

    int getSavedSize() const {
        return int(sizeof(data.mortonIndex) + sizeof(data.coordinates[0]) +
                sizeof(data.coordinates[1]) + sizeof(data.coordinates[2]) +
                sizeof(data.dataDown) + sizeof(data.dataUp));
    }

    /////////////////////////////////////////////////

    /** Serialize only up data in a buffer */
    template <class BufferWriterClass>
    void serializeUp(BufferWriterClass& buffer) const {
        buffer << data.dataUp;
    }
    /** Deserialize only up data in a buffer */
    template <class BufferReaderClass>
    void deserializeUp(BufferReaderClass& buffer){
        buffer >> data.dataUp;
    }

    /** Serialize only down data in a buffer */
    template <class BufferWriterClass>
    void serializeDown(BufferWriterClass& buffer) const {
        buffer << data.dataDown;
    }
    /** Deserialize only up data in a buffer */
    template <class BufferReaderClass>
    void deserializeDown(BufferReaderClass& buffer){
        buffer >> data.dataDown;
    }

    int getSavedSizeDown() {
        return int(sizeof(long long int));
    }

    int getSavedSizeUp() {
        return int(sizeof(long long int));
    }
};

static_assert(sizeof(FTestCellPODCore) == sizeof(FTestCellPOD), "Core should be equal to cell class size");

#endif // FTESTCELLPOD_HPP

