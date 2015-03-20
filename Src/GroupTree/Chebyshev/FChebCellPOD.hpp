#ifndef FCHEBCELLPOD_HPP
#define FCHEBCELLPOD_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../StarPUUtils/FStarPUDefaultAlign.hpp"
#include "../../Kernels/Chebyshev/FChebTensor.hpp"

struct alignas(FStarPUDefaultAlign::StructAlign) FChebCellPODCore {
    MortonIndex mortonIndex;
    int coordinates[3];
};

template <int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
struct alignas(FStarPUDefaultAlign::StructAlign) FChebCellPODPole {
    static const int VectorSize = TensorTraits<ORDER>::nnodes * 2;
    FReal multipole_exp[NRHS * NVALS * VectorSize];
};

template <int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
struct alignas(FStarPUDefaultAlign::StructAlign) FChebCellPODLocal {
    static const int VectorSize = TensorTraits<ORDER>::nnodes * 2;
    FReal local_exp[NLHS * NVALS * VectorSize]; //< Local expansion
};


template <int ORDER, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FChebCellPOD
{
    // nnodes = ORDER^3
    // we multiply by 2 because we store the  Multipole expansion end the compressed one.
    static const int VectorSize = TensorTraits<ORDER>::nnodes * 2;

    FChebCellPODCore* symb;
    FChebCellPODPole<ORDER,NRHS,NLHS,NVALS>* up;
    FChebCellPODLocal<ORDER,NRHS,NLHS,NVALS>* down;

public:
    FChebCellPOD(FChebCellPODCore* inSymb, FChebCellPODPole<ORDER,NRHS,NLHS,NVALS>* inUp,
              FChebCellPODLocal<ORDER,NRHS,NLHS,NVALS>* inDown): symb(inSymb), up(inUp), down(inDown){
    }

    FChebCellPOD()
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
    const FReal* getMultipole(const int inRhs) const
    {	return up->multipole_exp + inRhs*VectorSize;
    }
    /** Get Local */
    const FReal* getLocal(const int inRhs) const{
        return down->local_exp + inRhs*VectorSize;
    }

    /** Get Multipole */
    FReal* getMultipole(const int inRhs){
        return up->multipole_exp + inRhs*VectorSize;
    }
    /** Get Local */
    FReal* getLocal(const int inRhs){
        return down->local_exp + inRhs*VectorSize;
    }

    /** To get the leading dim of a vec */
    int getVectorSize() const{
        return VectorSize;
    }

    /** Make it like the begining */
    void resetToInitialState(){
        memset(up->multipole_exp, 0, sizeof(FReal) * NRHS * NVALS * VectorSize);
        memset(down->local_exp,         0, sizeof(FReal) * NLHS * NVALS * VectorSize);
    }

    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    template <class BufferWriterClass>
    void serializeUp(BufferWriterClass& buffer) const{
        buffer.write(up->multipole_exp, VectorSize*NVALS*NRHS);
    }
    template <class BufferReaderClass>
    void deserializeUp(BufferReaderClass& buffer){
        buffer.fillArray(up->multipole_exp, VectorSize*NVALS*NRHS);
    }

    template <class BufferWriterClass>
    void serializeDown(BufferWriterClass& buffer) const{
        buffer.write(down->local_exp, VectorSize*NVALS*NLHS);
    }
    template <class BufferReaderClass>
    void deserializeDown(BufferReaderClass& buffer){
        buffer.fillArray(down->local_exp, VectorSize*NVALS*NLHS);
    }

    ///////////////////////////////////////////////////////
    // to extend Serializable
    ///////////////////////////////////////////////////////
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        buffer << symb->mortonIndex << symb->coordinates[0]
               << symb->coordinates[1] << symb->coordinates[2];
        buffer.write(up->multipole_exp, VectorSize*NVALS*NRHS);
        buffer.write(down->local_exp, VectorSize*NVALS*NLHS);
    }
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        buffer >> symb->mortonIndex >> symb->coordinates[0]
               >> symb->coordinates[1] >> symb->coordinates[2];
        buffer.fillArray(up->multipole_exp, VectorSize*NVALS*NRHS);
        buffer.fillArray(down->local_exp, VectorSize*NVALS*NLHS);
    }

    int getSavedSize() const {
        return int(sizeof(FReal) * VectorSize*(NRHS+NLHS)*NVALS + sizeof(symb->mortonIndex) + sizeof(symb->coordinates[0]) +
                sizeof(symb->coordinates[1]) + sizeof(symb->coordinates[2]));
    }

    int getSavedSizeUp() const {
        return int(sizeof(FReal)) * VectorSize*(NRHS)*NVALS;
    }

    int getSavedSizeDown() const {
        return int(sizeof(FReal)) * VectorSize*(NLHS)*NVALS;
    }

};


#endif // FCHEBCELLPOD_HPP

