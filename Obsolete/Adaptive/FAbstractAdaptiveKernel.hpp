#ifndef FABSTRACTADAPTIVEKERNEL
#define FABSTRACTADAPTIVEKERNEL
// Keep in private GIT
// @SCALFMM_PRIVATE
/**
 * This class represent the method that an adaptive kernel must implement.
 * There are two kinds of operators, the first one represent computation and the others
 * should return the cretiria to know when the P2M should be performed.
 */
template <class CellClass, class ContainerClass>
class FAbstractAdaptiveKernel {
public:
    virtual ~FAbstractAdaptiveKernel(){
    }

    virtual void P2M(CellClass* const pole, const int cellLevel, const ContainerClass* const particles) = 0;

    virtual void M2M(CellClass* const pole, const int poleLevel, const CellClass* const subCell, const int subCellLevel) = 0;

    virtual void P2L(CellClass* const local, const int localLevel, const ContainerClass* const particles) = 0;

    virtual void M2L(CellClass* const local, const int localLevel, const CellClass* const aNeighbor, const int neighborLevel) = 0;

    virtual void M2P(const CellClass* const pole, const int poleLevel, ContainerClass* const particles) = 0;

    virtual void L2L(const CellClass* const local, const int localLevel, CellClass* const subCell, const int subCellLevel) = 0;

    virtual void L2P(const CellClass* const local, const int cellLevel, ContainerClass* const particles) = 0;

    virtual void P2P(ContainerClass* target, const ContainerClass* sources) = 0;

    virtual bool preferP2M(const ContainerClass* const particles) = 0;
    virtual bool preferP2M(const int atLevel, const ContainerClass*const particles[], const int nbContainers) = 0;
};

#endif
