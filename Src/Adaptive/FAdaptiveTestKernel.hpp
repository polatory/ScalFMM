#ifndef FADAPTIVETESTKERNEL_HPP
#define FADAPTIVETESTKERNEL_HPP
// Keep in private GIT
// @SCALFMM_PRIVATE
#include "FAbstractAdaptiveKernel.hpp"
#include "../Components/FTestKernels.hpp"

template< class CellClass, class ContainerClass>
class FAdaptiveTestKernel : public FTestKernels<CellClass, ContainerClass>, public FAbstractAdaptiveKernel<CellClass, ContainerClass> {
    const int p2mThresh;

public:
    FAdaptiveTestKernel(const int inP2mThresh = 10)
        : p2mThresh(inP2mThresh){
    }

    using FTestKernels<CellClass, ContainerClass>::P2M;
    using FTestKernels<CellClass, ContainerClass>::M2M;
    using FTestKernels<CellClass, ContainerClass>::M2L;
    using FTestKernels<CellClass, ContainerClass>::L2L;
    using FTestKernels<CellClass, ContainerClass>::L2P;
    using FTestKernels<CellClass, ContainerClass>::P2P;

    void P2M(CellClass* const pole, const int /*cellLevel*/, const ContainerClass* const particles) override {
        pole->setDataUp(pole->getDataUp() + particles->getNbParticles());
    }

    void M2M(CellClass* const pole, const int /*poleLevel*/, const CellClass* const subCell, const int /*subCellLevel*/) override {
        pole->setDataUp(pole->getDataUp() + subCell->getDataUp());
    }

    void P2L(CellClass* const local, const int /*localLevel*/, const ContainerClass* const particles) override {
        local->setDataDown(local->getDataDown() + particles->getNbParticles());
    }

    void M2L(CellClass* const local, const int /*localLevel*/, const CellClass* const aNeighbor, const int /*neighborLevel*/) override {
        local->setDataDown(local->getDataDown() + aNeighbor->getDataUp());
    }

    void M2P(const CellClass* const pole, const int /*poleLevel*/, ContainerClass* const particles) override {
        long long int*const particlesAttributes = particles->getDataDown();
        for(FSize idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += pole->getDataUp();
        }
    }

    void L2L(const CellClass* const local, const int /*localLevel*/, CellClass* const subCell, const int /*subCellLevel*/) override {
        subCell->setDataDown(local->getDataDown() + subCell->getDataDown());
    }

    void L2P(const CellClass* const local, const int /*cellLevel*/, ContainerClass* const particles)  override {
        long long int*const particlesAttributes = particles->getDataDown();
        for(FSize idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += local->getDataDown();
        }
    }

    void P2P(ContainerClass* target, const ContainerClass* sources)  override {
        long long int*const particlesAttributes = target->getDataDown();
        for(FSize idxPart = 0 ; idxPart < target->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += sources->getNbParticles();
        }
    }

    bool preferP2M(const ContainerClass* const particles) override {
        return particles->getNbParticles() > p2mThresh;
    }
    bool preferP2M(const int /*atLevel*/, const ContainerClass*const particles[], const int nbContainers) override {
        FSize counterParticles = 0;
        for(int idxContainer = 0 ; idxContainer < nbContainers ; ++idxContainer){
            counterParticles += particles[idxContainer]->getNbParticles();
        }
        return counterParticles > p2mThresh;
    }
};

#endif // FADAPTIVETESTKERNEL_HPP
