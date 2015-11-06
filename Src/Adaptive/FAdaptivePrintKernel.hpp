#ifndef FADAPTIVEPRINTKERNEL_HPP
#define FADAPTIVEPRINTKERNEL_HPP
// Keep in private GIT
// @SCALFMM_PRIVATE
#include <iostream>

#include "FAbstractAdaptiveKernel.hpp"
#include "../Components/FAbstractKernels.hpp"

/**
 * This class simply print the interactions
 */
template <class CellClass, class ContainerClass>
class FAdaptivePrintKernel : public FAbstractKernels<CellClass,ContainerClass>, public FAbstractAdaptiveKernel<CellClass,ContainerClass>  {
    const int p2mTresh;
public:
    FAdaptivePrintKernel(const int inP2mTresh = 10)
        : p2mTresh(inP2mTresh) {
    }

    virtual ~FAdaptivePrintKernel(){
    }

    void P2M(CellClass* const pole, const ContainerClass* const particles) override {
        std::cout << "Usual] P2M Idx = " << pole->getMortonIndex() << " and @" << particles << " with " << particles->getNbParticles() << "\n";
    }

    void M2M(CellClass* const FRestrict pole, const CellClass*const FRestrict *const FRestrict child , const int level) override {
        std::cout << "Usual] M2M Idx = " << pole->getMortonIndex() << " at level " << level  << " with\n";
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            if(child[idxChild]){
                std::cout << "\t Idx " << child[idxChild]->getMortonIndex() << "\n";
            }
        }
    }

    void M2L(CellClass* const FRestrict local, const CellClass* interactions[], const int /*positions*/[], const int counter, const int level) override {
        std::cout << "Usual] M2L Idx = " << local->getMortonIndex() << " at level " << level  << " with\n";
        for(int idxInter = 0 ; idxInter < counter ; ++idxInter){
                std::cout << "\t Idx " << interactions[idxInter]->getMortonIndex() << "\n";
        }
    }

    void L2L(const CellClass* const FRestrict local, CellClass* FRestrict *const FRestrict child , const int level) override {
        std::cout << "Usual] L2L Idx = " << local->getMortonIndex() << " at level " << level  << " with\n";
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            if(child[idxChild]){
                std::cout << "\t Idx " << child[idxChild]->getMortonIndex() << "\n";
            }
        }
    }

    void L2P(const CellClass* const local, ContainerClass* const particles) override {
        std::cout << "Usual] L2P Idx = " << local->getMortonIndex() << " and @" << particles << " with " << particles->getNbParticles() << "\n";
    }

    void P2P(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict ,
                     ContainerClass* const neighs[], const int /*positions*/[], const int size) override {
        std::cout << "Usual] P2P @" << targets << " has " << targets->getNbParticles() << " with\n";
        for(int idxNeigh = 0 ; idxNeigh < size ; ++idxNeigh){
                std::cout << "\t @" << neighs[idxNeigh]<< " has " << neighs[idxNeigh]->getNbParticles() << "\n";
        }
    }

    void P2PRemote(const FTreeCoordinate& ,
                     ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict ,
                     const ContainerClass* const neighs[], const int /*positions*/[], const int size) override {
        std::cout << "Usual] P2P remote @" << targets << " has " << targets->getNbParticles() << " with\n";
        for(int idxNeigh = 0 ; idxNeigh < size ; ++idxNeigh){
                std::cout << "\t @" << neighs[idxNeigh]<< " has " << neighs[idxNeigh]->getNbParticles() << "\n";
        }
    }

    void P2M(CellClass* const pole, const int cellLevel, const ContainerClass* const particles) override{
        std::cout << "Adaptive] P2M Idx = " << pole->getMortonIndex() << " and @" << particles << " with " << particles->getNbParticles() << ", cell at level " << cellLevel << "\n";
    }

    void M2M(CellClass* const pole, const int poleLevel, const CellClass* const subCell, const int subCellLevel) override{
        std::cout << "Adaptive] M2M Idx = " << pole->getMortonIndex() << " at level " << poleLevel
                  << " with " << subCell->getMortonIndex() << " from level " << subCellLevel << "\n";
    }

    void P2L(CellClass* const local, const int localLevel, const ContainerClass* const particles) override{
        std::cout << "Adaptive] P2L Local Idx = " << local->getMortonIndex() << " with " << particles->getNbParticles() << ", cell at level " << localLevel << "\n";
    }

    void M2L(CellClass* const local, const int localLevel, const CellClass* const aNeighbor, const int neighborLevel)  override{
        std::cout << "Adaptive] L2L Idx = " << local->getMortonIndex() << " at level " << localLevel
                  << " with " << aNeighbor->getMortonIndex() << " from level " << neighborLevel << "\n";
    }

    void M2P(const CellClass* const pole, const int poleLevel, ContainerClass* const particles) override{
        std::cout << "Adaptive] M2P Idx = " << pole->getMortonIndex() << " with " << particles->getNbParticles() << ", cell at level " << poleLevel << "\n";
    }

    void L2L(const CellClass* const local, const int localLevel, CellClass* const subCell, const int subCellLevel)  override{
        std::cout << "Adaptive] L2L Idx = " << local->getMortonIndex() << " at level " << localLevel
                  << " with " << subCell->getMortonIndex() << " from level " << subCellLevel << "\n";
    }

    void L2P(const CellClass* const local, const int cellLevel, ContainerClass* const particles) override{
        std::cout << "Adaptive] L2P Idx = " << local->getMortonIndex() << " with " << particles->getNbParticles() << ", cell at level " << cellLevel << "\n";
    }

    void P2P(ContainerClass* targets, const ContainerClass* sources)  override{
        std::cout << "Adaptive] P2P between @" << targets << " with " << targets->getNbParticles()
                     << " and @" << sources << " with " << sources->getNbParticles() << "\n";
    }

    bool preferP2M(const ContainerClass* const particles)  override{
        return particles->getNbParticles() > p2mTresh ;
    }

    bool preferP2M(const int atLevel, const ContainerClass*const particles[], const int nbContainers)  override{
        int count = 0;
        for(int idx = 0 ; idx < nbContainers ; ++idx){
            count += particles[idx]->getNbParticles();
        }
        return count > p2mTresh;
    }
};


#endif // FADAPTIVEPRINTKERNEL_HPP
