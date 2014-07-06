#ifndef FADAPTATIVEKERNELWRAPPER_HPP
#define FADAPTATIVEKERNELWRAPPER_HPP

#include "../Components/FAbstractKernels.hpp"
#include "../Containers/FVector.hpp"

#include "FAdaptativeCell.hpp"



/**
 * This class is a wrapper to use the usual algorithm but in a adaptative way.
 * The real computational kernel should be given in template.
 * It must propose the usual FMM Kernel operators and also the FAbstractKernels operators.
 *
 * The cells in the tree should respect the FAdaptativeCell.
 */
template <class RealAdaptativeKernel, class CellClass, class ContainerClass >
class FAdaptativeKernelWrapper : public FAbstractKernels<FAdaptativeCell<CellClass, ContainerClass>, ContainerClass>{
protected:
    RealAdaptativeKernel kernel;
public:
    template<typename... KernelParams>
    FAdaptativeKernelWrapper(KernelParams... parameters) : kernel(parameters...){
    }

    RealAdaptativeKernel& getKernel(){
        return kernel;
    }

    const RealAdaptativeKernel& getKernel() const{
        return kernel;
    }

    /** P2M is performed only if the kernel says that is it better than waiting for future development
      * A leaf cell contains the particles container if the development is not made,
      * else it has a multipole component.
      */
    void P2M(FAdaptativeCell<CellClass, ContainerClass>* const pole, const ContainerClass* const particles)  override {
        // Leaf is always adaptative
        pole->setAdaptative(true);
        if( kernel.preferP2M(particles) ){
            // If it is better to compute the P2M at this level
            pole->setHaveDevelopment(true);
            kernel.P2M(pole->getRealCell(), particles);
        }
        else{
            // Else simply keep the current leaf
            pole->setHaveDevelopment(false);
            pole->addSubLeaf(particles);
        }
    }

    /** The M2M need to manage all the possibilities and to propagate the data from the leaf.
      * If there is one child, we are not adaptative and keep the lower adaptive cell as target
      * If at least one child has a development (at any level) the current cell store the development of all the children
      * using M2M or P2M from the leaves.
      * If no children have developments we need to test if it is better to merge them in a development
      * or to store all the particles in order to perform the P2M later.
      * Finally if all the children have developments at level just under we can perform normal M2M
      */
    void M2M(FAdaptativeCell<CellClass, ContainerClass>* const FRestrict pole,
             const FAdaptativeCell<CellClass, ContainerClass>*const FRestrict *const FRestrict child, const int inLevel)  override {
        int nbChild       = 0;
        int lastChild     = 0;
        bool onlyParticlesCells         = true;
        FVector<const ContainerClass*> subLeaves;

        // Test all the children
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            if(child[idxChild]){
                nbChild  += 1;
                lastChild = idxChild;
                // We agragate all the leaves from the current child or its adaptative cell
                if(child[idxChild]->isAdaptative() && !child[idxChild]->hasDevelopment()){
                    subLeaves.memocopy( child[idxChild]->getSubLeaves(), child[idxChild]->getNbSubLeaves());
                }
                else if(!child[idxChild]->isAdaptative() && !child[idxChild]->getSubAdaptativeCell()->hasDevelopment()){
                    subLeaves.memocopy( child[idxChild]->getSubAdaptativeCell()->getSubLeaves(),
                                        child[idxChild]->getSubAdaptativeCell()->getNbSubLeaves());
                }
                else{
                    // If a child is made of development
                    onlyParticlesCells = false;
                }
            }
        }
        // We need to agregate if there are only particles and if the kernel says so
        const bool continueToAgregate = (onlyParticlesCells && kernel.preferP2M(inLevel, subLeaves.data(), subLeaves.getSize()));
        if(nbChild == 1){
            // One child means that the cell is not adaptative
            pole->setAdaptative(false);
            pole->setHaveDevelopment(false);
            if(child[lastChild]->isAdaptative()){
                pole->setSubAdaptativeCell(child[lastChild], inLevel + 1);
            }
            else{
                pole->setSubAdaptativeCell(child[lastChild]->getSubAdaptativeCell(),
                                           child[lastChild]->getSubAdaptativeLevel());
            }
        }
        else if(onlyParticlesCells && continueToAgregate){
            // There are only particles and no development to do
            pole->setAdaptative(true);
            pole->setHaveDevelopment(false);
            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if(child[idxChild]){
                    pole->addSubLeaves(child[idxChild]->getSubLeaves(), child[idxChild]->getNbSubLeaves());
                }
            }
        }
        else{
            // There development to do from developments or particles
            pole->setAdaptative(true);
            pole->setHaveDevelopment(true);

            const CellClass* realChild[8] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
            int counterRealLowerCell      = 0;
            // Test each child
            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if(child[idxChild]){
                    if(child[idxChild]->isAdaptative()){
                        if(child[idxChild]->hasDevelopment()){
                            // If it is adaptative and has development than we compute is using usual M2M
                            realChild[idxChild] = child[idxChild]->getRealCell();
                            counterRealLowerCell += 1;
                        }
                        else{
                            // If it is adaptative and has not development than we compute is using P2M
                            for(int idxLeaf = 0 ; idxLeaf < child[idxChild]->getNbSubLeaves() ; ++idxLeaf){
                                kernel.P2M(pole->getRealCell(), inLevel, child[idxChild]->getSubLeaf(idxLeaf));
                            }
                        }
                    }
                    else{
                        // Else we need the adaptative cell
                        const FAdaptativeCell<CellClass, ContainerClass>* lowerAdaptativeCell = child[idxChild]->getSubAdaptativeCell();
                        const int lowerAdaptativeLevel = child[idxChild]->getSubAdaptativeLevel();
                        if(lowerAdaptativeCell->hasDevelopment()){
                            // If it has development we perform a M2M
                            kernel.M2M(pole->getRealCell(), inLevel, lowerAdaptativeCell->getRealCell(),
                                       lowerAdaptativeLevel);
                        }
                        else{
                            // Else we perform P2M
                            for(int idxLeaf = 0 ; idxLeaf < child[idxChild]->getNbSubLeaves() ; ++idxLeaf){
                                kernel.P2M(pole->getRealCell(), inLevel, lowerAdaptativeCell->getSubLeaf(idxLeaf));
                            }
                        }
                    }
                }
            }
            // If there are usual M2M to do
            if( counterRealLowerCell ){
                kernel.M2M(pole->getRealCell(), realChild, inLevel);
            }
        }
    }


    /** The M2L should take into account if the current cell is adaptative or not.
      * Else we should work on the sub adaptative.
      * If it is composed or particles we have to perform M2P or P2P
      * Else we have to perform M2L or P2L.
      */
    void M2L(FAdaptativeCell<CellClass, ContainerClass>* const FRestrict local,
                const FAdaptativeCell<CellClass, ContainerClass>* distantNeighbors[343],
                const int /*size*/, const int inLevel)  override {
        // In case usual M2L can be done
        const CellClass* normalDistantNeighbors[343];
        int normalSize = 0;
        // The current adaptative cell
        FAdaptativeCell<CellClass, ContainerClass>* currentAdaptativeCell = nullptr;
        int currentAdaptativeLevel       = -1;
        // If the current adaptative cell is the current cell
        if(local->isAdaptative()){
            currentAdaptativeCell  = local;
            currentAdaptativeLevel = inLevel;
            // Then we may have some M2L to do
            memset(normalDistantNeighbors, 0, 343*sizeof(CellClass*));
        }
        else{
            // Else we are working with a lower cell
            currentAdaptativeCell = local->getSubAdaptativeCell();
            currentAdaptativeLevel= local->getSubAdaptativeLevel();
        }

        for(int idxNeigh = 0 ; idxNeigh < 343 ; ++idxNeigh){
            if(distantNeighbors[idxNeigh]){
                // If the current cell is adaptative and the neighbor too
                if(distantNeighbors[idxNeigh]->isAdaptative() && local->isAdaptative()){
                    if(distantNeighbors[idxNeigh]->hasDevelopment() && currentAdaptativeCell->hasDevelopment()){
                        // If both have development than we can use usual M2L
                        normalDistantNeighbors[idxNeigh] = distantNeighbors[idxNeigh]->getRealCell();
                        normalSize += 1;
                    }
                    else if(currentAdaptativeCell->hasDevelopment()){
                        // If only current cell has development the neighbor has particles
                        for(int idxLeafSrc = 0 ; idxLeafSrc < distantNeighbors[idxNeigh]->getNbSubLeaves() ; ++idxLeafSrc){
                            kernel.P2L(currentAdaptativeCell->getRealCell(), currentAdaptativeLevel, distantNeighbors[idxNeigh]->getSubLeaf(idxLeafSrc));
                        }
                    }
                    else if(distantNeighbors[idxNeigh]->hasDevelopment()){
                        // If only current cell has particles the neighbor has development
                        for(int idxLeafTgt = 0 ; idxLeafTgt < currentAdaptativeCell->getNbSubLeaves() ; ++idxLeafTgt){
                            kernel.M2P(distantNeighbors[idxNeigh]->getRealCell(), currentAdaptativeLevel, currentAdaptativeCell->getSubLeaf(idxLeafTgt));
                        }
                    }
                    else{
                        // If both have particles
                        for(int idxLeafTgt = 0 ; idxLeafTgt < currentAdaptativeCell->getNbSubLeaves() ; ++idxLeafTgt){
                            for(int idxLeafSrc = 0 ; idxLeafSrc < distantNeighbors[idxNeigh]->getNbSubLeaves() ; ++idxLeafSrc){
                                kernel.P2P(currentAdaptativeCell->getSubLeaf(idxLeafTgt), distantNeighbors[idxNeigh]->getSubLeaf(idxLeafSrc));
                            }
                        }
                    }
                }
                else{
                    const FAdaptativeCell<CellClass, ContainerClass>* lowerAdaptativeCell = distantNeighbors[idxNeigh];
                    int lowerAdaptativeLevel       = inLevel;
                    // If we need to look at lower level to find the adaptative cell
                    if(!distantNeighbors[idxNeigh]->isAdaptative()){
                        lowerAdaptativeCell  = distantNeighbors[idxNeigh]->getSubAdaptativeCell();
                        lowerAdaptativeLevel = distantNeighbors[idxNeigh]->getSubAdaptativeLevel();
                    }

                    if(lowerAdaptativeCell->hasDevelopment() && currentAdaptativeCell->hasDevelopment()){
                        // We are doing a M2L with distant interaction
                        kernel.M2L(currentAdaptativeCell->getRealCell(), currentAdaptativeLevel,
                            lowerAdaptativeCell->getRealCell(), lowerAdaptativeLevel);
                    }
                    else if(currentAdaptativeCell->hasDevelopment()){
                        // If only current cell has development the neighbor has particles
                        for(int idxLeafSrc = 0 ; idxLeafSrc < lowerAdaptativeCell->getNbSubLeaves() ; ++idxLeafSrc){
                            kernel.P2L(currentAdaptativeCell->getRealCell(), currentAdaptativeLevel, lowerAdaptativeCell->getSubLeaf(idxLeafSrc));
                        }
                    }
                    else if(lowerAdaptativeCell->hasDevelopment()){
                        // If only current cell has particles the neighbor has development
                        for(int idxLeafTgt = 0 ; idxLeafTgt < currentAdaptativeCell->getNbSubLeaves() ; ++idxLeafTgt){
                            kernel.M2P(lowerAdaptativeCell->getRealCell(), currentAdaptativeLevel, currentAdaptativeCell->getSubLeaf(idxLeafTgt));
                        }
                    }
                    else{
                        // If both have particles
                        for(int idxLeafTgt = 0 ; idxLeafTgt < currentAdaptativeCell->getNbSubLeaves() ; ++idxLeafTgt){
                            for(int idxLeafSrc = 0 ; idxLeafSrc < lowerAdaptativeCell->getNbSubLeaves() ; ++idxLeafSrc){
                                kernel.P2P(currentAdaptativeCell->getSubLeaf(idxLeafTgt), lowerAdaptativeCell->getSubLeaf(idxLeafSrc));
                            }
                        }
                    }
                }
            }
        }

        // If we work on the current cell and it has development
        if(normalSize){
            kernel.M2L(local->getRealCell(), normalDistantNeighbors, normalSize, inLevel);
        }
    }


    /** Nothing special */
    void finishedLevelM2L(const int level){
        kernel.finishedLevelM2L(level);
    }


    /** If the current cell is not adaptative we have nothing to do.
      * If it is adaptative we have to test the children.
      * If a child is adaptative than
      *     it has development and it is a normal L2L that should be performed
      *     or it has a list of leaves (particles) and we need L2P
      * Else it points to the lower adaptative cell and we need a L2L or L2P as in the previous case
      */
    void L2L(const FAdaptativeCell<CellClass, ContainerClass>* const FRestrict local,
             FAdaptativeCell<CellClass, ContainerClass>* FRestrict * const FRestrict child, const int inLevel)  override {
        // If there is something on this cell
        if(local->isAdaptative()){
            // We store the usual cell
            CellClass* realChild[8] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
            int counterRealChild    = 0;
            // For each child
            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
                if(child[idxChild]){
                    // If it is adaptative (it holds development or particles)
                    if(child[idxChild]->isAdaptative()){
                        if(child[idxChild]->hasDevelopment()){
                            // We need to perform a usual L2L on it
                            realChild[idxChild] = child[idxChild]->getRealCell();
                            counterRealChild   += 1;
                        }
                        else {
                            // We need to propagate on the particles
                            for(int idxLeafSrc = 0 ; idxLeafSrc < child[idxChild]->getNbSubLeaves() ; ++idxLeafSrc){
                                kernel.L2P(local->getRealCell(), inLevel, child[idxChild]->getSubLeaf(idxLeafSrc));
                            }
                        }
                    }
                    else{
                        // Get the lower adaptative cell
                        FAdaptativeCell<CellClass, ContainerClass>* lowerAdaptativeCell = child[idxChild]->getSubAdaptativeCell();
                        const int lowerAdaptativeLevel = child[idxChild]->getSubAdaptativeLevel();
                        if(lowerAdaptativeCell->hasDevelopment()){
                            // If it has a development we do a L2L with more than 1 level difference
                            kernel.L2L(local->getRealCell(), inLevel, lowerAdaptativeCell->getRealCell(), lowerAdaptativeLevel);
                        }
                        else{
                            // Else we propagate on the particles
                            for(int idxLeafSrc = 0 ; idxLeafSrc < child[idxChild]->getNbSubLeaves() ; ++idxLeafSrc){
                                kernel.L2P(local->getRealCell(), inLevel, child[idxChild]->getSubLeaf(idxLeafSrc));
                            }
                        }
                    }
                }
            }
            // Perform the usual L2L
            if(counterRealChild){
                kernel.L2L(local->getRealCell(), realChild, inLevel);
            }
        }
    }

    /** We do a Local to Particles only if the local (leaf) cell has some development */
    void L2P(const FAdaptativeCell<CellClass, ContainerClass>* const local, ContainerClass* const particles)  override {
        if(local->hasDevelopment()){
            kernel.L2P(local->getRealCell(), particles);
        }
    }

    /** This is a normal P2P */
    void P2P(const FTreeCoordinate& inLeafPosition,
             ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
             ContainerClass* const directNeighborsParticles[27], const int size)  override {
        kernel.P2P(inLeafPosition, targets, sources, directNeighborsParticles, size);
    }

    /** This is a normal P2P */
    void P2PRemote(const FTreeCoordinate& inLeafPosition,
                   ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
                   ContainerClass* const directNeighborsParticles[27], const int size) override {
        kernel.P2PRemote(inLeafPosition, targets, sources, directNeighborsParticles, size);
    }

};

#endif
