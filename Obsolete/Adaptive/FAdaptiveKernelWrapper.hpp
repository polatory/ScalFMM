#ifndef FAdaptiveKERNELWRAPPER_HPP
#define FAdaptiveKERNELWRAPPER_HPP
// Keep in private GIT
// @SCALFMM_PRIVATE
#include "Components/FAbstractKernels.hpp"
#include "Containers/FVector.hpp"

#include "FAdaptiveCell.hpp"



/**
 * This class is a wrapper to use the usual algorithm but in a Adaptive way.
 * The real computational kernel should be given in template.
 * It must propose the usual FMM Kernel operators and also the FAbstractKernels operators.
 *
 */
template <class RealAdaptiveKernel, class CellClass, class ContainerClass >
class FAdaptiveKernelWrapper : public FAbstractKernels<FAdaptiveCell<CellClass, ContainerClass>, ContainerClass>{
protected:
	RealAdaptiveKernel kernel;

public:

	template<typename... KernelParams>
	FAdaptiveKernelWrapper(KernelParams... parameters) : kernel(parameters...){
	}

	RealAdaptiveKernel& getKernel(){
		return kernel;
	}

	const RealAdaptiveKernel& getKernel() const{
		return kernel;
	}

	/** P2M is performed only if the kernel says that is it better than waiting for future development
	 * A leaf cell contains the particles container if the development is not made,
	 * else it has a multipole component.
	 */
	void P2M(FAdaptiveCell<CellClass, ContainerClass>* const pole, const ContainerClass* const particles)  override {
		// Leaf is always Adaptive
		pole->setAdaptive(true);
        pole->resetSubLeaves();

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
	 * If there is one child, we are not Adaptive and keep the lower adaptive cell as target
	 * If at least one child has a development (at any level) the current cell store the development of all the children
	 * using M2M or P2M from the leaves.
	 * If no children have developments we need to test if it is better to merge them in a development
	 * or to store all the particles in order to perform the P2M later.
	 * Finally if all the children have developments at level just under we can perform normal M2M
	 */
	void M2M(FAdaptiveCell<CellClass, ContainerClass>* const FRestrict pole,
			const FAdaptiveCell<CellClass, ContainerClass>*const FRestrict *const FRestrict child, const int inLevel)  override {
		int nbChild       = 0;
		int lastChild     = 0;
		bool onlyParticlesCells         = true;
		FVector<const ContainerClass*> subLeaves;
        pole->resetSubLeaves();
		// Test all the children
		for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
			if(child[idxChild]){
				nbChild  += 1;
				lastChild = idxChild;
				// We aggregate all the leaves from the current child or its Adaptive cell
				if(child[idxChild]->isAdaptive() && !child[idxChild]->hasDevelopment()){
					subLeaves.memocopy( child[idxChild]->getSubLeaves(), child[idxChild]->getNbSubLeaves());
				}
				else if(!child[idxChild]->isAdaptive() && !child[idxChild]->getSubAdaptiveCell()->hasDevelopment()){
					subLeaves.memocopy( child[idxChild]->getSubAdaptiveCell()->getSubLeaves(),
							child[idxChild]->getSubAdaptiveCell()->getNbSubLeaves());
				}
				else{
					// If a child is made of development
					onlyParticlesCells = false;
				}
			}
		}
        // We need to aggregate if there are only particles and if the kernel says so
        const bool continueToAgregate = (onlyParticlesCells && (kernel.preferP2M(inLevel, subLeaves.data(), int(subLeaves.getSize())) == false));
		if(nbChild == 1){
			// One child means that the cell is not Adaptive
			pole->setAdaptive(false);
			pole->setHaveDevelopment(false);
            if(child[lastChild]->isAdaptive()){
				pole->setSubAdaptiveCell(child[lastChild], inLevel + 1);
			}
			else{
				pole->setSubAdaptiveCell(child[lastChild]->getSubAdaptiveCell(),
						child[lastChild]->getSubAdaptiveLevel());
			}
		}
		else if(onlyParticlesCells && continueToAgregate){ // nbChild > 1
			// There are only particles and no development to do
            pole->setAdaptive(true);
			pole->setHaveDevelopment(false);
			for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
				if(child[idxChild]){
                    if(child[idxChild]->isAdaptive()){
						pole->addSubLeaves(child[idxChild]->getSubLeaves(), child[idxChild]->getNbSubLeaves());
					}
					else{
						const FAdaptiveCell<CellClass, ContainerClass>* lowerAdaptiveCell = child[idxChild]->getSubAdaptiveCell();
						pole->addSubLeaves(lowerAdaptiveCell->getSubLeaves(), lowerAdaptiveCell->getNbSubLeaves());
					}
				}
			}
		}
		else{
			// There development to do from developments or particles
			pole->setAdaptive(true);
			pole->setHaveDevelopment(true);

			const CellClass* realChild[8] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
			int counterRealLowerCell      = 0;
			// Test each child
			for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
				if(child[idxChild]){
                    if(child[idxChild]->isAdaptive()){
						if(child[idxChild]->hasDevelopment()){
							// If it is Adaptive and has development than we compute is using usual M2M
							realChild[idxChild] = child[idxChild]->getRealCell();
							counterRealLowerCell += 1;
						}
						else{
							// If it is Adaptive and has not development than we compute is using P2M
							for(int idxLeaf = 0 ; idxLeaf < child[idxChild]->getNbSubLeaves() ; ++idxLeaf){
								kernel.P2M(pole->getRealCell(), inLevel, child[idxChild]->getSubLeaf(idxLeaf));
                            }
						}
					}
					else{
						// Else we need the Adaptive cell
						const FAdaptiveCell<CellClass, ContainerClass>* lowerAdaptiveCell = child[idxChild]->getSubAdaptiveCell();
						const int lowerAdaptiveLevel = child[idxChild]->getSubAdaptiveLevel();
						if(lowerAdaptiveCell->hasDevelopment()){
							// If it has development we perform a M2M
							kernel.M2M(pole->getRealCell(), inLevel, lowerAdaptiveCell->getRealCell(),
									lowerAdaptiveLevel);
                        }
						else{
							// Else we perform P2M
							for(int idxLeaf = 0 ; idxLeaf < lowerAdaptiveCell->getNbSubLeaves() ; ++idxLeaf){
								kernel.P2M(pole->getRealCell(), inLevel, lowerAdaptiveCell->getSubLeaf(idxLeaf));
                            }
						}
					}
				}
			}
			// If there are usual M2M to do
			if( counterRealLowerCell ){
                kernel.M2M(pole->getRealCell(), realChild, inLevel);
			}
        } // End of 	if(nbChild == 1)
	}


	/** The M2L should take into account if the current cell is Adaptive or not.
	 * Else we should work on the sub Adaptive.
	 * If it is composed or particles we have to perform M2P or P2P
	 * Else we have to perform M2L or P2L.
	 */
	void M2L(FAdaptiveCell<CellClass, ContainerClass>* const FRestrict local,
            const FAdaptiveCell<CellClass, ContainerClass>* distantNeighbors[],
             const int positions[],
            const int size, const int inLevel)  override {
		// In case usual M2L can be done
        const CellClass* normalDistantNeighbors[189];
        int normalPositions[189];
		int normalSize = 0;
		// The current Adaptive cell
		FAdaptiveCell<CellClass, ContainerClass>* currentAdaptiveCell = nullptr;
		int currentAdaptiveLevel       = -1;
		// If the current Adaptive cell is the current cell
        if(local->isAdaptive()){
			currentAdaptiveCell  = local;
			currentAdaptiveLevel = inLevel;
			// Then we may have some M2L to do
            memset(normalDistantNeighbors, 0, 189*sizeof(CellClass*));
		}
		else{
			// Else we are working with a lower cell
			currentAdaptiveCell = local->getSubAdaptiveCell();
			currentAdaptiveLevel= local->getSubAdaptiveLevel();
		}

        for(int idxNeighExisting = 0 ; idxNeighExisting < size ; ++idxNeighExisting){
            const int idxNeigh = positions[idxNeighExisting];
            // If the current cell is Adaptive and the neighbor too
            if(distantNeighbors[idxNeighExisting]->isAdaptive() && local->isAdaptive()){
                if(distantNeighbors[idxNeighExisting]->hasDevelopment() && currentAdaptiveCell->hasDevelopment()){
                    // If both have development than we can use usual M2L
                    normalDistantNeighbors[normalSize] = distantNeighbors[idxNeighExisting]->getRealCell();
                    normalPositions[normalSize] = idxNeigh;
                    normalSize += 1;
                }
                else if(currentAdaptiveCell->hasDevelopment()){
                    // If only current cell has development the neighbor has particles
                    for(int idxLeafSrc = 0 ; idxLeafSrc < distantNeighbors[idxNeighExisting]->getNbSubLeaves() ; ++idxLeafSrc){
                        kernel.P2L(currentAdaptiveCell->getRealCell(), currentAdaptiveLevel, distantNeighbors[idxNeighExisting]->getSubLeaf(idxLeafSrc));
                    }
                }
                else if(distantNeighbors[idxNeighExisting]->hasDevelopment()){
                    // If only current cell has particles the neighbor has development
                    for(int idxLeafTgt = 0 ; idxLeafTgt < currentAdaptiveCell->getNbSubLeaves() ; ++idxLeafTgt){
                        kernel.M2P(distantNeighbors[idxNeighExisting]->getRealCell(), currentAdaptiveLevel, currentAdaptiveCell->getSubLeaf(idxLeafTgt));
                    }
                }
                else{
                    // If both have particles
                    for(int idxLeafTgt = 0 ; idxLeafTgt < currentAdaptiveCell->getNbSubLeaves() ; ++idxLeafTgt){
                        for(int idxLeafSrc = 0 ; idxLeafSrc < distantNeighbors[idxNeighExisting]->getNbSubLeaves() ; ++idxLeafSrc){
                            kernel.P2P(currentAdaptiveCell->getSubLeaf(idxLeafTgt), distantNeighbors[idxNeighExisting]->getSubLeaf(idxLeafSrc));
                        }
                    }
                }
            }
            else{
                const FAdaptiveCell<CellClass, ContainerClass>* lowerAdaptiveCell = distantNeighbors[idxNeighExisting];
                int lowerAdaptiveLevel       = inLevel;
                // If we need to look at lower level to find the Adaptive cell
                if(!distantNeighbors[idxNeighExisting]->isAdaptive()){
                    lowerAdaptiveCell  = distantNeighbors[idxNeighExisting]->getSubAdaptiveCell();
                    lowerAdaptiveLevel = distantNeighbors[idxNeighExisting]->getSubAdaptiveLevel();
                }

                if(lowerAdaptiveCell->hasDevelopment() && currentAdaptiveCell->hasDevelopment()){
                    // We are doing a M2L with distant interaction
                    kernel.M2L(currentAdaptiveCell->getRealCell(), currentAdaptiveLevel,
                            lowerAdaptiveCell->getRealCell(), lowerAdaptiveLevel);
                }
                else if(currentAdaptiveCell->hasDevelopment()){
                    // If only current cell has development the neighbor has particles
                    for(int idxLeafSrc = 0 ; idxLeafSrc < lowerAdaptiveCell->getNbSubLeaves() ; ++idxLeafSrc){
                        kernel.P2L(currentAdaptiveCell->getRealCell(), currentAdaptiveLevel, lowerAdaptiveCell->getSubLeaf(idxLeafSrc));
                    }
                }
                else if(lowerAdaptiveCell->hasDevelopment()){
                    // If only current cell has particles the neighbor has development
                    for(int idxLeafTgt = 0 ; idxLeafTgt < currentAdaptiveCell->getNbSubLeaves() ; ++idxLeafTgt){
                        kernel.M2P(lowerAdaptiveCell->getRealCell(), lowerAdaptiveLevel, currentAdaptiveCell->getSubLeaf(idxLeafTgt));
                    }
                }
                else{
                    // If both have particles
                    for(int idxLeafTgt = 0 ; idxLeafTgt < currentAdaptiveCell->getNbSubLeaves() ; ++idxLeafTgt){
                        for(int idxLeafSrc = 0 ; idxLeafSrc < lowerAdaptiveCell->getNbSubLeaves() ; ++idxLeafSrc){
                            kernel.P2P(currentAdaptiveCell->getSubLeaf(idxLeafTgt), lowerAdaptiveCell->getSubLeaf(idxLeafSrc));
                        }
                    }
                }
			}
		}

		// If we work on the current cell and it has development
		if(normalSize){
            kernel.M2L(local->getRealCell(), normalDistantNeighbors, normalPositions, normalSize, inLevel);
        }
	}


	/** Nothing special */
	void finishedLevelM2L(const int level){
		kernel.finishedLevelM2L(level);
	}


	/** If the current cell is not Adaptive we have nothing to do.
	 * If it is Adaptive we have to test the children.
	 * If a child is Adaptive than
	 *     it has development and it is a normal L2L that should be performed
	 *     or it has a list of leaves (particles) and we need L2P
	 * Else it points to the lower Adaptive cell and we need a L2L or L2P as in the previous case
	 */
	void L2L(const FAdaptiveCell<CellClass, ContainerClass>* const FRestrict local,
			FAdaptiveCell<CellClass, ContainerClass>* FRestrict * const FRestrict child, const int inLevel)  override {
		// If there is something on this cell
        if(local->isAdaptive() && local->hasDevelopment()){
            // We store the usual cell
			CellClass* realChild[8] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
			int counterRealChild    = 0;
			// For each child
			for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
				if(child[idxChild]){
					// If it is Adaptive (it holds development or particles)
					if(child[idxChild]->isAdaptive()){
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
						// Get the lower Adaptive cell
						FAdaptiveCell<CellClass, ContainerClass>* lowerAdaptiveCell = child[idxChild]->getSubAdaptiveCell();
						const int lowerAdaptiveLevel = child[idxChild]->getSubAdaptiveLevel();
						if(lowerAdaptiveCell->hasDevelopment()){
							// If it has a development we do a L2L with more than 1 level difference
                            kernel.L2L(local->getRealCell(), inLevel, lowerAdaptiveCell->getRealCell(), lowerAdaptiveLevel);
						}
						else{
                            // Else we propagate on the particles
							for(int idxLeafSrc = 0 ; idxLeafSrc < lowerAdaptiveCell->getNbSubLeaves() ; ++idxLeafSrc){
                                kernel.L2P(local->getRealCell(), inLevel, lowerAdaptiveCell->getSubLeaf(idxLeafSrc));
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
	void L2P(const FAdaptiveCell<CellClass, ContainerClass>* const local, ContainerClass* const particles)  override {
        if(local->hasDevelopment()){
			kernel.L2P(local->getRealCell(), particles);
        }
	}

	/** This is a normal P2P */
	void P2P(const FTreeCoordinate& inLeafPosition,
			ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
            ContainerClass* const directNeighborsParticles[], const int positions[], const int size)  override {
        kernel.P2P(inLeafPosition, targets, sources, directNeighborsParticles, positions, size);
	}

    /** This is a normal P2P */
    void P2POuter(const FTreeCoordinate& inLeafPosition,
            ContainerClass* const FRestrict targets,
            ContainerClass* const directNeighborsParticles[], const int positions[], const int size)  override {
        kernel.P2POuter(inLeafPosition, targets, directNeighborsParticles, positions, size);
    }

	/** This is a normal P2P */
	void P2PRemote(const FTreeCoordinate& inLeafPosition,
			ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
            const ContainerClass* const directNeighborsParticles[], const int positions[], const int size) override {
        kernel.P2PRemote(inLeafPosition, targets, sources, directNeighborsParticles, positions, size);
	}

};

#endif
