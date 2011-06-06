#ifndef FTESTKERNELS_HPP
#define FTESTKERNELS_HPP
// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include "FAbstractKernels.hpp"
#include "../Containers/FList.hpp"
#include "../Containers/FOctree.hpp"

#include "../Utils/FTrace.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class AbstractKernels
* @brief
* Please read the license
*
* This kernels is a virtual kernels to validate that the fmm algorithm is
* correctly done on particles.
* It used FTestCell and FTestParticle
*/
template< class ParticleClass, class CellClass>
class FTestKernels : public FAbstractKernels<ParticleClass,CellClass> {
public:
    /** Default destructor */
    virtual ~FTestKernels(){
    }

    // Before upward
    void P2M(CellClass* const pole, const FList<ParticleClass>* const particles) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        // the pole represents all particles under
        pole->setDataUp(particles->getSize());
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    // During upward
    void M2M(CellClass* const FRestrict pole, const CellClass *const FRestrict *const FRestrict child, const int ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        // A parent represents the sum of the child
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child[idx]){
                pole->setDataUp(pole->getDataUp() + child[idx]->getDataUp());
            }
        }
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    // Before Downward
    void M2L(CellClass* const FRestrict pole, const CellClass* distantNeighbors[208], const FTreeCoordinate& , FTreeCoordinate [208], const int size, const int ) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        // The pole is impacted by what represent other poles
        for(int idx = 0 ; idx < size ; ++idx){
            pole->setDataDown(pole->getDataDown() + distantNeighbors[idx]->getDataUp());
        }
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    // During Downward
    void L2L(const CellClass*const FRestrict local, CellClass* FRestrict *const FRestrict child, const int) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        // Each child is impacted by the father
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child[idx]){
                child[idx]->setDataDown(local->getDataDown() + child[idx]->getDataDown());
            }
        }
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    // After Downward
    void L2P(const CellClass* const  local, FList<ParticleClass>*const particles){
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        // The particles is impacted by the parent cell
        typename FList<ParticleClass>::BasicIterator iter(*particles);
        while( iter.hasNotFinished() ){
            iter.data().setDataDown(iter.data().getDataDown() + local->getDataDown());
            iter.gotoNext();
        }
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

    // After Downward
    void P2P(FList<ParticleClass>* const FRestrict targets, const FList<ParticleClass>* const FRestrict sources,
             const FList<ParticleClass>* const directNeighborsParticles[26], const int size) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        // Each particles targeted is impacted by the particles sources
        long inc = sources->getSize();
        if(targets == sources){
            inc -= 1;
        }
        for(int idx = 0 ; idx < size ; ++idx){
            inc += directNeighborsParticles[idx]->getSize();
        }

        typename FList<ParticleClass>::BasicIterator iter(*targets);
        while( iter.hasNotFinished() ){
            iter.data().setDataDown(iter.data().getDataDown() + inc);
            iter.gotoNext();
        }
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }


    // After Downward
    void P2P(const MortonIndex ,
             FList<ParticleClass>* const FRestrict targets, const FList<ParticleClass>* const FRestrict sources,
             FList<ParticleClass>* const directNeighborsParticles[26], const MortonIndex [26], const int size) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        // Each particles targeted is impacted by the particles sources
        long inc = sources->getSize();
        if(targets == sources){
            inc -= 1;
        }
        for(int idx = 0 ; idx < size ; ++idx){
            inc += directNeighborsParticles[idx]->getSize();
        }

        typename FList<ParticleClass>::BasicIterator iter(*targets);
        while( iter.hasNotFinished() ){
            iter.data().setDataDown(iter.data().getDataDown() + inc);
            iter.gotoNext();
        }
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }
};


/** This function test the octree to be sure that the fmm algorithm
  * has worked completly.
  */
template< class ParticleClass, class CellClass, template<class ParticleClass> class LeafClass>
void ValidateFMMAlgo(FOctree<ParticleClass, CellClass, LeafClass>* const tree){
    std::cout << "Check Result\n";
    const int TreeHeight = tree->getHeight();
    int NbPart = 0;
    { // Check that each particle has been summed with all other
        typename FOctree<ParticleClass, CellClass, LeafClass>::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            if(octreeIterator.getCurrentCell()->getDataUp() != octreeIterator.getCurrentListSrc()->getSize() ){
                    std::cout << "Problem P2M : " << (octreeIterator.getCurrentCell()->getDataUp() - octreeIterator.getCurrentListSrc()->getSize()) << "\n";
            }
            NbPart += octreeIterator.getCurrentListSrc()->getSize();
        } while(octreeIterator.moveRight());
    }
    { // Ceck if there is number of NbPart summed at level 1
        typename FOctree<ParticleClass, CellClass, LeafClass>::Iterator octreeIterator(tree);
        octreeIterator.moveDown();
        long res = 0;
        do{
            res += octreeIterator.getCurrentCell()->getDataUp();
        } while(octreeIterator.moveRight());
        if(res != NbPart){
            std::cout << "Problem M2M at level 1 : " << res << "\n";
        }
    }
    { // Ceck if there is number of NbPart summed at level 1
        typename FOctree<ParticleClass, CellClass, LeafClass>::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        for(int idxLevel = TreeHeight - 1 ; idxLevel > 1 ; --idxLevel ){
            long res = 0;
            do{
                res += octreeIterator.getCurrentCell()->getDataUp();
            } while(octreeIterator.moveRight());
            if(res != NbPart){
                std::cout << "Problem M2M at level " << idxLevel << " : " << res << "\n";
            }
            octreeIterator.moveUp();
            octreeIterator.gotoLeft();
        }
    }
    { // Check that each particle has been summed with all other
        typename FOctree<ParticleClass, CellClass, LeafClass>::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename FList<ParticleClass>::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());

            while( iter.hasNotFinished() ){
                // If a particles has been impacted by less than NbPart - 1 (the current particle)
                // there is a problem
                if( (!isUsingTsm && iter.data().getDataDown() != NbPart - 1) ||
                    (isUsingTsm && iter.data().getDataDown() != NbPart) ){
                    std::cout << "Problem L2P + P2P : " << iter.data().getDataDown() << "\n";
                }
                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());
    }

    std::cout << "Done\n";
}



#endif //FTESTKERNELS_HPP

// [--LICENSE--]
