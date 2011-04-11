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
* correctly done on particules.
* It used FTestCell and FTestParticule
*/
template< class ParticuleClass, class CellClass, int TreeHeight>
class FTestKernels : public FAbstractKernels<ParticuleClass,CellClass,TreeHeight> {
public:
    /** Default destructor */
    virtual ~FTestKernels(){
    }

    // Before upward
    void P2M(CellClass* const pole, const FList<ParticuleClass*>* const particules) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        // the pole represents all particules under
        pole->setDataUp(particules->getSize());
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
    void M2L(CellClass* const FRestrict pole, const CellClass*const FRestrict *const distantNeighbors, const int size, const int ) {
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
    void L2P(const CellClass* const  local, FList<ParticuleClass*>*const particules){
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        // The particules is impacted by the parent cell
        typename FList<ParticuleClass*>::BasicIterator iter(*particules);
        while( iter.isValide() ){
            iter.value()->setDataDown(iter.value()->getDataDown() + local->getDataDown());
            iter.progress();
        }
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }
    // After Downward
    void P2P(FList<ParticuleClass*>* const FRestrict targets, const FList<ParticuleClass*>* const FRestrict sources,
             const FList<ParticuleClass*>* FRestrict const* FRestrict directNeighbors, const int size) {
        FTRACE( FTrace::Controller.enterFunction(FTrace::KERNELS, __FUNCTION__ , __FILE__ , __LINE__) );
        // Each particules targeted is impacted by the particules sources
        long inc = sources->getSize();
        if(targets == sources){
            inc -= 1;
        }
        for(int idx = 0 ; idx < size ; ++idx){
            inc += directNeighbors[idx]->getSize();
        }

        typename FList<ParticuleClass*>::BasicIterator iter(*targets);
        while( iter.isValide() ){
            iter.value()->setDataDown(iter.value()->getDataDown() + inc);
            iter.progress();
        }
        FTRACE( FTrace::Controller.leaveFunction(FTrace::KERNELS) );
    }

};


/** This function test the octree to be sure that the fmm algorithm
  * has worked completly.
  */
template< class ParticuleClass, class CellClass, template<class ParticuleClass> class LeafClass, int TreeHeight, int SizeSubLevels>
void ValidateFMMAlgo(FOctree<ParticuleClass, CellClass, LeafClass, TreeHeight , SizeSubLevels>* const tree){
    std::cout << "Check Result\n";
    int NbPart = 0;
    { // Check that each particule has been summed with all other
        typename FOctree<ParticuleClass, CellClass, LeafClass, TreeHeight, SizeSubLevels>::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            if(octreeIterator.getCurrentCell()->getDataUp() != octreeIterator.getCurrentListSources()->getSize() ){
                    std::cout << "Problem P2M : " << (octreeIterator.getCurrentCell()->getDataUp() - octreeIterator.getCurrentListSources()->getSize()) << "\n";
            }
            NbPart += octreeIterator.getCurrentListSources()->getSize();
        } while(octreeIterator.moveRight());
    }
    { // Ceck if there is number of NbPart summed at level 1
        typename FOctree<ParticuleClass, CellClass, LeafClass, TreeHeight, SizeSubLevels>::Iterator octreeIterator(tree);
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
        typename FOctree<ParticuleClass, CellClass, LeafClass, TreeHeight, SizeSubLevels>::Iterator octreeIterator(tree);
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
    { // Check that each particule has been summed with all other
        typename FOctree<ParticuleClass, CellClass, LeafClass, TreeHeight, SizeSubLevels>::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename FList<ParticuleClass*>::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            const bool isUsingToS = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSources());

            while( iter.isValide() ){
                // If a particules has been impacted by less than NbPart - 1 (the current particule)
                // there is a problem
                if( (!isUsingToS && iter.value()->getDataDown() != NbPart - 1) ||
                    (isUsingToS && iter.value()->getDataDown() != NbPart) ){
                    std::cout << "Problem L2P + P2P : " << iter.value()->getDataDown() << "\n";
                }
                iter.progress();
            }
        } while(octreeIterator.moveRight());
    }

    std::cout << "Done\n";
}



#endif //FTESTKERNELS_HPP

// [--LICENSE--]
