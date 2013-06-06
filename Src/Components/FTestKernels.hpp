// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FTESTKERNELS_HPP
#define FTESTKERNELS_HPP


#include <iostream>

#include "FAbstractKernels.hpp"
#include "../Containers/FOctree.hpp"
#include "../Utils/FGlobal.hpp"
#include "../Utils/FTrace.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class AbstractKernels
* @brief
* Please read the license
*
* This kernels is a virtual kernels to validate that the fmm algorithm is
* correctly done on particles.
* It should use FTestCell and FTestParticle.
* A the end of a computation, the particles should host then number of particles
* in the simulation (-1).
*/
template< class CellClass, class ContainerClass>
class FTestKernels  : public FAbstractKernels<CellClass,ContainerClass> {
public:
    /** Default destructor */
    virtual ~FTestKernels(){
    }

    /** Before upward */
    void P2M(CellClass* const pole, const ContainerClass* const particles) {
        // the pole represents all particles under
        pole->setDataUp(pole->getDataUp() + particles->getNbParticles());
    }

    /** During upward */
    void M2M(CellClass* const FRestrict pole, const CellClass *const FRestrict *const FRestrict child, const int /*level*/) {
        // A parent represents the sum of the child
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child[idx]){
                pole->setDataUp(pole->getDataUp() + child[idx]->getDataUp());
            }
        }
    }

    /** Before Downward */
    void M2L(CellClass* const FRestrict pole, const CellClass* distantNeighbors[343], const int /*size*/, const int /*level*/) {
        // The pole is impacted by what represent other poles
        for(int idx = 0 ; idx < 343 ; ++idx){
            if(distantNeighbors[idx]){
                pole->setDataDown(pole->getDataDown() + distantNeighbors[idx]->getDataUp());
            }
        }
    }

    /** During Downward */
    void L2L(const CellClass*const FRestrict local, CellClass* FRestrict *const FRestrict child, const int /*level*/) {
        // Each child is impacted by the father
        for(int idx = 0 ; idx < 8 ; ++idx){
            if(child[idx]){
                child[idx]->setDataDown(local->getDataDown() + child[idx]->getDataDown());
            }
        }

    }

    /** After Downward */
    void L2P(const CellClass* const  local, ContainerClass*const particles){
        // The particles is impacted by the parent cell      
        long long int*const particlesAttributes = particles->getDataDown();
        for(int idxPart = 0 ; idxPart < particles->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += local->getDataDown();
        }
    }


    /** After Downward */
    void P2P(const FTreeCoordinate& ,
                 ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
                 ContainerClass* const directNeighborsParticles[27], const int ){
        // Each particles targeted is impacted by the particles sources
        long long int inc = sources->getNbParticles();
        if(targets == sources){
            inc -= 1;
        }
        for(int idx = 0 ; idx < 27 ; ++idx){
            if( directNeighborsParticles[idx] ){
                inc += directNeighborsParticles[idx]->getNbParticles();
            }
        }

        long long int*const particlesAttributes = targets->getDataDown();
        for(int idxPart = 0 ; idxPart < targets->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += inc;
        }
    }

    /** After Downward */
    void P2PRemote(const FTreeCoordinate& ,
                 ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
                 ContainerClass* const directNeighborsParticles[27], const int ){
        // Each particles targeted is impacted by the particles sources
        long long int inc = 0;
        for(int idx = 0 ; idx < 27 ; ++idx){
            if( directNeighborsParticles[idx] ){
                inc += directNeighborsParticles[idx]->getNbParticles();
            }
        }

        long long int*const particlesAttributes = targets->getDataDown();
        for(int idxPart = 0 ; idxPart < targets->getNbParticles() ; ++idxPart){
            particlesAttributes[idxPart] += inc;
        }
    }
};


/** This function test the octree to be sure that the fmm algorithm
  * has worked completly.
  */
template< class OctreeClass, class CellClass, class ContainerClass, class LeafClass>
void ValidateFMMAlgo(OctreeClass* const tree){
    std::cout << "Check Result\n";
    const int TreeHeight = tree->getHeight();
    long long int NbPart = 0;
    { // Check that each particle has been summed with all other
        tree->forEachCellLeaf([&](CellClass* cell, LeafClass* leaf){
            if(cell->getDataUp() != leaf->getSrc()->getNbParticles() ){
                    std::cout << "Problem P2M : " << cell->getDataUp() <<
                                 " (should be " << leaf->getSrc()->getNbParticles() << ")\n";
            }
            NbPart += leaf->getSrc()->getNbParticles();
        });
    }
    { // Ceck if there is number of NbPart summed at level 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.moveDown();
        long long int res = 0;
        do{
            res += octreeIterator.getCurrentCell()->getDataUp();
        } while(octreeIterator.moveRight());
        if(res != NbPart){
            std::cout << "Problem M2M at level 1 : " << res << "\n";
        }
    }
    { // Ceck if there is number of NbPart summed at level 1
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        for(int idxLevel = TreeHeight - 1 ; idxLevel > 1 ; --idxLevel ){
            long long int res = 0;
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
        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());
            const long long int* dataDown = octreeIterator.getCurrentListTargets()->getDataDown();
            for(int idxPart = 0 ; idxPart < octreeIterator.getCurrentListTargets()->getNbParticles() ; ++idxPart){
                if( (!isUsingTsm && dataDown[idxPart] != NbPart - 1) ||
                    (isUsingTsm && dataDown[idxPart] != NbPart) ){
                    std::cout << "Problem L2P + P2P : " << dataDown[idxPart] <<
                                 "(" << octreeIterator.getCurrentGlobalIndex() << ")\n";
                }
            }
        } while(octreeIterator.moveRight());
    }

    std::cout << "Done\n";
}



#endif //FTESTKERNELS_HPP


