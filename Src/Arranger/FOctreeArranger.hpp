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
#ifndef FOCTREEARRANGER_HPP
#define FOCTREEARRANGER_HPP

#include "../Utils/FGlobal.hpp"
#include "../Utils/FPoint.hpp"
#include "../Containers/FVector.hpp"
#include "../Utils/FAssertable.hpp"


/** This class is an arranger, it move the particles that need
  * to be hosted in a different leaf
  * \n
  * For example, if a simulation has been executed and the position
  * of the particles have been changed, then it may be better
  * to move the particles in the tree instead of building a new
  * tree.
  */
template <class OctreeClass, class ContainerClass, class ParticleClass>
class FOctreeArranger : FAssertable {
    OctreeClass* const tree; //< The tree to work on

public:
    /** Basic constructor */
    FOctreeArranger(OctreeClass* const inTree) : tree(inTree) {
        fassert(tree, "Tree cannot be null", __LINE__ , __FILE__ );
    }

    /** Arrange */
    void rearrange(const bool isPeriodic = false){
        // This vector is to keep the moving particles
        FVector<ParticleClass> tomove;

        // For periodic
        const FReal boxWidth = tree->getBoxWidth();
        const FPoint min(tree->getBoxCenter(),-boxWidth/2);
        const FPoint max(tree->getBoxCenter(),boxWidth/2);

        { // iterate on the leafs and found particle to remove
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            do{
                const MortonIndex currentIndex = octreeIterator.getCurrentGlobalIndex();

                typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
                while( iter.hasNotFinished() ){
                    if(isPeriodic){
                        FPoint partPos = iter.data().getPosition();
                        while(partPos.getX() < min.getX()){
                            partPos.incX(boxWidth);
                        }
                        while(partPos.getX() >= max.getX()){
                            partPos.incX(-boxWidth);
                        }
                        while(partPos.getY() < min.getY()){
                            partPos.incY(boxWidth);
                        }
                        while(partPos.getY() >= max.getY()){
                            partPos.incY(-boxWidth);
                        }
                        while(partPos.getZ() < min.getZ()){
                            partPos.incZ(boxWidth);
                        }
                        while(partPos.getZ() >= max.getZ()){
                            partPos.incZ(-boxWidth);
                        }
                        iter.data().setPosition(partPos);
                    }
                    const MortonIndex particuleIndex = tree->getMortonFromPosition(iter.data().getPosition());
                    if(particuleIndex != currentIndex){
                        tomove.push(iter.data());
                        iter.remove();
                    }
                    else {
                        iter.gotoNext();
                    }
                }
            } while(octreeIterator.moveRight());
        }

        { // insert particles that moved
            for(int idxPart = 0 ; idxPart < tomove.getSize() ; ++idxPart){
                tree->insert(tomove[idxPart]);
            }
        }

        { // Remove empty leaves
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            bool workOnNext = true;
            do{
                // Empty leaf
                if( octreeIterator.getCurrentListTargets()->getSize() == 0 ){
                    const MortonIndex currentIndex = octreeIterator.getCurrentGlobalIndex();
                    workOnNext = octreeIterator.moveRight();
                    tree->removeLeaf( currentIndex );
                }
                // Not empty, just continue
                else {
                    workOnNext = octreeIterator.moveRight();
                }
            } while( workOnNext );
        }
    }

};

#endif // FOCTREEARRANGER_HPP
