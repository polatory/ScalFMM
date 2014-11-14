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
#include "../Utils/FAssert.hpp"

#include "../Utils/FGlobalPeriodic.hpp"

/**
* This example show how to use the FOctreeArranger.
* @example testOctreeRearrange.cpp
*/

/**
* @brief This class is an arranger, it moves the particles that need to be hosted in a different leaf.
*
* For example, if a simulation has been executed and the position
* of the particles have been changed, then it may be better
* to move the particles in the tree instead of building a new
* tree.
*/
template <class OctreeClass, class ContainerClass, class ParticleClass, class ConverterClass >
class FOctreeArranger {
    OctreeClass* const tree; //< The tree to work on

public:
    /** Basic constructor */
    FOctreeArranger(OctreeClass* const inTree) : tree(inTree) {
        FAssertLF(tree, "Tree cannot be null" );
    }

    /** Arrange */
    void rearrange(const int isPeriodic = DirNone){
        // This vector is to keep the moving particles
        FVector<ParticleClass> tomove;

        // For periodic
        const FReal boxWidth = tree->getBoxWidth();
        const FPoint min(tree->getBoxCenter(),-boxWidth/2);
        const FPoint max(tree->getBoxCenter(),boxWidth/2);

        { // iterate on the leafs and found particle to remove
            FVector<int> indexesToExtract;

            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            do{
                const MortonIndex currentIndex = octreeIterator.getCurrentGlobalIndex();
                ContainerClass* particles = octreeIterator.getCurrentLeaf()->getSrc();
                for(int idxPart = 0 ; idxPart < particles->getNbParticles(); ++idxPart){
                    FPoint partPos( particles->getPositions()[0][idxPart],
                            particles->getPositions()[1][idxPart],
                            particles->getPositions()[2][idxPart] );
                    // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    if( TestPeriodicCondition(isPeriodic, DirPlusX) ){
                        while(partPos.getX() >= max.getX()){
                            partPos.incX(-boxWidth);
                        }
                    }
                    else if(partPos.getX() >= max.getX()){
                        printf("Error, particle out of Box in +X, index %lld\n", currentIndex);
                        printf("Application is exiting...\n");
                    }
                    if( TestPeriodicCondition(isPeriodic, DirMinusX) ){
                        while(partPos.getX() < min.getX()){
                            partPos.incX(boxWidth);
                        }
                    }
                    else if(partPos.getX() < min.getX()){
                        printf("Error, particle out of Box in -X, index %lld\n", currentIndex);
                        printf("Application is exiting...\n");
                    }
                    // YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
                    if( TestPeriodicCondition(isPeriodic, DirPlusY) ){
                        while(partPos.getY() >= max.getY()){
                            partPos.incY(-boxWidth);
                        }
                    }
                    else if(partPos.getY() >= max.getY()){
                        printf("Error, particle out of Box in +Y, index %lld\n", currentIndex);
                        printf("Application is exiting...\n");
                    }
                    if( TestPeriodicCondition(isPeriodic, DirMinusY) ){
                        while(partPos.getY() < min.getY()){
                            partPos.incY(boxWidth);
                        }
                    }
                    else if(partPos.getY() < min.getY()){
                        printf("Error, particle out of Box in -Y, index %lld\n", currentIndex);
                        printf("Application is exiting...\n");
                    }
                    // ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
                    if( TestPeriodicCondition(isPeriodic, DirPlusX) ){
                        while(partPos.getZ() >= max.getZ()){
                            partPos.incZ(-boxWidth);
                        }
                    }
                    else if(partPos.getZ() >= max.getZ()){
                        printf("Error, particle out of Box in +Z, index %lld\n", currentIndex);
                        printf("Application is exiting...\n");
                    }
                    if( TestPeriodicCondition(isPeriodic, DirMinusX) ){
                        while(partPos.getZ() < min.getZ()){
                            partPos.incZ(boxWidth);
                        }
                    }
                    else if(partPos.getZ() < min.getZ()){
                        printf("Error, particle out of Box in -Z, index %lld\n", currentIndex);
                        printf("Application is exiting...\n");
                    }
                    // set pos
                    particles->getWPositions()[0][idxPart] = partPos.getX();
                    particles->getWPositions()[1][idxPart] = partPos.getY();
                    particles->getWPositions()[2][idxPart] = partPos.getZ();

                    const MortonIndex particuleIndex = tree->getMortonFromPosition(partPos);
                    if(particuleIndex != currentIndex){
                        tomove.push(ConverterClass::GetParticleAndRemove(particles,idxPart));
                    }
                }
            } while(octreeIterator.moveRight());
        }

        { // insert particles that moved
            for(int idxPart = 0 ; idxPart < tomove.getSize() ; ++idxPart){
                ConverterClass::Insert( tree , tomove[idxPart]);
            }
        }

        { // Remove empty leaves
            typename OctreeClass::Iterator octreeIterator(tree);
            octreeIterator.gotoBottomLeft();
            bool workOnNext = true;
            do{
                // Empty leaf
                if( octreeIterator.getCurrentListTargets()->getNbParticles() == 0 ){
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
