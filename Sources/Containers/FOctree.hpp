#ifndef FOCTREE_HPP
#define FOCTREE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FSubOctree.hpp"

#include "../Utils/F3DPosition.hpp"
#include "../Utils/FMath.hpp"
#include "FTreeCoordinate.hpp"

#include "../Utils/FAssertable.hpp"

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FOctree
 * Please read the license
 *
 * This class is an octree container.
 *
 * Please refere to testOctree.cpp to see an example.
 * <code>
 * // can be used as : <br>
 * FOctree<TestParticule, TestCell, 10, 3> tree(1.0,F3DPosition(0.5,0.5,0.5));
 * </code>
 *
 * Particules and cells has to respect the Abstract class definition.
 * Particule can extend {FExtendPosition}
 * Cell can extend {FExtendPosition,FExtendMortonIndex}
 */
template< class ParticuleClass, class CellClass, int OctreeHeight, int SubtreeHeight = 3>
class FOctree {
        const int height;		//< tree height
	const int leafIndex;		//< index of leaf int array

        const double boxWidth;          //< the space system width
	const F3DPosition boxCenter;	//< the space system center
	const F3DPosition boxCorner;	//< the space system corner (used to compute morton index)

	double boxWidthAtLevel[OctreeHeight];		//< to store the width of each boxs at all levels

        FSubOctree< ParticuleClass, CellClass > root;//< root suboctree

	/** Forbiden copy operator */
	FOctree& operator=(const FOctree&) {
		return *this;
	}

	/** Forbiden copy constructor */
	FOctree(const FOctree&) {
	}

        /**
        * Get morton index from a position for the leaf leavel
        * @param inPosition position to compute
        * @return the morton index
        */
        MortonIndex getLeafMortonFromPosition(const F3DPosition& inPosition) const {
                // box coordinate to host the particule
                FTreeCoordinate host;
                // position has to be relative to corner not center
                host.setX( getTreeCoordinate( inPosition.getX() - this->boxCorner.getX() ));
                host.setY( getTreeCoordinate( inPosition.getY() - this->boxCorner.getY() ));
                host.setZ( getTreeCoordinate( inPosition.getZ() - this->boxCorner.getZ() ));
                return host.getMortonIndex(leafIndex);
        }

        /**
        * Get the box number from a position
        * at a position POS with a leaf level box width of WI, the position is RELATIVE_TO_CORNER(POS)/WI
        * @param inRelativePosition a position from the corner of the box
        * @return the box num at the leaf level that contains inRelativePosition
        */
        long getTreeCoordinate(const double inRelativePosition) const {
                const double indexDouble = inRelativePosition / this->boxWidthAtLevel[this->leafIndex];
                const long index = FMath::dfloor(indexDouble);
                if( index && FMath::LookEqual(inRelativePosition, this->boxWidthAtLevel[this->leafIndex] * index ) ){
                        return index - 1;
                }
                return index;
        }

public:
	/**
	* Constructor	
	* @param inBoxWidth box width for this simulation
	* @param inBoxCenter box center for this simulation
	*/
	FOctree(const double inBoxWidth, const F3DPosition& inBoxCenter)
                        : boxWidth(inBoxWidth) , boxCenter(inBoxCenter), boxCorner(inBoxCenter - (inBoxWidth/2)),
                        height(OctreeHeight) , leafIndex(OctreeHeight-1),
                        root(0, 0, SubtreeHeight, 1){
		double tempWidth = this->boxWidth;
                // pre compute box width for each level
                for(int indexLevel = 0; indexLevel < this->height; ++indexLevel ){
                        this->boxWidthAtLevel[indexLevel] = tempWidth;
			tempWidth /= 2.0;
		}

	}

	/** Desctructor */
        virtual ~FOctree() {
	}

	/**
	* Insert a particule on the tree
	* algorithm is :
	* Compute morton index for the particule
	* ask node to insert this particule
	* @param inParticule the particule to insert (must inherite from FAbstractParticule)
	*/
	void insert(ParticuleClass* const inParticule){
                const MortonIndex particuleIndex = getLeafMortonFromPosition( inParticule->getPosition() );
                root.insert( particuleIndex, inParticule, this->height, this->boxWidthAtLevel);
	}


        /**
          * The class works on suboctree. Most of the resources needed
          * are avaiblable by using FAbstractSubOctree. But when accessing
          * to the leaf we have to use FSubOctree or FSubOctreeWithLeafs
          * depending if we are working on the bottom of the tree.
          */
        union SubOctreeTypes {
            FAbstractSubOctree<ParticuleClass,CellClass>* tree;     //< Usual pointer to work
            FSubOctree<ParticuleClass,CellClass>* middleTree;       //< To access to sub-octree under
            FSubOctreeWithLeafs<ParticuleClass,CellClass>* leafTree;//< To access to particules lists
        };

        /**
          * The class works on suboctree. Most of the resources needed
          * are avaiblable by using FAbstractSubOctree. But when accessing
          * to the leaf we have to use FSubOctree or FSubOctreeWithLeafs
          * depending if we are working on the bottom of the tree.
          */
        union SubOctreeTypesConst {
            const FAbstractSubOctree<ParticuleClass,CellClass>* tree;     //< Usual pointer to work
            const FSubOctree<ParticuleClass,CellClass>* middleTree;       //< To access to sub-octree under
            const FSubOctreeWithLeafs<ParticuleClass,CellClass>* leafTree;//< To access to particules lists
        };

        /**
          * This has to be used to iterate on an octree
          * It simply stores an pointer on a suboctree and moves to right/left/up/down.
          * Please refere to testOctreeIter file to see how it works.
          *
          * <code>
          * FOctree<TestParticule, TestCell, NbLevels, NbSubLevels>::Iterator octreeIterator(&tree); <br>
          * octreeIterator.gotoBottomLeft(); <br>
          * for(int idx = 0 ; idx < NbLevels - 1; ++idx ){ <br>
          *     do{ <br>
          *         // ...  <br>
          *     } while(octreeIterator.moveRight()); <br>
          *     octreeIterator.moveUp(); <br>
          *     octreeIterator.gotoLeft(); <br>
          * } <br>
          * <code>
          * Remark :
          * It uses the left right limit on each suboctree and their morton index.
          * Please have a look to the move functions to understand how the system is working.
          */
        class Iterator : protected FAssertable {
            SubOctreeTypes current; //< Current suboctree

            int currentLocalLevel;  //< Current level in the current suboctree
            long currentLocalIndex; //< Current index (array position) in the current_suboctree.cells[ currentLocalLevel ]

            /**
              * To know what is the left limit on the current level on the current subtree
              * @retrun suboctree.left_limit >> 3 * diff(leaf_index, current_index).
              */
            static long TransposeIndex(const long indexInLeafLevel, const int distanceFromLeafLevel) {
                return indexInLeafLevel >> (3 * distanceFromLeafLevel);
            }


        public:
            /**
            * Constructor
            * @param inTarget the octree to iterate on
            * After building a iterator, this one is positioned at the level 0
            * of the root (level 1 of octree) at the left limit index
            */
            Iterator(FOctree* const inTarget)
                    : currentLocalLevel(0), currentLocalIndex(0) {
                assert(inTarget, "Target for Octree::Iterator cannot be null", __LINE__, __FILE__);
                assert(inTarget->root.getRightLeafIndex() >= 0, "Octree seems to be empty, getRightLeafIndex == 0", __LINE__, __FILE__);

                // Start by the root
                this->current.tree = &inTarget->root;
                // On the left limit
                this->currentLocalIndex = TransposeIndex(this->current.tree->getLeftLeafIndex(), (this->current.tree->getSubOctreeHeight() - this->currentLocalLevel - 1) );
            }

            Iterator() : currentLocalLevel(0), currentLocalIndex(0) {
                current.tree = 0;
            }

            /** Copy constructor
              * @param other source iterator to copy
              */
            Iterator(const Iterator& other){
                this->current = other.current ;
                this->currentLocalLevel = other.currentLocalLevel ;
                this->currentLocalIndex = other.currentLocalIndex ;
            }

            /** Copy operator
              * @param other source iterator to copy
              * @return this after copy
              */
            Iterator& operator=(const Iterator& other){
                this->current = other.current ;
                this->currentLocalLevel = other.currentLocalLevel ;
                this->currentLocalIndex = other.currentLocalIndex ;
                return *this;
            }

            /**
              * Move iterator to the top! (level 0 of root, level 1 of octree)
              * index = left limit at root level
              */
            void gotoTop(){
                while(this->current.tree->hasParent()){
                    this->current.tree = this->current.tree->getParent();
                }
                this->currentLocalLevel = 0;
                this->currentLocalIndex = TransposeIndex(this->current.tree->getLeftLeafIndex(), (this->current.tree->getSubOctreeHeight() - 1) );
            }

            /**
              * Move iterator to the bottom left place
              * We are on a leaf
              */
            void gotoBottomLeft(){
                gotoTop();
                while(1) {
                    this->currentLocalLevel = this->current.tree->getSubOctreeHeight() - 1;
                    this->currentLocalIndex = this->current.tree->getLeftLeafIndex();
                    if( isAtLeafLevel() ){
                        return;
                    }
                    this->current.tree = this->current.middleTree->leafs( this->currentLocalIndex );
                }
            }

            /**
              * Move iterator to the left place at the same level
              */
            void gotoLeft(){
                //  Function variables
                const int currentLevel = level();

                // Goto root sutoctree
                while( this->current.tree->hasParent() ){
                    this->current.tree = this->current.tree->getParent();
                }

                // Go down on the left until arriving on the same global level
                while( this->current.tree->getSubOctreeHeight() + this->current.tree->getSubOctreePosition() - 1 < currentLevel ) {
                    this->current.tree = this->current.middleTree->leafs(this->current.tree->getLeftLeafIndex());
                }

                // Level still unchanged we only go to the left
                // The left limit on this tree at the level we want to stay
                this->currentLocalIndex = TransposeIndex(this->current.tree->getLeftLeafIndex(), (this->current.tree->getSubOctreeHeight() - this->currentLocalLevel - 1));
            }

            /**
              * Move iterator to the left place at the same level
              */
            void gotoRight(){
                //  Function variables
                const int currentLevel = level();
                // Goto root sutoctree
                while( this->current.tree->hasParent() ){
                    this->current.tree = this->current.tree->getParent();
                }
                // Go down on the left until arriving on the same global level
                while( this->current.tree->getSubOctreeHeight() + this->current.tree->getSubOctreePosition() - 1 < currentLevel ) {
                    this->current.tree = this->current.middleTree->leafs(this->current.tree->getRightLeafIndex());
                }
                // Level still unchanged we only go to the left
                // The left limit on this tree at the level we want to stay
                this->currentLocalIndex = TransposeIndex(this->current.tree->getRightLeafIndex(), (this->current.tree->getSubOctreeHeight() - this->currentLocalLevel - 1));
            }

            /**
              * Goto the next value on the right at the same level
              *
              * The algorithm here is :
              * As long as we are on the right limit, go to the parent suboctree
              * if we are on the root and on the right then return (there is no more data on the right)
              *
              * After that point we do not know where we are but we know that there is some data
              * on the right (without knowing our position!)
              *
              * We progress on the brother to find an allocated cell (->)
              * for example if we are on index 2 we will look until 8 = 2 | 7 + 1
              * if we arrive a 8 without finding a cell we go upper and do the same
              * we know we will find something because we are not at the right limit
              *
              * We find an allocated cell.
              * We have to go down, we go on the left child of this cells
              * until : the current level if we did not have change the current suboctree
              * or : the leaf level
              *
              * In the second case, it meanse we need to change octree downward
              * but it is easy because we can use the left limit!
              */
            bool moveRight(){
                //  Function variables
                SubOctreeTypes workingTree = this->current;    // To cover parent other sub octre
                long workingLevel = this->currentLocalLevel;        // To know where we are
                long workingIndex = this->currentLocalIndex;        // To know where we are

                // -- First we have to go in a tree where we can move on the right --
                // Number of time we go into parent subtree
                int countUpward = 0;
                // We stop when we can right move or if there is no more parent (root)
                while( workingIndex == TransposeIndex(workingTree.tree->getRightLeafIndex(), (workingTree.tree->getSubOctreeHeight() - workingLevel - 1) )
                        && workingTree.tree->hasParent() ){
                    // Goto the leaf level into parent at current_tree.position_into_parent_array
                    workingIndex        = workingTree.tree->getIndexInParent();
                    workingTree.tree    = workingTree.tree->getParent();
                    workingLevel        = workingTree.tree->getSubOctreeHeight() - 1;
                    // inc counter
                    ++countUpward;
                }

                // Do we stop because we are on the root (means we cannot move right?)
                if( workingIndex < TransposeIndex(workingTree.tree->getRightLeafIndex(), (workingTree.tree->getSubOctreeHeight() - workingLevel - 1) ) ){
                    // Move to the first right cell pointer(!)
                    ++workingIndex;

                    // Maybe it is null, but we know there is almost one cell on the right
                    // we need to find it
                    if( !workingTree.tree->cellsAt(workingLevel)[workingIndex] ){
                        // While we are not on a allocated cell
                        while( true ){
                            // Test element on the right (test brothers)
                            const long rightLimite = (workingIndex | 7) + 1;
                            while( workingIndex < rightLimite && !workingTree.tree->cellsAt(workingLevel)[workingIndex]){
                                ++workingIndex;
                            }
                            // Stop if we are on an allocated cell
                            if( workingTree.tree->cellsAt(workingLevel)[workingIndex] ){
                                break;
                            }
                            // Else go to the upper level
                            --workingLevel;
                            workingIndex >>= 3;
                        }
                    }

                    // if wokring tree != current tree => working tree leafs level ; else current level
                    const int objectiveLevel = (countUpward ? workingTree.tree->getSubOctreeHeight() - 1 : this->currentLocalLevel );

                    // We need to go down as left as possible
                    while( workingLevel != objectiveLevel ){
                        ++workingLevel;
                        workingIndex <<= 3;
                        const long rightLimite = (workingIndex | 7); // not + 1 because if the 7th first are null it must be the 8th!
                        while( workingIndex < rightLimite && !workingTree.tree->cellsAt(workingLevel)[workingIndex]){
                            ++workingIndex;
                        }
                    }

                    // Do we change from the current sub octree?
                    if( countUpward ){
                        // Then we simply need to go down the same number of time
                        workingTree.tree = workingTree.middleTree->leafs(workingIndex);
                        while( --countUpward ){
                            workingTree.tree = workingTree.middleTree->leafs(workingTree.tree->getLeftLeafIndex());
                        }
                        // level did not change, simpli set octree and get left limit of this octree at the current level
                        this->current = workingTree;
                        this->currentLocalIndex = TransposeIndex(workingTree.tree->getLeftLeafIndex(), (workingTree.tree->getSubOctreeHeight() - this->currentLocalLevel - 1) );
                    }
                    else{
                        // We are on the right tree
                        this->currentLocalIndex = workingIndex;
                    }

                    return true;
                }
                return false;
            }
            /**
              * Move to the upper level
              * It may cause to change the suboctree we are working on
              */
            bool moveUp() {
                // It is on the top level?
                if( this->currentLocalLevel ){
                    // No so simply go up
                    --this->currentLocalLevel;
                    this->currentLocalIndex >>= 3;
                }
                // Yes need to change suboctree
                else if( this->current.tree->hasParent() ){
                    this->currentLocalIndex = this->current.tree->getIndexInParent();
                    this->current.tree = this->current.tree->getParent();
                    this->currentLocalLevel =  this->current.tree->getSubOctreeHeight() - 1;
                }
                else{
                    return false;
                }
                return true;
            }
            /**
              * Move down
              * It may cause to change the suboctree we are working on
              * We point on the first child found from left to right in the above
              * level
              */
            bool moveDown(){
                if( !isAtLeafLevel() ){
                    // We are on the leaf of the current suboctree?
                    if(this->currentLocalLevel + 1 == this->current.tree->getSubOctreeHeight()){
                        // Yes change sub octree
                        this->current.tree = this->current.middleTree->leafs(this->currentLocalIndex);
                        this->currentLocalIndex = 0;
                        this->currentLocalLevel = 0;
                    }
                    // No simply go down
                    else{
                        ++this->currentLocalLevel;
                        this->currentLocalIndex <<= 3;
                    }
                    // Find the first allocated cell from left to right
                    while(!this->current.tree->cellsAt(this->currentLocalLevel)[this->currentLocalIndex]){
                        ++this->currentLocalIndex;
                    }
                    return true;
                }
                return false;
            }


            /**
              * To know if we are not on the root level
              * @return true if we can move up
              */
            bool canProgressToUp() const {
                return this->currentLocalLevel || this->current.tree->hasParent();
            }

            /**
              * To know if we are not on the leafs level
              * @return true if we can move down
              */
            bool canProgressToDown() const {
                return !isAtLeafLevel();
            }

            /**
              * To know if we are on the leafs level
              * @return true if we are at the bottom of the tree
              */
            bool isAtLeafLevel() const {
                return level() + 1 == OctreeHeight;
            }

            /**
              * To know the current level (not local but global)
              * @return the level in the entire octree
              */
            int level() const {
                return this->currentLocalLevel + this->current.tree->getSubOctreePosition();
            }

            /** To access the current particules list
              * You have to be at the leaf level to call this function!
              * @return current element list
              */
            FList<ParticuleClass*>* getCurrentList() {
                return this->current.leafTree->getLeaf(this->currentLocalIndex);
            }

            /** Get the current pointed cell
              * @return current cell element
              */
            CellClass* getCurrentCell() {
                return this->current.tree->cellsAt(this->currentLocalLevel)[this->currentLocalIndex];
            }

            /** Get the child of the current cell
              * This function return an array of CellClass (array size = 8)
              * User has to test each case to know if there is a cell
              * @return the child array
              */
            CellClass** getCurrentChild() {
                // are we at the bottom of the suboctree
                if(this->current.tree->getSubOctreeHeight() - 1 == this->currentLocalLevel ){
                    // then return first level of the suboctree under
                    return &this->current.middleTree->leafs(this->currentLocalIndex)->cellsAt(0)[0];
                }
                else{
                    // else simply return the array at the right position
                    return &this->current.tree->cellsAt(this->currentLocalLevel + 1)[this->currentLocalIndex << 3];
                }
            }

            /** Get the Morton index of the current cell pointed by the iterator
              * @return The global morton index
              * <code>iter.getCurrentGlobalIndex();<br>
              * // is equivalent to :<br>
              * iter.getCurrentCell()->getMortonIndex();</code>
              */
            MortonIndex getCurrentGlobalIndex() const{
                return this->current.tree->cellsAt(this->currentLocalLevel)[this->currentLocalIndex]->getMortonIndex();
            }

        };

        // To be able to access octree root
        friend class Iterator;

        ///////////////////////////////////////////////////////////////////////////
        // This part is related to the FMM algorithm (needed by M2M,M2L,etc.)
        ///////////////////////////////////////////////////////////////////////////

        /** This function return a cell (if it exists) from a morton index and a level
          * @param inIndex the index of the desired cell
          * @param inLevel the level of the desired cell (cannot be infered from the index)
          * @return the cell if it exist or null (0)
          * This function starts from the root until it find a missing cell or the right cell
          */
        CellClass* getCell(const MortonIndex inIndex, const int inLevel) const{
           SubOctreeTypesConst workingTree;
           workingTree.tree = &this->root;
           const MortonIndex treeSubLeafMask = ~(~0x00LL << (3 *  workingTree.tree->getSubOctreeHeight() ));

            // Find the suboctree a the correct level
            while(inLevel >= workingTree.tree->getSubOctreeHeight() + workingTree.tree->getSubOctreePosition()) {
                // compute the leaf index
                const MortonIndex fullIndex = inIndex >> (3 * (inLevel + 1 - (workingTree.tree->getSubOctreeHeight() + workingTree.tree->getSubOctreePosition()) ) );
                // point to next suboctree
                workingTree.tree = workingTree.middleTree->leafs(treeSubLeafMask & fullIndex);
                if(!workingTree.tree) return 0;
            }

            // compute correct index in the array
            const MortonIndex treeLeafMask = ~(~0x00LL << (3 *  (inLevel + 1 - workingTree.tree->getSubOctreePosition()) ));
            return workingTree.tree->cellsAt(inLevel - workingTree.tree->getSubOctreePosition())[treeLeafMask & inIndex];
        }


        /** This function fill a array with the neighbors of a cell
          * it does not put the brothers in the array (brothers are cells
          * at the same level with the same parent) because they are of course
          * direct neighbors.
          * There is a maximum of 26 (3*3*3-1) direct neighbors
          * @param inNeighbors the array to store the elements
          * @param inIndex the index of the element we want the neighbors
          * @param inLevel the level of the element
          * @return the number of neighbors
          */
        int getNeighborsNoBrothers(CellClass* inNeighbors[26], const MortonIndex inIndex, const int inLevel) const {
            FTreeCoordinate center;
            center.setPositionFromMorton(inIndex, inLevel);

            const long limite = FMath::pow(2,inLevel);

            int idxNeighbors = 0;

            // We test all cells around
            for(long idxX = -1 ; idxX <= 1 ; ++idxX){
                if(!FMath::Between(center.getX() + idxX,0l,limite)) continue;

                for(long idxY = -1 ; idxY <= 1 ; ++idxY){
                    if(!FMath::Between(center.getY() + idxY,0l,limite)) continue;

                    for(long idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                        if(!FMath::Between(center.getZ() + idxZ,0l,limite)) continue;

                        // if we are not on the current cell
                        if( !(!idxX && !idxY && !idxZ) ){
                            const FTreeCoordinate other(center.getX() + idxX,center.getY() + idxY,center.getZ() + idxZ);
                            const MortonIndex mortonOther = other.getMortonIndex(inLevel);
                            // if not a brother
                            if( mortonOther>>3 != inIndex>>3 ){
                                // get cell
                                CellClass* const cell = getCell(mortonOther, inLevel);
                                // add to list if not null
                                if(cell) inNeighbors[idxNeighbors++] = cell;
                            }
                        }
                    }
                }
            }

            return idxNeighbors;
        }


        /** This function return an adresse of cell array from a morton index and a level
          *
          * @param inIndex the index of the desired cell array has to contains
          * @param inLevel the level of the desired cell (cannot be infered from the index)
          * @return the cell if it exist or null (0)
          *
          */
        CellClass** getCellPt(const MortonIndex inIndex, const int inLevel) const{
           SubOctreeTypesConst workingTree;
           workingTree.tree = &this->root;

           const MortonIndex treeMiddleMask = ~(~0x00LL << (3 *  workingTree.tree->getSubOctreeHeight() ));

            // Find the suboctree a the correct level
            while(inLevel >= workingTree.tree->getSubOctreeHeight() + workingTree.tree->getSubOctreePosition()) {
                // compute the leaf index
                const MortonIndex fullIndex = inIndex >> 3 * (inLevel + 1 - (workingTree.tree->getSubOctreeHeight() + workingTree.tree->getSubOctreePosition()) );
                // point to next suboctree
                workingTree.tree = workingTree.middleTree->leafs(treeMiddleMask & fullIndex);
                if(!workingTree.tree) return 0;
            }

            // Be sure there is a parent allocated
            const int levelInTree = inLevel - workingTree.tree->getSubOctreePosition();
            if( levelInTree && !workingTree.tree->cellsAt(levelInTree - 1)[~(~0x00LL << (3 * levelInTree )) & (inIndex>>3)]){
                return 0;
            }

            // compute correct index in the array and return the @ in array
            const MortonIndex treeLeafMask = ~(~0x00LL << (3 * (levelInTree + 1 ) ));
            return &workingTree.tree->cellsAt(levelInTree)[treeLeafMask & inIndex];
        }

        /** This function fill an array with the distant neighbors of a cell
          * @param inNeighbors the array to store the elements
          * @param inIndex the index of the element we want the neighbors
          * @param inLevel the level of the element
          * @return the number of neighbors
          */
        int getDistantNeighbors(CellClass* inNeighbors[208], const MortonIndex inIndex, const int inLevel) const{
            // Take the neighbors != brothers
            CellClass* directNeighbors[26];
            const int nbDirectNeighbors = getNeighborsNoBrothers(directNeighbors,inIndex,inLevel);
            // Then take each child of the parent's neighbors if not in directNeighbors
            // Father coordinate
            FTreeCoordinate center;
            center.setPositionFromMorton(inIndex>>3, inLevel-1);
            // Limite at parent level
            const long limite = FMath::pow(2,inLevel-1);

            int idxNeighbors = 0;
            // We test all cells around
            for(long idxX = -1 ; idxX <= 1 ; ++idxX){
                if(!FMath::Between(center.getX() + idxX,0l,limite)) continue;

                for(long idxY = -1 ; idxY <= 1 ; ++idxY){
                    if(!FMath::Between(center.getY() + idxY,0l,limite)) continue;

                    for(long idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                        if(!FMath::Between(center.getZ() + idxZ,0l,limite)) continue;

                        // if we are not on the current cell
                        if( !(!idxX && !idxY && !idxZ) ){

                            const FTreeCoordinate other(center.getX() + idxX,center.getY() + idxY,center.getZ() + idxZ);
                            const MortonIndex mortonOther = other.getMortonIndex(inLevel-1);
                            // Get child
                            CellClass** const cells = getCellPt(mortonOther<<3, inLevel);

                            // If there is one or more child
                            if(cells){
                                // For each child
                                for(int idxCousin = 0 ; idxCousin < 8 ; ++idxCousin){
                                    if(cells[idxCousin]){
                                        // Test if it is a direct neighbor
                                        bool existInDirectNeigh = false;
                                        for(int idxDirectNeigh = 0 ; idxDirectNeigh < nbDirectNeighbors ; ++idxDirectNeigh){
                                            if( cells[idxCousin] == directNeighbors[idxDirectNeigh] ){
                                                existInDirectNeigh = true;
                                                break;
                                            }
                                        }
                                        // add to neighbors
                                        if(!existInDirectNeigh){
                                            inNeighbors[idxNeighbors++] = cells[idxCousin];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            return idxNeighbors;
        }


        /** This function return a cell (if it exists) from a morton index and a level
          * @param inIndex the index of the desired cell
          * @param inLevel the level of the desired cell (cannot be infered from the index)
          * @return the cell if it exist or null (0)
          *
          */
        FList<ParticuleClass*>* getLeaf(const MortonIndex inIndex){
           SubOctreeTypes workingTree;
           workingTree.tree = &this->root;
           const MortonIndex treeSubLeafMask = ~(~0x00LL << (3 *  workingTree.tree->getSubOctreeHeight() ));

            // Find the suboctree a the correct level
            while(leafIndex >= workingTree.tree->getSubOctreeHeight() + workingTree.tree->getSubOctreePosition()) {
                // compute the leaf index
                const MortonIndex fullIndex = inIndex >> (3 * (leafIndex + 1  - (workingTree.tree->getSubOctreeHeight() + workingTree.tree->getSubOctreePosition()) ) );
                // point to next suboctree
                workingTree.tree = workingTree.middleTree->leafs(treeSubLeafMask & fullIndex);
                if(!workingTree.tree) return 0;
            }

            // compute correct index in the array
            const MortonIndex treeLeafMask = ~(~0x00LL << (3 *  (leafIndex + 1 - workingTree.tree->getSubOctreePosition()) ));
            return workingTree.leafTree->getLeaf(treeLeafMask & inIndex);
        }

        /** This function fill an array with the neighbors of a cell
          * @param inNeighbors the array to store the elements
          * @param inIndex the index of the element we want the neighbors
          * @param inLevel the level of the element
          * @return the number of neighbors
          */
        int getLeafsNeighbors(FList<ParticuleClass*>* inNeighbors[26], const MortonIndex inIndex, const int inLevel){
            FTreeCoordinate center;
            center.setPositionFromMorton(inIndex, inLevel);

            const long limite = FMath::pow(2,inLevel);

            int idxNeighbors = 0;

            // We test all cells around
            for(long idxX = -1 ; idxX <= 1 ; ++idxX){
                if(!FMath::Between(center.getX() + idxX,0l,limite)) continue;

                for(long idxY = -1 ; idxY <= 1 ; ++idxY){
                    if(!FMath::Between(center.getY() + idxY,0l,limite)) continue;

                    for(long idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                        if(!FMath::Between(center.getZ() + idxZ,0l,limite)) continue;

                        // if we are not on the current cell
                        if( !(!idxX && !idxY && !idxZ) ){
                            const FTreeCoordinate other(center.getX() + idxX,center.getY() + idxY,center.getZ() + idxZ);
                            const MortonIndex mortonOther = other.getMortonIndex(inLevel);
                            // get cell
                            FList<ParticuleClass*>* const leaf = getLeaf(mortonOther);
                            // add to list if not null
                            if(leaf) inNeighbors[idxNeighbors++] = leaf;
                        }
                    }
                }
            }

            return idxNeighbors;
        }


};

#endif //FOCTREE_HPP
// [--LICENSE--]
