#ifndef FOCTREE_HPP
#define FOCTREE_HPP
// [--License--]

#include "FSubOctree.hpp"

#include "../Utils/FDebug.hpp"
#include "../Utils/FGlobal.hpp"
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
 * FOctree<TestParticle, TestCell> tree(1.0,F3DPosition(0.5,0.5,0.5));
 * </code>
 *
 * Particles and cells has to respect the Abstract class definition.
 * Particle must extend {FExtendPosition}
 * Cell must extend extend {FExtendPosition,FExtendMortonIndex}
 */
template< class ParticleClass, class CellClass, class ContainerClass, class LeafClass>
class FOctree : protected FAssertable {
    FReal*const boxWidthAtLevel;		//< to store the width of each boxs at all levels

    const int height;		//< tree height
    const int subHeight;		//< tree height
    const int leafIndex;		//< index of leaf int array

    FSubOctree< ParticleClass, CellClass , ContainerClass, LeafClass> root;   //< root suboctree

    const F3DPosition boxCenter;	//< the space system center
    const F3DPosition boxCorner;	//< the space system corner (used to compute morton index)

    const FReal boxWidth;          //< the space system width


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
    FTreeCoordinate getCoordinateFromPosition(const F3DPosition& inPosition) const {
        // box coordinate to host the particle
        FTreeCoordinate host;
        // position has to be relative to corner not center
        host.setX( getTreeCoordinate( inPosition.getX() - this->boxCorner.getX() ));
        host.setY( getTreeCoordinate( inPosition.getY() - this->boxCorner.getY() ));
        host.setZ( getTreeCoordinate( inPosition.getZ() - this->boxCorner.getZ() ));
        return host;
    }

    /**
        * Get the box number from a position
        * at a position POS with a leaf level box width of WI, the position is RELATIVE_TO_CORNER(POS)/WI
        * @param inRelativePosition a position from the corner of the box
        * @return the box num at the leaf level that contains inRelativePosition
        */
    int getTreeCoordinate(const FReal inRelativePosition) const {
        FDEBUG( fassert(inRelativePosition >= 0 && inRelativePosition < this->boxWidth, "Particle out of box", __LINE__, __FILE__) );
        const FReal indexFReal = inRelativePosition / this->boxWidthAtLevel[this->leafIndex];
        const int index = int(FMath::dfloor(indexFReal));
        if( index && FMath::LookEqual(inRelativePosition, this->boxWidthAtLevel[this->leafIndex] * FReal(index) ) ){
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
    FOctree(const int inHeight, const int inSubHeight,
            const FReal inBoxWidth, const F3DPosition& inBoxCenter)
        : boxWidthAtLevel(new FReal[inHeight]),
          height(inHeight) , subHeight(inSubHeight), leafIndex(this->height-1),
          root(0, 0, this->subHeight, 1), boxCenter(inBoxCenter), boxCorner(inBoxCenter - (inBoxWidth/2)), boxWidth(inBoxWidth)
    {
        FReal tempWidth = this->boxWidth;
        // pre compute box width for each level
        for(int indexLevel = 0; indexLevel < this->height; ++indexLevel ){
            this->boxWidthAtLevel[indexLevel] = tempWidth;
            tempWidth /= FReal(2.0);
        }
    }

    /** Use a loader to be filled
      * @param the loader to fill the current tree
      */
    template <class LoaderClass>
    void fillWithLoader(LoaderClass& loader){
        ParticleClass particleToFill;
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(particleToFill);
            insert(particleToFill);
        }
    }

    /** Desctructor */
    virtual ~FOctree() {
        delete [] boxWidthAtLevel;
    }

    /** To get the tree height */
    int getHeight() const {
        return this->height;
    }

    /** To get the tree subheight */
    int getSubHeight() const{
        return this->subHeight;
    }

    /** To get the box width */
    FReal getBoxWidth() const{
        return this->boxWidth;
    }

    /** To get the center of the box */
    const F3DPosition& getBoxCenter() const{
        return this->boxCenter;
    }

    /**
 * Insert a particle on the tree
 * algorithm is :
 * Compute morton index for the particle
 * ask node to insert this particle
 * @param inParticle the particle to insert (must inherite from FAbstractParticle)
 */
    void insert(const ParticleClass& inParticle){
        const FTreeCoordinate host = getCoordinateFromPosition( inParticle.getPosition() );
        const MortonIndex particleIndex = host.getMortonIndex(leafIndex);
        root.insert( particleIndex, host, inParticle, this->height, this->boxWidthAtLevel);
    }

    /** Remove a leaf from its morton index
          * @param the index of the leaf to remove
          */
    void removeLeaf(const MortonIndex indexToRemove ){
        root.removeLeaf( indexToRemove , this->height);
    }

    /**
        * Get a morton index from a real position
        * @param a position to compute MI
        * @return the morton index
        */
    MortonIndex getMortonFromPosition(const F3DPosition& position) const {
        return getCoordinateFromPosition(position).getMortonIndex(leafIndex);
    }

    /**
          * The class works on suboctree. Most of the resources needed
          * are avaiblable by using FAbstractSubOctree. But when accessing
          * to the leaf we have to use FSubOctree or FSubOctreeWithLeafs
          * depending if we are working on the bottom of the tree.
          */
    union SubOctreeTypes {
        FAbstractSubOctree<ParticleClass,CellClass,ContainerClass,LeafClass>* tree;     //< Usual pointer to work
        FSubOctree<ParticleClass,CellClass,ContainerClass,LeafClass>* middleTree;       //< To access to sub-octree under
        FSubOctreeWithLeafs<ParticleClass,CellClass,ContainerClass,LeafClass>* leafTree;//< To access to particles lists
    };

    /**
          * This class is a const SubOctreeTypes
          */
    union SubOctreeTypesConst {
        const FAbstractSubOctree<ParticleClass,CellClass,ContainerClass,LeafClass>* tree;     //< Usual pointer to work
        const FSubOctree<ParticleClass,CellClass,ContainerClass,LeafClass>* middleTree;       //< To access to sub-octree under
        const FSubOctreeWithLeafs<ParticleClass,CellClass,ContainerClass,LeafClass>* leafTree;//< To access to particles lists
    };

    /**
          * This has to be used to iterate on an octree
          * It simply stores an pointer on a suboctree and moves to right/left/up/down.
          * Please refere to testOctreeIter file to see how it works.
          *
          * <code>
          * FOctree<TestParticle, TestCell, NbLevels, NbSubLevels>::Iterator octreeIterator(&tree); <br>
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

        int currentLocalIndex; //< Current index (array position) in the current_suboctree.cells[ currentLocalLevel ]
        int currentLocalLevel;  //< Current level in the current suboctree

        /**
              * To know what is the left limit on the current level on the current subtree
              * @retrun suboctree.left_limit >> 3 * diff(leaf_index, current_index).
              */
        static int TransposeIndex(const int indexInLeafLevel, const int distanceFromLeafLevel) {
            return indexInLeafLevel >> (3 * distanceFromLeafLevel);
        }


    public:
        /**
            * Constructor
            * @param inTarget the octree to iterate on
            * After building a iterator, this one is positioned at the level 0
            * of the root (level 1 of octree) at the left limit index
            */
        explicit Iterator(FOctree* const inTarget)
            : currentLocalIndex(0) , currentLocalLevel(0) {
            fassert(inTarget, "Target for Octree::Iterator cannot be null", __LINE__, __FILE__);
            fassert(inTarget->root.getRightLeafIndex() >= 0, "Octree seems to be empty, getRightLeafIndex == 0", __LINE__, __FILE__);

            // Start by the root
            this->current.tree = &inTarget->root;
            // On the left limit
            this->currentLocalIndex = TransposeIndex(this->current.tree->getLeftLeafIndex(), (this->current.tree->getSubOctreeHeight() - this->currentLocalLevel - 1) );
        }

        Iterator() : currentLocalIndex(0),currentLocalLevel(0) {
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
              * Move iterator to the top! (level 0 of root suboctree, level 1 of octree)
              * after this function : index = left limit at root level
              * the Algorithm is :
              *     going to root suboctree
              *     going to the first level and most left node
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
              * We are on a leaf a the most left node
              * the Algorithm is :
              *     first go to top
              *     then stay on the left and go downward
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
              * if needed we go on another suboctree but we stay on at the same level
              * the Algorithm is :
              *     go to top
              *     go downard until we are a the same level
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
              * Move iterator to the right place at the same level
              * if needed we go on another suboctree but we stay on at the same level
              * the Algorithm is :
              *     go to top
              *     go downard until we are a the same level
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
              * We gotoNext on the brother to find an allocated cell (->)
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
              *
              * @return true if we succeed to go to the right, else false
              */
        bool moveRight(){
            //  Function variables
            SubOctreeTypes workingTree = this->current;    // To cover parent other sub octre
            int workingLevel = this->currentLocalLevel;        // To know where we are
            int workingIndex = this->currentLocalIndex;        // To know where we are

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
                        const int rightLimite = (workingIndex | 7) + 1;
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
                    const int rightLimite = (workingIndex | 7); // not + 1 because if the 7th first are null it must be the 8th!
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
              * but we are using the same morton index >> 3
              * @return true if succeed
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
              * @return true if succeed
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
            return this->current.tree->isLeafPart() && this->currentLocalLevel + 1 == this->current.tree->getSubOctreeHeight();
        }

        /**
              * To know the current level (not local but global)
              * @return the level in the entire octree
              */
        int level() const {
            return this->currentLocalLevel + this->current.tree->getSubOctreePosition();
        }

        /** Get the current pointed leaf
              * @return current leaf element
              */
        LeafClass* getCurrentLeaf() const {
            return this->current.leafTree->getLeaf(this->currentLocalIndex);
        }

        /** To access the current particles list
              * You have to be at the leaf level to call this function!
              * @return current element list
              */
        ContainerClass* getCurrentListSrc() const {
            return this->current.leafTree->getLeafSrc(this->currentLocalIndex);
        }

        /** To access the current particles list
              * You have to be at the leaf level to call this function!
              * @return current element list
              */
        ContainerClass* getCurrentListTargets() const {
            return this->current.leafTree->getLeafTargets(this->currentLocalIndex);
        }

        /** Get the current pointed cell
              * @return current cell element
              */
        CellClass* getCurrentCell() const {
            return this->current.tree->cellsAt(this->currentLocalLevel)[this->currentLocalIndex];
        }

        /** Get the child of the current cell
              * This function return an array of CellClass (array size = 8)
              * User has to test each case to know if there is a cell
              * @return the child array
              */
        CellClass** getCurrentChild() const {
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

        /** Get the part of array that contains all the pointers
              *
              */
        CellClass** getCurrentBox() const {
            return &this->current.tree->cellsAt(this->currentLocalLevel)[this->currentLocalIndex & ~7];
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

        /** To get the tree coordinate of the current working cell
              *
              */
        const FTreeCoordinate& getCurrentGlobalCoordinate() const{
            return this->current.tree->cellsAt(this->currentLocalLevel)[this->currentLocalIndex]->getCoordinate();
        }

    };

    // To be able to access octree root & data
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
          *  // Take the neighbors != brothers
          *  CellClass* directNeighbors[26];
          *  const int nbDirectNeighbors = getNeighborsNoBrothers(directNeighbors,inIndex,inLevel);
          * @param inNeighbors the array to store the elements
          * @param inIndex the index of the element we want the neighbors
          * @param inLevel the level of the element
          * @return the number of neighbors
          */
    int getNeighborsNoBrothers(CellClass* inNeighbors[26], const MortonIndex inIndex, const int inLevel) const {
        FTreeCoordinate center;
        center.setPositionFromMorton(inIndex, inLevel);

        const int limite = FMath::pow(2,inLevel);

        int idxNeighbors = 0;

        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(center.getX() + idxX,0,limite)) continue;

            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(center.getY() + idxY,0,limite)) continue;

                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(center.getZ() + idxZ,0,limite)) continue;

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
            workingTree.tree = workingTree.middleTree->leafs(int(treeMiddleMask & fullIndex));
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
          * @param inNeighborsIndex the array to store morton index of the neighbors
          * @param inIndex the index of the element we want the neighbors
          * @param inLevel the level of the element
          * @return the number of neighbors
          */
    int getDistantNeighbors(const CellClass* inNeighbors[189],
                            const FTreeCoordinate& workingCell,
                            const int inLevel) const{

        // Then take each child of the parent's neighbors if not in directNeighbors
        // Father coordinate
        const FTreeCoordinate parentCell(workingCell.getX()>>1,workingCell.getY()>>1,workingCell.getZ()>>1);

        // Limite at parent level number of box (split by 2 by level)
        const int limite = FMath::pow(2,inLevel-1);

        int idxNeighbors = 0;
        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(parentCell.getX() + idxX,0,limite)) continue;

            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(parentCell.getY() + idxY,0,limite)) continue;

                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(parentCell.getZ() + idxZ,0,limite)) continue;

                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){

                        const FTreeCoordinate otherParent(parentCell.getX() + idxX,parentCell.getY() + idxY,parentCell.getZ() + idxZ);
                        const MortonIndex mortonOtherParent = otherParent.getMortonIndex(inLevel-1) << 3;
                        // Get child
                        CellClass** const cells = getCellPt(mortonOtherParent, inLevel);

                        // If there is one or more child
                        if(cells){
                            // For each child
                            for(int idxCousin = 0 ; idxCousin < 8 ; ++idxCousin){
                                if(cells[idxCousin]){
                                    //FTreeCoordinate potentialNeighbor;
                                    //potentialNeighbor.setPositionFromMorton(mortonOtherParent | idxCousin, inLevel);
                                    const FTreeCoordinate potentialNeighbor((otherParent.getX()<<1) | (idxCousin>>2 & 1),
                                                                            (otherParent.getY()<<1) | (idxCousin>>1 & 1),
                                                                            (otherParent.getZ()<<1) | (idxCousin&1));

                                    // Test if it is a direct neighbor
                                    if(FMath::Abs(workingCell.getX() - potentialNeighbor.getX()) > 1 ||
                                            FMath::Abs(workingCell.getY() - potentialNeighbor.getY()) > 1 ||
                                            FMath::Abs(workingCell.getZ() - potentialNeighbor.getZ()) > 1){
                                        // add to neighbors
                                        inNeighbors[idxNeighbors] = cells[idxCousin];
                                        ++idxNeighbors;
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


    /** This function fill an array with the distant neighbors of a cell
          * it respects the periodic condition and will give the relative distance
          * between the working cell and the neighbors
          * @param inNeighbors the array to store the elements
          * @param inRelativePosition the array to store the relative position of the neighbors
          * @param workingCell the index of the element we want the neighbors
          * @param inLevel the level of the element
          * @return the number of neighbors
          */
    int getDistantNeighbors(const CellClass* inNeighbors[189], FTreeCoordinate inRelativePosition[189],
                            const FTreeCoordinate& workingCell,
                            const int inLevel) const{

        // Then take each child of the parent's neighbors if not in directNeighbors
        // Father coordinate
        const FTreeCoordinate parentCell(workingCell.getX()>>1,workingCell.getY()>>1,workingCell.getZ()>>1);

        // Limite at parent level number of box (split by 2 by level)
        const int limite = FMath::pow(2,inLevel-1);


        int idxNeighbors = 0;
        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){

                        const FTreeCoordinate otherParent(parentCell.getX() + idxX,parentCell.getY() + idxY,parentCell.getZ() + idxZ);
                        FTreeCoordinate otherParentInBox(otherParent);
                        // periodic
                        if( otherParentInBox.getX() < 0 ) otherParentInBox.setX( otherParentInBox.getX() + limite );
                        else if( limite <= otherParentInBox.getX() ) otherParentInBox.setX( otherParentInBox.getX() - limite );

                        if( otherParentInBox.getY() < 0 ) otherParentInBox.setY( otherParentInBox.getY() + limite );
                        else if( limite <= otherParentInBox.getY() ) otherParentInBox.setY( otherParentInBox.getY() - limite );

                        if( otherParentInBox.getZ() < 0 ) otherParentInBox.setZ( otherParentInBox.getZ() + limite );
                        else if( limite <= otherParentInBox.getZ() ) otherParentInBox.setZ( otherParentInBox.getZ() - limite );


                        const MortonIndex mortonOtherParent = otherParentInBox.getMortonIndex(inLevel-1) << 3;
                        // Get child
                        CellClass** const cells = getCellPt(mortonOtherParent, inLevel);

                        // If there is one or more child
                        if(cells){
                            // For each child
                            for(int idxCousin = 0 ; idxCousin < 8 ; ++idxCousin){
                                if(cells[idxCousin]){
                                    //FTreeCoordinate potentialNeighbor;
                                    //potentialNeighbor.setPositionFromMorton(mortonOtherParent | idxCousin, inLevel);
                                    const FTreeCoordinate potentialNeighbor((otherParent.getX()<<1) | (idxCousin>>2 & 1),
                                                                            (otherParent.getY()<<1) | (idxCousin>>1 & 1),
                                                                            (otherParent.getZ()<<1) | (idxCousin&1));

                                    const FTreeCoordinate relativePosition(potentialNeighbor.getX() - workingCell.getX(),
                                                                           potentialNeighbor.getY() - workingCell.getY(),
                                                                           potentialNeighbor.getZ() - workingCell.getZ());

                                    // Test if it is a direct neighbor
                                    if(FMath::Abs(relativePosition.getX()) > 1 ||
                                            FMath::Abs(relativePosition.getY()) > 1 ||
                                            FMath::Abs(relativePosition.getZ()) > 1){
                                        // add to neighbors
                                        inNeighbors[idxNeighbors] = cells[idxCousin];
                                        inRelativePosition[idxNeighbors] = relativePosition;
                                        ++idxNeighbors;
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
    ContainerClass* getLeafSrc(const MortonIndex inIndex){
        SubOctreeTypes workingTree;
        workingTree.tree = &this->root;
        const MortonIndex treeSubLeafMask = ~(~0x00LL << (3 *  workingTree.tree->getSubOctreeHeight() ));

        // Find the suboctree a the correct level
        while(leafIndex >= workingTree.tree->getSubOctreeHeight() + workingTree.tree->getSubOctreePosition()) {
            // compute the leaf index
            const MortonIndex fullIndex = inIndex >> (3 * (leafIndex + 1  - (workingTree.tree->getSubOctreeHeight() + workingTree.tree->getSubOctreePosition()) ) );
            // point to next suboctree
            workingTree.tree = workingTree.middleTree->leafs(int(treeSubLeafMask & fullIndex));
            if(!workingTree.tree) return 0;
        }

        // compute correct index in the array
        const MortonIndex treeLeafMask = ~(~0x00LL << (3 *  (leafIndex + 1 - workingTree.tree->getSubOctreePosition()) ));
        return workingTree.leafTree->getLeafSrc(int(treeLeafMask & inIndex));
    }

    /** This function fill an array with the neighbors of a cell
          * @param inNeighbors the array to store the elements
          * @param inIndex the index of the element we want the neighbors
          * @param inLevel the level of the element
          * @return the number of neighbors
          */
    int getLeafsNeighbors(ContainerClass* inNeighbors[26], const MortonIndex inIndex, const int inLevel){
        MortonIndex inNeighborsIndex[26];
        return getLeafsNeighborsWithIndex(inNeighbors, inNeighborsIndex, inIndex, inLevel);
    }

    /** This function fill an array with the neighbors of a cell
          * @param inNeighbors the array to store the elements
          * @param inIndex the index of the element we want the neighbors
          * @param inLevel the level of the element
          * @return the number of neighbors
          */
    int getLeafsNeighborsWithIndex(ContainerClass* inNeighbors[26], MortonIndex inNeighborsIndex[26], const MortonIndex inIndex, const int inLevel){
        FTreeCoordinate center;
        center.setPositionFromMorton(inIndex, inLevel);

        const int limite = FMath::pow(2,inLevel);

        int idxNeighbors = 0;

        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){
            if(!FMath::Between(center.getX() + idxX,0,limite)) continue;

            for(int idxY = -1 ; idxY <= 1 ; ++idxY){
                if(!FMath::Between(center.getY() + idxY,0,limite)) continue;

                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){
                    if(!FMath::Between(center.getZ() + idxZ,0,limite)) continue;

                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){
                        const FTreeCoordinate other(center.getX() + idxX,center.getY() + idxY,center.getZ() + idxZ);
                        const MortonIndex mortonOther = other.getMortonIndex(inLevel);
                        // get cell
                        ContainerClass* const leaf = getLeafSrc(mortonOther);
                        // add to list if not null
                        if(leaf){
                            inNeighborsIndex[idxNeighbors] = mortonOther;
                            inNeighbors[idxNeighbors++] = leaf;
                        }
                    }
                }
            }
        }

        return idxNeighbors;
    }

    /** This function fill an array with the neighbors of a cell
          * @param inNeighbors the array to store the elements
          * @param inIndex the index of the element we want the neighbors
          * @param inLevel the level of the element
          * @return the number of neighbors
          */
    int getLeafsNeighborsWithIndex(ContainerClass* inNeighbors[26], FTreeCoordinate inNeighborsPosition[26], const MortonIndex inIndex, const int inLevel){
        FTreeCoordinate center;
        center.setPositionFromMorton(inIndex, inLevel);

        const int limite = FMath::pow(2,inLevel);

        int idxNeighbors = 0;

        // We test all cells around
        for(int idxX = -1 ; idxX <= 1 ; ++idxX){

            for(int idxY = -1 ; idxY <= 1 ; ++idxY){

                for(int idxZ = -1 ; idxZ <= 1 ; ++idxZ){

                    // if we are not on the current cell
                    if( idxX || idxY || idxZ ){
                        FTreeCoordinate other(center.getX() + idxX,center.getY() + idxY,center.getZ() + idxZ);

                        // To give the orientation of the neighbors
                        FTreeCoordinate offset;

                        if( other.getX() < 0 ){
                            other.setX( other.getX() + limite );
                            offset.setX(-1);
                        }
                        else if( limite <= other.getX() ){
                            other.setX( other.getX() - limite );
                            offset.setX(1);
                        }
                        if( other.getY() < 0 ){
                            other.setY( other.getY() + limite );
                            offset.setY(-1);
                        }
                        else if( limite <= other.getY() ){
                            other.setY( other.getY() - limite );
                            offset.setY(1);
                        }
                        if( other.getZ() < 0 ){
                            other.setZ( other.getZ() + limite );
                            offset.setZ(-1);
                        }
                        else if( limite <= other.getZ() ){
                            other.setZ( other.getZ() - limite );
                            offset.setZ(1);
                        }

                        const MortonIndex mortonOther = other.getMortonIndex(inLevel);
                        // get cell
                        ContainerClass* const leaf = getLeafSrc(mortonOther);
                        // add to list if not null
                        if(leaf){
                            inNeighborsPosition[idxNeighbors] = offset;
                            inNeighbors[idxNeighbors++] = leaf;
                        }
                    }
                }
            }
        }

        return idxNeighbors;
    }

};

#endif //FOCTREE_HPP
// [--END--]
