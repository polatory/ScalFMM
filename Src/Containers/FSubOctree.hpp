#ifndef FSUBOCTREE_HPP
#define FSUBOCTREE_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Utils/FGlobal.hpp"
#include "../Utils/F3DPosition.hpp"
#include "../Utils/FAssertable.hpp"
#include "../Utils/FMath.hpp"

#include "FTreeCoordinate.hpp"
#include "FList.hpp"


/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FAbstractSubOctree
 * Please read the license
 *
 * This class is a sub-octree container.
 * This class is the main component of the octree.
 * 
 * This class does not have a main root. In fact, the root
 * is contained in the parent sub-octree.
 *
 * If the sub-octree is a middle suboctree, the leaf level is pointers
 * on other suboctrees.
 *
 * If the sub-octree is a bottom subtree then the leaf level is pointers on
 * lists<particle>
 *
 * This two situations are implemented in two different classes that inherite of FAbstractSubOctree.
 *
 * Please refere to testOctree.cpp to see an example
 * @warning Give the particleClass & cellClass
 */
template< class ParticleClass, class CellClass , template <class ParticleClass> class LeafClass>
class FAbstractSubOctree : protected FAssertable{
protected:

    CellClass*** cells;		            //< Potential cells, cells are allocated only if needed
    FAbstractSubOctree* const parent;       //< Parent suboctree (null for root)

    const long indexInParent;               //< This is the index of the current octree in the parent's array

    long leftLeafIndex;                     //< The leaf at the left position (this is the array index to start when iterate)
    long rightLeafIndex;                    //< The leaf at the right position (this is the last array index when iterate)

    const int subOctreeHeight;              //< Height of this suboctree
    const int subOctreePosition;	    //< Level of the current suboctree in the global tree (0 if node)

    const bool isLeafSubtree;               //< To know if a subtree is leaf or not (we prefere that to a virtual method)


    /**
     * This function compute the morton index for the last level of this suboctree.
     * suppose we have an index like : 000.010.111.001.101.110
     * and now considere that suboctree's height is 2, the morton index have to be cut as :
     * [000.010].[111.001].[101.110] each part correspond to a morton index a leaf level for the
     * different suboctrees.
     * This is why depending on the level of the octree we need to remove extra part on the left and
     * on the right.
     */
    MortonIndex getLeafIndex(const MortonIndex index, const int inTreeHeight) const {
        // Remove right useless part - used by child
        const MortonIndex fullIndex = index >> (3 * (inTreeHeight - (this->subOctreeHeight + this->subOctreePosition) ) );
        // Remove left extra data part - used by parent
        const MortonIndex treeLeafMask = ~(~0x00LL << (3 *  this->subOctreeHeight ));
        return treeLeafMask & fullIndex;
    }

    /**
     * This function creats all the intermediates cells between
     * a leaf and the root.
     * It is used after inserting a new leaf to have cells from leaf to root
     * when computing
     * @param arrayIndex the index at the leaf index of the new element
     */
    void createPreviousCells(MortonIndex arrayIndex, MortonIndex inLeafCellIndex, const FTreeCoordinate& treePosition, const FReal* const inBoxWidthAtLevel){
        int indexLevel = this->subOctreeHeight - 1;
        int bottomToTop = 0;
        while(indexLevel >= 0 && !this->cells[indexLevel][arrayIndex]){
            this->cells[indexLevel][arrayIndex] = new CellClass();
            this->cells[indexLevel][arrayIndex]->setMortonIndex(inLeafCellIndex);

            this->cells[indexLevel][arrayIndex]->setCoordinate(treePosition.getX() >> bottomToTop,
                                                               treePosition.getY() >> bottomToTop,
                                                               treePosition.getZ() >> bottomToTop);

            const int realLevel = indexLevel + this->getSubOctreePosition();
            const FReal widthAtLevel = inBoxWidthAtLevel[realLevel];
            this->cells[indexLevel][arrayIndex]->setPosition(F3DPosition(
                            treePosition.getX() * widthAtLevel + widthAtLevel/2.0,
                            treePosition.getY() * widthAtLevel + widthAtLevel/2.0,
                            treePosition.getZ() * widthAtLevel + widthAtLevel/2.0));

            --indexLevel;
            ++bottomToTop;

            arrayIndex >>= 3;
            inLeafCellIndex >>= 3;
        }
    }

    /**
      * This function is initializing variables when a new leaf is inserted in the tree
      * for example it updates the leaf array marges and calls createPreviousCells()
      * @param arrayIndex the position of the new leaf in the leafs array
      */
    void newLeafInserted(const long arrayIndex, const MortonIndex inLeafCellIndex,  const FTreeCoordinate& host, const FReal* const inBoxWidthAtLevel){
        createPreviousCells(arrayIndex,inLeafCellIndex, host, inBoxWidthAtLevel);
        // Update if this is the bottom left
        if(arrayIndex < this->leftLeafIndex) this->leftLeafIndex = arrayIndex;
        if(arrayIndex > this->rightLeafIndex) this->rightLeafIndex = arrayIndex;
    }

    /** Disable copy */
private:
    FAbstractSubOctree(const FAbstractSubOctree&){}
    FAbstractSubOctree& operator=(const FAbstractSubOctree&){return *this;}

public:
    /**
    * Constructor
    * Allocate the cells arrays to be able to create every potential cells
    * @param inParent the SubOctree parent (0 if node)
    * @param inSubOctreeHeight Height of this suboctree
    * @param inSubOctreePosition Level of the current suboctree in the global tree (1 if upper tree)
    */
    FAbstractSubOctree(FAbstractSubOctree* const inParent, const long inIndexInParent,
                       const int inSubOctreeHeight, const int inSubOctreePosition, const bool inIsLeafSubtree) :
                        cells(0), parent( inParent ), indexInParent(inIndexInParent), leftLeafIndex(1 << (3 * inSubOctreeHeight)), rightLeafIndex(-1),
                        subOctreeHeight( inSubOctreeHeight ), subOctreePosition( inSubOctreePosition ), isLeafSubtree(inIsLeafSubtree) {

        this->cells = new CellClass**[this->subOctreeHeight];
        assert(this->cells, "Allocation failled", __LINE__, __FILE__);

        // We start at a sub-level - 8^1
        long cellsAtlevel = 8;
        for( int indexLevel = 0 ; indexLevel < this->subOctreeHeight ; ++indexLevel ){
            this->cells[indexLevel] = new CellClass*[cellsAtlevel];
            assert(this->cells[indexLevel], "Allocation failled", __LINE__, __FILE__);

            for( int indexCells = 0 ; indexCells < cellsAtlevel ; ++indexCells ){
                this->cells[indexLevel][indexCells] = 0;
            }

            cellsAtlevel <<= 3; // => * 8 >> 8^indexLevel
        }
    }

    /**
    * Destructor
    * Delete cells arrays and allocated cells
    */
    virtual ~FAbstractSubOctree(){
        long cellsAtlevel = 8;
        for( int indexLevel = 0 ; indexLevel < this->subOctreeHeight ; ++indexLevel ){
            for( int indexCells = 0 ; indexCells < cellsAtlevel ; ++indexCells ){
                if(this->cells[indexLevel][indexCells]){
                    delete this->cells[indexLevel][indexCells];
                }
            }

            delete [] this->cells[indexLevel];
            cellsAtlevel <<= 3; // => * 8 >> 8^indexLevel
        }

        delete [] this->cells;
    }


    /**
    * Insert a particle on the subtree
    * @param index the morton index of the particle to insert
    * @param inParticle the particle to insert (must inherite from FAbstractParticle)
    * @param inParticle the inTreeHeight the height of the tree
    */
    virtual void insert(const MortonIndex index, const FTreeCoordinate& host, const ParticleClass& inParticle, const int inTreeHeight, const FReal* const inBoxWidthAtLevel) = 0;

    ///////////////////////////////////////
    // This is the FOctree::Iterator Part
    ///////////////////////////////////////

    /** Suboctree height accessor (leaf level + 1)
      * @return subOctreeHeight */
    int getSubOctreeHeight() const{
        return subOctreeHeight;
    }

    /** Suboctree position in the real tree
      * @return subOctreePosition */
    int getSubOctreePosition() const {
        return subOctreePosition;
    }

    /** Return the more left leaf index
      * the smallest index on the leafs array
      * @return leftLeafIndex */
    long getLeftLeafIndex() const {
        return leftLeafIndex;
    }

    /** Return the more right leaf index
      * the biggest index on the leafs array
      * @return rightLeafIndex */
    long getRightLeafIndex() const {
        return rightLeafIndex;
    }

    /** Return the array of cells at a specious index
      * @param level the level to access cells array (must be < subOctreeHeight)
      * @return cells[level] */
    CellClass** cellsAt(const int level) const{
        assert(level < subOctreeHeight, "Level out of memory", __LINE__, __FILE__);
        return cells[level];
    }

    /** To know if it is the root suboctree
      * @return true if has parent otherwise return false */
    bool hasParent() const {
        return parent;
    }

    /** To get access to the parent suboctree
      * @return parent */
    FAbstractSubOctree* getParent(){
        return parent;
    }

    /** To get access to the parent suboctree (const version)
      * @return parent */
    const FAbstractSubOctree* getParent() const{
        return parent;
    }

    /** To get the index of the current suboctree in the parent leafs array
      * This index can is part of the morton index of the cells contains
      * in this suboctree
      * @return indexInParent */
    long getIndexInParent() const{
        return indexInParent;
    }

    /** To know if this class has been inherited by a leaf subtree or middle subtree
      * Of course we can do a virtual method to do that but virtual method are slower
      * if we point to parent class to call the method
      */
    bool isLeafPart() const{
        return this->isLeafSubtree;
    }
};



/////////////////////////////////////////////////////////////////////////////////////
// Last level sub octree
/////////////////////////////////////////////////////////////////////////////////////

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FSubOctreeWithLeafs
 * Please read the license
 *
 * This class is an sub-octree container.
 * But this is the specialized bottom part, in fact the last level is composed by
 * a cells array (managing by abstract class) and lists of particles.
 *
 * Please refere to testOctree.cpp to see an example.
 * @warning Give the particleClass & cellClass
 */
template< class ParticleClass, class CellClass , template <class ParticleClass> class LeafClass>
class FSubOctreeWithLeafs : public FAbstractSubOctree<ParticleClass,CellClass,LeafClass> {
private:

    LeafClass<ParticleClass>** leafs;            //< Leafs array

    /** Disable copy */
    FSubOctreeWithLeafs(const FSubOctreeWithLeafs&){}
    FSubOctreeWithLeafs& operator=(const FSubOctreeWithLeafs&){return *this;}

public:     
    /**
    * Constructor
    * Allocate the leafs array
    * @param inParent the SubOctree parent (0 if node)
    * @param inSubOctreeHeight Height of this suboctree
    * @param inSubOctreePosition Level of the current suboctree in the global tree (1 if upper tree)
    */
    FSubOctreeWithLeafs(FAbstractSubOctree<ParticleClass,CellClass,LeafClass>* const inParent, const long inIndexInParent,
                        const int inSubOctreeHeight, const int inSubOctreePosition) :
                        FAbstractSubOctree<ParticleClass,CellClass,LeafClass>(inParent, inIndexInParent, inSubOctreeHeight, inSubOctreePosition, true) {

        const long cellsAtLeafLevel = 1 << (3 * inSubOctreeHeight);

        this->leafs = new LeafClass<ParticleClass>*[cellsAtLeafLevel];
        assert(this->leafs, "Allocation failled", __LINE__, __FILE__);

        for( long indexLeaf = 0 ; indexLeaf < cellsAtLeafLevel ; ++indexLeaf ){
            this->leafs[indexLeaf] = 0;
        }
    }

    /**
    * Destructor dealloc all leafs & the leaf array
    */
    virtual ~FSubOctreeWithLeafs(){
        const int cellsAtLeafLevel = 1 << (3 * this->subOctreeHeight );
        for( int indexLeaf = 0 ; indexLeaf < cellsAtLeafLevel ; ++indexLeaf ){
            if(this->leafs[indexLeaf]){
                delete this->leafs[indexLeaf];
            }
        }
        delete [] this->leafs;
    }

    /**
    * Refer to FAbstractSubOctree::insert
    */
    void insert(const MortonIndex index, const FTreeCoordinate& host, const ParticleClass& inParticle, const int inTreeHeight, const FReal* const inBoxWidthAtLevel){
        // Get the morton index for the leaf level
        const MortonIndex arrayIndex = FAbstractSubOctree<ParticleClass,CellClass,LeafClass>::getLeafIndex(index,inTreeHeight);
        // is there already a leaf?
        if( !this->leafs[arrayIndex] ){
            this->leafs[arrayIndex] = new LeafClass<ParticleClass>();

            FAbstractSubOctree<ParticleClass,CellClass,LeafClass>::newLeafInserted( arrayIndex , index, host, inBoxWidthAtLevel);
        }
        // add particle to leaf list
        this->leafs[arrayIndex]->push(inParticle);
    }

    /** To get access to leafs elements
      * @param index the position of the leaf
      * @return the list of particles at this index */
    FList<ParticleClass>* getLeafSrc(const int index){
        LeafClass<ParticleClass>* const leaf = this->leafs[index];
        return (leaf ? leaf->getSrc(): 0);
    }

    /** To get access to leafs elements
      * @param index the position of the leaf
      * @return the list of particles at this index */
    FList<ParticleClass>* getLeafTargets(const int index){
        LeafClass<ParticleClass>* const leaf = this->leafs[index];
        return (leaf ? leaf->getTargets(): 0);
    }

    /** To get access to leafs elements
      * @param index the position of the leaf
      * @return the list of particles at this index */
    const FList<ParticleClass>* getLeafSrc(const int index) const {
        LeafClass<ParticleClass>* const leaf = this->leafs[index];
        return (leaf ? leaf->getSrc(): 0);
    }

    /** To get access to leafs elements
      * @param index the position of the leaf
      * @return the list of particles at this index */
    const FList<ParticleClass>* getLeafTargets(const int index) const {
        LeafClass<ParticleClass>* const leaf = this->leafs[index];
        return (leaf ? leaf->getTargets() : 0);
    }

};



/////////////////////////////////////////////////////////////////////////////////////
// Middle level sub octree
/////////////////////////////////////////////////////////////////////////////////////

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FSubOctree
 * Please read the license
 *
 * This class is an sub-octree container.
 * This is the middle level specialized suboctree, it means that it does not contain
 * leaf but pointers to other suboctree.
 * These suboctrees at the last level can be FSubOctree of FSubOctreeWithLeafs depending
 * if they are at the bottom of the tree or not
 *
 * @warning Give the particleClass & cellClass
 */
template< class ParticleClass, class CellClass , template <class ParticleClass> class LeafClass>
class FSubOctree : public FAbstractSubOctree<ParticleClass,CellClass,LeafClass> {
private:
    FAbstractSubOctree<ParticleClass,CellClass,LeafClass>** subleafs;    //< Last levels is composed of suboctree

    /** Disable copy */
    FSubOctree(const FSubOctree&){}
    FSubOctree& operator=(const FSubOctree&){return *this;}

public:	
    /**
    * Constructor
    * Allocate the subleafs array
    * @param inParent the SubOctree parent (0 if node)
    * @param inSubOctreeHeight Height of this suboctree
    * @param inSubOctreePosition Level of the current suboctree in the global tree (0 if upper tree)
    */
    FSubOctree(FAbstractSubOctree<ParticleClass,CellClass,LeafClass>* const inParent,  const long inIndexInParent,
               const int inSubOctreeHeight, const int inSubOctreePosition) :
            FAbstractSubOctree<ParticleClass,CellClass,LeafClass>(inParent, inIndexInParent, inSubOctreeHeight, inSubOctreePosition, false) {

        const long cellsAtLeafLevel = 1 << (3 * inSubOctreeHeight);
        this->subleafs = new FAbstractSubOctree<ParticleClass,CellClass,LeafClass>*[cellsAtLeafLevel];
        assert(this->subleafs, "Allocation failled", __LINE__, __FILE__);

        for( int indexLeaf = 0 ; indexLeaf < cellsAtLeafLevel ; ++indexLeaf ){
            this->subleafs[indexLeaf] = 0;
        }
    }

    /**
    * Destructor dealloc all suboctrees leafs & leafs array
    */
    virtual ~FSubOctree(){
        const long cellsAtLeafLevel = 1 << (3 * this->subOctreeHeight);
        for( int indexLeaf = 0 ; indexLeaf < cellsAtLeafLevel ; ++indexLeaf ){
            if(this->subleafs[indexLeaf]) delete this->subleafs[indexLeaf];
        }
        delete [] this->subleafs;
    }

    /**
    * Refer to FAbstractSubOctree::insert
    */
    void insert(const MortonIndex index, const FTreeCoordinate& host, const ParticleClass& inParticle, const int inTreeHeight, const FReal* const inBoxWidthAtLevel){
        // We need the morton index at the bottom level of this sub octree
        // so we remove the right side
        const MortonIndex arrayIndex = FAbstractSubOctree<ParticleClass,CellClass,LeafClass>::getLeafIndex(index,inTreeHeight);
        // Is there already a leaf?
        if( !this->subleafs[arrayIndex] ){
            // We need to create leaf sub octree
            const int nextSubOctreePosition = this->subOctreePosition + this->subOctreeHeight;
            const int nextSubOctreeHeight = FMath::Min(inTreeHeight - nextSubOctreePosition, this->subOctreeHeight);

            // Next suboctree is a middle suboctree
            if(inTreeHeight > nextSubOctreeHeight + nextSubOctreePosition){
                this->subleafs[arrayIndex] = new FSubOctree(this,arrayIndex,nextSubOctreeHeight,nextSubOctreePosition);
            }
            // Or next suboctree contains the reail leaf!
            else{
                this->subleafs[arrayIndex] = new FSubOctreeWithLeafs<ParticleClass,CellClass,LeafClass>(this,arrayIndex,nextSubOctreeHeight,nextSubOctreePosition);
            }

            const FTreeCoordinate hostAtLevel(
                        host.getX() >> (inTreeHeight - nextSubOctreePosition ),
                        host.getY() >> (inTreeHeight - nextSubOctreePosition ),
                        host.getZ() >> (inTreeHeight - nextSubOctreePosition ));

            // We need to inform parent class
            FAbstractSubOctree<ParticleClass,CellClass,LeafClass>::newLeafInserted( arrayIndex, index >> (3 * (inTreeHeight-nextSubOctreePosition) ), hostAtLevel, inBoxWidthAtLevel);
        }
        // Ask next suboctree to insert the particle
        this->subleafs[arrayIndex]->insert( index, host, inParticle, inTreeHeight, inBoxWidthAtLevel);
    }

    /** To get access to leafs elements (child suboctree)
      * @param index the position of the leaf/child suboctree
      * @return child at this index */
    FAbstractSubOctree<ParticleClass,CellClass,LeafClass>* leafs(const int index) {
        return this->subleafs[index];
    }

    /** To get access to leafs elements (child suboctree)
      * @param index the position of the leaf/child suboctree
      * @return child at this index */
    const FAbstractSubOctree<ParticleClass,CellClass,LeafClass>* leafs(const int index) const {
        return this->subleafs[index];
    }
};


#endif //FSUBOCTREE_HPP
// [--LICENSE--]
