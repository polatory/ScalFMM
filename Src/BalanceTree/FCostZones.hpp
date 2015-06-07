// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _COSTZONES_HPP_
#define _COSTZONES_HPP_

#include "FCostCell.hpp"
#include "FCoordColour.hpp"

#include <vector>
#include <stdexcept>
#include <sstream>

/**
 * \brief The costzones algorithm implementation.
 * \author Quentin Khan <quentin.khan@inria.fr>
 *
 * This class is an implementation of the costzones algorithm described in "A
 * Parallel Adaptive Fast Multipole Method" (1993). The algorithm consists in an
 * in-order traversal of the octree where cell costs are accumulated. When an
 * accumulation is too big, a new zone is created.
 *
 * It is possible to set the levels of the tree that must be considered to
 * compute the costzones.
 *
 * \tparam OctreeClass The type of the octree to work on.
 * \tparam CellClass   The type of the cells to work with.
 * \parblock
 * This class must provide a typedef named CostType that identifies
 * the type used to store data.
 * \endparblock
 */
template<typename OctreeClass, typename CellClass>
class FCostZones {
public:
    using CostType = typename CellClass::costtype;
    using TreeIterator = typename OctreeClass::Iterator;

    /**
     * \brief Class used to store the bounds of a zone.
     * The bounds consist in the Morton index of the first node and the number
     * of subsequent nodes.
     */ 
    using BoundClass = std::pair<TreeIterator, int>;
    /// Initial value for empty bounds.
    const BoundClass _boundInit {TreeIterator(),0};


protected:
    /// Alias to FCoordColour#coord2colour -*- not optimised
    constexpr static int(*coord2colour)(const FTreeCoordinate&)
    = FCoordColour::coord2colour;

    /**
     * \brief Enumeration to specify the children to move to during the in-order
     * traversal.
     */
    enum ChildrenSide {LEFT, RIGHT};


    /// The iterator to move through the tree.
    typename OctreeClass::Iterator _it;
    /// The number of zones to create.
    int _nbZones;
    /// The tree height.
    int _treeHeight;
    /// The highest level in the tree in which to consider cells (inclusive)
    int _topMostLevel = 0;
    /// The lowest level in the tree in which to consider cells (exclusive)
    int _bottomMostLevel = 1;


    /// The current cumulative cost of visited cells.
    CostType _internalCurrentCost = 0;
    /// The total cost of internal cell.
    CostType _internalTotalCost = 0;
    /**
     * \brief The vector containing the boundaries of the internal node zones
     *
     * Sorted by zone > level > bounds.
     * \details This means there are _treeHeight elements in the inner vectors.
     */
    std::vector< std::vector< BoundClass > > _internalZoneBounds;


    /// The current cumulative cost of the visited leaves.
    std::vector<CostType> _leafCurrentCost;
    /// The total cost of the leaves.
    std::vector<CostType> _leafTotalCost;
    /**
     * \brief The vector containing the boundaries of the leaf zones
     *
     * Sorted by zone > colour > bounds.
     * \details This means there are FCoordColour::range elements in the inner
     * vectors.
     */
    std::vector< std::vector< BoundClass > > _leafZoneBounds;

    /// The vector containing the costzone cells sorted by zone > level > cells.
    std::vector< std::vector< std::pair<int, CellClass*> > > _zones;



public:

    /**
     * \brief Constructor
     * \param tree    The tree to work on.
     * \param nbZones The number of zones to create.
     */
    FCostZones(OctreeClass* tree, int nbZones) : 
        _it( tree ),
        _nbZones( nbZones ),
        _treeHeight( tree->getHeight() ),
        _bottomMostLevel(_treeHeight),
        _internalZoneBounds(
            _nbZones,
            std::vector< BoundClass >(_treeHeight, _boundInit )),
        _leafCurrentCost( FCoordColour::range, 0),
        _leafTotalCost( FCoordColour::range, 0),
        _leafZoneBounds(
            _nbZones,
            std::vector< BoundClass >(FCoordColour::range, _boundInit)),
        _zones(
            _nbZones,
            std::vector< std::pair< int, CellClass*> >( ))
        {
            _it.gotoBottomLeft();            
            typename OctreeClass::Iterator ittest(_it);
            ittest.gotoBottomLeft();
        }

    /**
     * \brief Gets the computed zones.
     *
     * See #_zones.
     * \return The computed zones.
     */
    const std::vector< std::vector< std::pair<int, CellClass*> > >& getZones() const {
        return _zones;
    }

    /**
     * \brief Gets the computed internal zone bounds.
     *
     * See #_internalZoneBounds.
     * \return The computed zone bounds.
     */
    const std::vector< std::vector< BoundClass > >& getZoneBounds() const {
        return _internalZoneBounds;
    }

    /**
     * \brief Gets the computed leaf zone bounds.
     *
     * See #_leafZoneBounds.
     * \return The computed zone bounds.
     */
    const std::vector< std::vector< BoundClass > >& getLeafZoneBounds() const {
        return _leafZoneBounds;
    }

    /// Gets the tree topmost level used. 
    int getTopMostLevel() const {
        return _topMostLevel;
    }

    /// Sets the tree topmost level used. 
    void setTopMostLevel(unsigned int level) {
        if( level > _treeHeight-1 ) {
            std::stringstream msgstream;
            msgstream << __FUNCTION__ << ": level is to deep. level=" << level
                      << " tree height=" << _treeHeight;
            throw std::out_of_range(msgstream.str());
        }
            
        _topMostLevel = level;
    }

    /// Gets the tree bottom most level that we use. 
    int getBottomMostLevel() const {
        return _bottomMostLevel;
    }

    /// Sets the tree bottom most level that we use. 
    void setBottomMostLevel(unsigned int level) {
        if( level > _treeHeight-1 ) {
            std::stringstream msgstream;
            msgstream << __FUNCTION__ << ": level is to deep. level=" << level
                      << " tree height=" << _treeHeight;
            throw std::out_of_range(msgstream.str());
        }
            
        _bottomMostLevel = level;
    }


    /**
     * \brief Runs the costzones algorithm.
     *
     * Since ScalFMM's implementation of the octree doesn't have a real root, we
     * run the algorithm on each of its children.
     */
    void run() {

        // Compute tree leaves total cost;
        computeLeavesCost();
        // Compute tree internal nodes total cost;
        computeInternalCost();

        // Count the root's children (the root is not stored in the tree)
        _it.gotoTop();
        _it.gotoLeft();
        do {
            this->costzones();
        } while( _it.moveRight() );

    }

protected:

    /**
     * \brief Main costzone algorithm.
     *
     * Moves through the tree in-order and assigns each cell to a zone. When a
     * zone's cumulative cost is too high, the new cells are inserted in the
     * next one.
     */
    void costzones() {
        
        std::pair<int,int> childrenCount;
        const int level = _it.level();
        const bool progressDown = _it.canProgressToDown()
            && (level < _bottomMostLevel);
        // Is this cell within the level range we consider
        const bool useCell = (level < _bottomMostLevel)
            && (level >= _topMostLevel);

        // When not on a leaf, apply to left children first
        if ( progressDown ) {
            childrenCount = countLeftRightChildren();
            callCostZonesOnChildren(LEFT, childrenCount);
        }
        
        // if the current cell is within the levels we consider, we add it
        if( useCell )
            addCurrentCell();
        
        // When not on a leaf, apply to right children
        if ( progressDown ) {
            callCostZonesOnChildren(RIGHT, childrenCount);
        }

    }

    
    /**
     * \brief Applies costzones to the left or right children of the current cell.
     *
     * The current cell is the one currently pointed at by the iterator _it.
     *
     * \warning You must check by yourself whether the cell is a leaf or not.
     *
     * \param side The children side we want to visit.
     * \param childrenCount The children count as returned by
     *                      countLeftRightChildren.
     */
    void callCostZonesOnChildren(const ChildrenSide side, const std::pair<int, int>& childrenCount) {
        
        const int& nbChildren = (side == LEFT ? childrenCount.first : childrenCount.second);

        // Don't move if there are no children on the right when we want to
        // visit them. We test this before moving in case one day moving in the
        // tree becomes expensive.
        if ( side == RIGHT && childrenCount.second == 0)
            return;

        // move down to the children level
        _it.moveDown();
        
        if ( side == RIGHT ) {
            // move to the first right child
            for ( int childIdx = 0; childIdx < childrenCount.first; childIdx++) {
                _it.moveRight();
            }
        }

        // Call costzones
        for ( int childIdx = 0; childIdx < nbChildren; childIdx++ ) {
            this->costzones();
            if(childIdx < nbChildren -1) // nbChildren-1 to avoid changing tree
                _it.moveRight();
        }

        // move up to the previous cell level
        _it.moveUp();

    }


    /**
     * \brief Adds the current cell to a zone.
     *
     * The choice of the zone is made according to the current cost accumulation
     * compared to the mean cost of a zone (_totalCost/_nbZones +1).
     *
     * This method uses the following attributes to choose the zone into which
     * the current cell must be stored :
     *
     *   - #_internalCurrentCost
     *   - #_leafCurrentCost
     *   - #_internalTotalCost
     *   - #_leafTotalCost
     *   - #_nbZones
     */
    void addCurrentCell() {

        const int& level = _it.level();
        CellClass* cell = _it.getCurrentCell();
        CostType cellCost = cell->getCost();
        bool isLeaf = (level == _treeHeight -1);

        if ( 0 == cellCost ) {
            return;
        }

        // find cell zone
        long long int cellZone = 0;

        // Near-field zone /////////////////
        if( isLeaf ) {
            CostType leafCost = cell->getNearCost();
            int colour = coord2colour(_it.getCurrentGlobalCoordinate());

            cellZone = _leafCurrentCost[colour] * _nbZones / (_leafTotalCost[colour]+1);
            _leafCurrentCost[colour] += leafCost;

            if( _leafZoneBounds.at(cellZone)[colour] == _boundInit ) {
                _leafZoneBounds.at(cellZone)[colour].first = _it;
                _leafZoneBounds.at(cellZone)[colour].second = 1;
            } else {
                _leafZoneBounds.at(cellZone)[colour].second++;
            }
        }
        ////////////////////////////////////

        // Far-field zone //////////////////
        cellZone = _internalCurrentCost * _nbZones / (_internalTotalCost+1);

        _internalCurrentCost += cellCost;

        if( _boundInit == _internalZoneBounds.at(cellZone)[level] ) {
            _internalZoneBounds.at(cellZone)[level].first = _it;
            _internalZoneBounds.at(cellZone)[level].second = 1;
        } else {
            _internalZoneBounds.at(cellZone)[level].second++;                
        }
        ///////////////////////////////

        // add cell to exhaustive zone vector
        (_zones.at(cellZone)).emplace_back(level, cell);
    }


    /**
     * \brief Computes and stores the leaves' total cost.
     *
     * The tree iterator (#_it) is moved to the bottom level of the
     * tree by this method. After the method returns, the iterator is left at
     * the rightmost leaf.
     */
    void computeLeavesCost() {

        // Reset colour costs.
        for( CostType& colourCost : _leafTotalCost ) {
            colourCost = 0;
        }

        _it.gotoBottomLeft();
        do {
            int leafColour = coord2colour(_it.getCurrentGlobalCoordinate());
            _leafTotalCost[leafColour] += _it.getCurrentCell()->getNearCost();
        } while(_it.moveRight());        

    }

    /**
     * \brief Computes and stores the internal cells' total cost.
     * \warning This method makes use of 
     *
     * The tree iterator (#_it) is moved to the bottom level of the
     * tree by this method. After the method returns, the iterator is left at
     * the rightmost leaf.
     */
    void computeInternalCost() {
        _it.gotoBottomLeft();
        //_it.moveUp();

        while( _it.level() >= _bottomMostLevel ) {
            _it.moveUp();
        }

        do {
            _it.gotoLeft();
            do {
                _internalTotalCost += _it.getCurrentCell()->getCost();
            } while(_it.moveRight());        
        } while(_it.moveUp());

    }

    /**
     * \brief Counts the left and right children of the current cell.
     *
     * The current cell is the one currently pointed at by the iterator _it.
     *
     * \warning It must be checked whether the current cell is a leaf or not
     * before calling this method.
     *
     * \return A pair of int containing the count of left (first) and right
     *         (second) children.
     */
    std::pair<int,int> countLeftRightChildren() {
        CellClass** children = _it.getCurrentChildren();
        int nbLeftChildren = 0, nbRightChildren = 0; 
        // Left children
        for ( int childIdx = 0; childIdx < 4; childIdx++) {
            if ( children[childIdx] != nullptr ) {
                ++ nbLeftChildren;
            }
        }
        // Right children
        for ( int childIdx = 4; childIdx < 8; childIdx++) {
            if ( children[childIdx] != nullptr) {
                ++ nbRightChildren;
            }
        }

        return std::pair<int,int> (nbLeftChildren, nbRightChildren);
    }



};

#endif
