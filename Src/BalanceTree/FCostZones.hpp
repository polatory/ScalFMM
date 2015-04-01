// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _COSTZONES_HPP_
#define _COSTZONES_HPP_

#include "FCostCell.hpp"

#include <vector>

/**
 * \brief The costzones algorithm implementation.
 * \author Quentin Khan <quentin.khan@inria.fr>
 *
 * This class is an implementation of the costzones algorithm described in "A
 * Parallel Adaptive Fast Multipole Method" (1993). The algorithm consists in an
 * in-order traversal of the octree where cell costs are accumulated. When an
 * accumulation is too big, a new zone is created.
 *
 * \tparam OctreeClass The type of the octree to work on.
 * \tparam CellClass   The type of the cells we work with.
 */
template<typename OctreeClass, typename CellClass>
class CostZones {
public:

    using BoundClass = std::pair<MortonIndex, int>;
    /// Initial value for empty bounds.
    const BoundClass _boundInit {-1,0};

protected:
    /// The iterator to move through the tree.
    typename OctreeClass::Iterator _it;
    /// The number of zones to create.
    int _nbZones;
    /// The tree height.
    int _treeHeight;

    /// The current cumulative cost of visited cells.
    unsigned long long _currentCost = 0;
    /// The total tree cost.
    unsigned long long _totalCost = 0;

    /// The vector containing the costzone cells. Sorted by zone > (level, cell).
    std::vector< std::vector< std::pair<int, CellClass*> > > _zones;

    /**
     * \brief The vector containing the Morton index boundaries of the zones.
     * Sorted by zone > level > bounds.
     * \details This means there are at most _treeHeight elements in the inner
     * vectors.
     */
    std::vector< std::vector< BoundClass > > _zonebounds;

    /// \brief Enumeration to specify the children to move to during the in-order
    /// traversal.
    enum ChildrenSide {LEFT, RIGHT};

public:

    /**
     * \brief Constructor
     * \param tree    The tree to work on.
     * \param nbZones The number of zones to create.
     */
    CostZones(OctreeClass* tree, int nbZones) : 
        _it( tree ),
        _nbZones( nbZones ),
        _treeHeight( tree->getHeight() ),
        _zones( 1, std::vector< std::pair< int, CellClass*> >( )),
        _zonebounds( 1, std::vector< BoundClass >(_treeHeight, _boundInit ))
        {}

    /**
     * \brief Gets the computed zones.
     *
     * See CostZones#_zones.
     * \return The computed zones.
     */
    const std::vector< std::vector< std::pair<int, CellClass*> > >& getZones() const {
        return _zones;
    }

    /**
     * \brief Gets the computed zone bounds.
     *
     * See CostZones#_zonebounds.
     * \return The computed zone bounds.
     */
    const std::vector< std::vector< BoundClass > >& getZoneBounds() const {
        return _zonebounds;
    }


    /**
     * \brief Runs the costzones algorithm.
     */
    void run() {
        _totalCost = 0;
        int nbRootChildren = 0;
        // Compute tree total cost;
        _it.gotoBottomLeft();
        do {
            _it.gotoLeft();
            nbRootChildren = 0; // waste of ressources only toplevel _iteration is kept
            do {
                _totalCost += _it.getCurrentCell()->getCost();
                nbRootChildren++;
            } while(_it.moveRight());        
        } while(_it.moveUp());
        
        _it.gotoLeft();
        
        // Compute costzones, we have to do the first level manualy
        for ( int i = 0; i < nbRootChildren; i++ ) {
            this->costzones();
            _it.moveRight();
        }
    }

private:

    /**
     * \brief Counts the left and right children of the current cell.
     *
     * The current cell is the one currently pointed at by the iterator _it.
     *
     * \warning You must check by yourself whether the cell is a leaf or not.
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

        // move up to the cell level
        _it.moveUp();
    }


    /**
     * \brief Adds the current cell to a zone.
     *
     * The choice of the zone is made according to the current cost accumulation.
     */
    void addCurrentCell() {
        CellClass* cell = _it.getCurrentCell();
        long int cellCost = cell->getCost();

        if ( cellCost != 0) {
            if ( _currentCost + cellCost < _zones.size() * _totalCost / _nbZones + 1 ) {
                _currentCost += cellCost;
                // Zones update
                _zones.back().push_back({_it.level(), cell});
                // Zone bounds update
                if( _zonebounds.back()[_it.level()] == _boundInit ) {
                    _zonebounds.back()[_it.level()].first = _it.getCurrentGlobalIndex();
                    _zonebounds.back()[_it.level()].second = 1;
                } else {
                    _zonebounds.back()[_it.level()].second++;
                    // This was to keep the end index
                    // _zonebounds.back()[_it.level()].second = _it.getCurrentGlobalIndex();
                }
                
            } else {
                // Add a new zone
                _zones.push_back(
                    std::vector< std::pair<int, CellClass*> >(1, {_it.level(), cell}));
                // Add a new inferior bound
                _zonebounds.push_back(
                    std::vector< BoundClass >(_treeHeight, _boundInit));
                _zonebounds.back()[_it.level()].first = _it.getCurrentGlobalIndex();
                _zonebounds.back()[_it.level()].second = 1;
            }
        }        
    }


    /**
     * \brief Main costzone algorithm.
     *
     * Moves through the tree in-order and assigns each cell to a zone. When a
     * zone's cumulative cost is too high, the new cells are insterted in the
     * next one.
     */
    void costzones() {

// DEBUG SECTION
#if 0
        CellClass* cell = _it.getCurrentCell();    
        std::cout << "in  lvl " << std::setw(2) << _it.level() << " |"
                  << "cellidx " << std::setw(4)
                  << cell->getMortonIndex() << " : "
                  << cell->getCoordinate() << " "
                  << ( _it.canProgressToDown() ? "internal" : "leaf") << " "
                  << std::endl;
//
        if (cell->_visited) {
            std::cerr << "Error : cell revisited..." << _it.level()
                      << ": " << cell->getCoordinate() << std::endl;
            throw std::exception();
            return;
        } else {
            cell->_visited = true;
        }
#endif
// END DEBUG SECTION

        std::pair<int,int> childrenCount;
        bool isNotLeaf = _it.canProgressToDown();

        // When not on a leaf, apply to left children first
        if ( isNotLeaf ) {
            childrenCount = countLeftRightChildren();
            callCostZonesOnChildren(LEFT, childrenCount);
        }

        addCurrentCell();

        // When not on a leaf, apply to right children
        if ( isNotLeaf ) {
            callCostZonesOnChildren(RIGHT, childrenCount);
        }
    }
};



#endif
