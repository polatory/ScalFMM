#ifndef _COSTZONES_HPP_
#define _COSTZONES_HPP_

//#include <utility>
#include <stdexcept>

template<typename OctreeClass, typename CellClass>
class CostZones {
    unsigned long long _currentCost = 0;
    unsigned long long _totalCost = 0;
    std::vector< std::pair< int, CellClass*> > _emptyzone;
    std::vector< std::vector< std::pair<int, CellClass*> > > _zones;

    typename OctreeClass::Iterator _it;
    int _nbZones;

    enum ChildrenSide {LEFT, RIGHT};

public:

    CostZones(OctreeClass* tree, int nbZones) : 
        _zones(1, _emptyzone),
        _it(tree),
        _nbZones(nbZones)
        {}

    std::vector< std::vector< std::pair<int, CellClass*> > >& getZones() {return _zones;}

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
            costzones();
            _it.moveRight();
        }
    }

private:

    /**
     * Counts the children to the left and to the right of the cell currently
     * pointed to by the iterator _it. You must check by yourself whether the
     * cell is a leaf or not.
     *
     * \return A pair of int containing the count of left (first) and right
     *         (second) children.
     */
    std::pair<int,int> countLeftRightChildren() {
        FCostCell** children = _it.getCurrentChildren();
        int nbLeftChildren = 0, nbRightChildren = 0; 
        // Left children
        for ( int childIdx = 0; childIdx < 4; childIdx++) {
            if ( children[childIdx] ) {
                ++ nbLeftChildren;
            }
        }
        // Right children
        for ( int childIdx = 4; childIdx < 8; childIdx++) {
            if ( children[childIdx] ) {
                ++ nbRightChildren;
            }
        }

        return std::pair<int,int> (nbLeftChildren, nbRightChildren);
    }


    /**
     * Calls the costzones function on the left or right children of the current
     * cell.
     *
     * \param side The children side we want to visit.
     * \param childrenCount The children count as returned by
     *                      countLeftRightChildren
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
            costzones();
            if(childIdx < nbChildren -1) // nbChildren-1 to avoid changing tree
                _it.moveRight();
        }

        // move up to the cell level
        _it.moveUp();
    }


    void costzones() {

        CellClass* cell = _it.getCurrentCell();    

// DEBUG SECTION
#if 0
        std::cout << "in  lvl " << std::setw(2) << _it.level() << " |"
                  << "cellidx " << std::setw(4)
                  << cell->getMortonIndex() << " : "
                  << cell->getCoordinate() << " "
                  << ( _it.canProgressToDown() ? "internal" : "leaf") << " "
                  << std::endl;
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

        int cellCost = cell->getCost();
        std::pair<int,int> childrenCount;

        // When not on a leaf, apply to left children first
        if ( _it.canProgressToDown() ) {
            childrenCount = countLeftRightChildren();
            callCostZonesOnChildren(LEFT, childrenCount);
        }

        if ( cellCost != 0) {
            // Add the current cell
            if ( _currentCost + cellCost < _zones.size() * _totalCost / _nbZones + 1 ) {
                _currentCost += cellCost;
                _zones.back().push_back({_it.level(), cell});
            } else {
                _zones.push_back(std::vector< std::pair<int, CellClass*> >(1, {_it.level(), cell}));
            }
        }
    
        // When not on a leaf, apply to right children
        if ( _it.canProgressToDown() ) {
            callCostZonesOnChildren(RIGHT, childrenCount);
        }
    }
};



#endif
