#ifndef _COSTZONES_HPP_
#define _COSTZONES_HPP_

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
     */
    void countLeftRightChildren(int& nbLeftChildren, int& nbRightChildren) {
        FCostCell** children = _it.getCurrentChildren();
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
    }

    void callCostZonesOnChildren(ChildrenSide side, int nbLeftChildren, int nbRightChildren) {
        int startIdx = (side == LEFT ? 0 : 4);
        int endIdx   = (side == LEFT ? 4 : 8);
        
        int nbChildren = (side == LEFT ? nbLeftChildren : nbRightChildren);


        _it.moveDown();
        
        // move to the first right child
        if ( side == RIGHT ) {
            if ( nbRightChildren == 0)
                return;
            
        }


        _it.moveUp();
    }

    void costzones() {
            // Current position is a leaf
            CellClass* cell = _it.getCurrentCell();    

            if (cell->_visited) {
                std::cerr << "Error : cell revisited..." << _it.level()
                          << ": " << cell->getCoordinate() << std::endl;
                return;
            }
            else
                cell->_visited = true;

            int cellCost = cell->getCost();

            int nbLeftChildren  = 0;
            int nbRightChildren = 0;

#if 0
            std::cout << "lvl " << std::setw(2) << _it.level() << " |"
                      << "cellidx " << std::setw(4)
                      << cell->getMortonIndex() << " : "
                      << cell->getCoordinate() << " "
                      << ( _it.canProgressToDown() ? "internal" : "leaf") << " "
                      << std::endl;
#endif
            // When not on a leaf, apply to left children first
            if ( _it.canProgressToDown() ) {

                FCostCell** children = _it.getCurrentChildren();
                // Count left children
                for ( int childIdx = 0; childIdx < 4; childIdx++) {
                    if ( children[childIdx] ) {
                        ++ nbLeftChildren;
                    }
                }
        
                _it.moveDown();
                // Apply costzones to left children
                for ( int childIdx = 0; childIdx < nbLeftChildren; childIdx++ ) {
                    costzones();
                    // avoid changing tree
                    if(childIdx < nbLeftChildren -1) 
                        _it.moveRight();
                }
                _it.moveUp();
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

                FCostCell** children = _it.getCurrentChildren();
                // Count right children
                for ( int childIdx = 4; childIdx < 8; childIdx++) {
                    if ( children[childIdx] ) {
                        ++ nbRightChildren;
                    }
                }
                
                if ( nbRightChildren == 0)
                    return;

                // Move to the first right child
                _it.moveDown();
                for ( int childIdx = 0; childIdx < nbLeftChildren; childIdx++) {
                    _it.moveRight();
                }
                // Apply costzones to the right children
                for ( int childIdx = 0; childIdx < nbRightChildren; childIdx++) {
                    costzones();
                    // avoid changing tree
                     if(childIdx < nbRightChildren -1)
                        _it.moveRight();
                }
                _it.moveUp();
            }
        }
};



#endif
