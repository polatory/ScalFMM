// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#ifndef FLIGHTOCTREE_HPP
#define FLIGHTOCTREE_HPP

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)This class is a light octree
* It is just a linked list with 8 pointers per node
* it is used to store small data in an octree way.
*/
template <class CellClass>
class FLightOctree {
    // The node class
    struct Node {
        Node* next[8];   // Child
        CellClass* cell; // Cell data

        Node() : cell(0) {
            memset(next, 0, sizeof(Node*)*8);
        }

        virtual ~Node(){
            for(int idxNext = 0 ; idxNext < 8 ; ++idxNext){
                delete next[idxNext];
            }
            delete cell;
        }
    };

    // Tree root
    Node root;

public:
    FLightOctree(){
    }

    // Insert a cell
    void insertCell(const MortonIndex& index, int level, CellClass*const inCell){
        Node* iter = &root;

        while(level){
            const int host = (index >> (3 * (level-1))) & 0x07;
            if(!iter->next[host]){
                iter->next[host] = new Node();
            }
            iter = iter->next[host];
            level -= 1;
        }

        iter->cell = inCell;
    }
    // Retreive a cell
    CellClass* getCell(const MortonIndex& index, int level) const{
        const Node* iter = &root;

        while(level){
            const int host = (index >> (3 * (level-1))) & 0x07;
            if(!iter->next[host]){
                return 0;
            }
            iter = iter->next[host];
            level -= 1;
        }

        return iter->cell;
    }
};

#endif // FLIGHTOCTREE_HPP
