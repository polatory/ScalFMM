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
class FLightOctree {
    // The node class
    class Node {
        Node* next[8];      // Child
        const void* data;   // Data in this cell
    public:
        Node(){
            memset(next, 0, sizeof(Node*)*8);
        }

        virtual ~Node(){
            for(int idxNext = 0 ; idxNext < 8 ; ++idxNext){
                delete next[idxNext];
            }
        }

        void insert(const MortonIndex& index, const void* const cell, const int level){
            if(level){
                const int host = (index >> (3 * (level-1))) & 0x07;
                if(!next[host]){
                    next[host] = new Node();
                }
                next[host]->insert(index, cell, level - 1);
            }
            else{
                data = cell;
            }
        }

        const void* getCell(const MortonIndex& index, const int level) const {
            if(level){
                const int host = (index >> (3 * (level-1))) & 0x07;
                if(next[host]){
                    return next[host]->getCell(index, level - 1);
                }
                return 0;
            }
            else{
                return data;
            }
        }
    };

    // Tree root
    Node root;

public:
    FLightOctree(){
    }
    // Insert a cell
    void insertCell(const MortonIndex& index, const void* const cell, const int level){
        root.insert(index, cell, level);
    }
    // Retreive a cell
    const void* getCell(const MortonIndex& index, const int level) const{
        return root.getCell(index, level);
    }
};

#endif // FLIGHTOCTREE_HPP
