// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
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
