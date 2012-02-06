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
#ifndef FTREECSVSAVER_HPP
#define FTREECSVSAVER_HPP

#include "../Utils/FGlobal.hpp"

#include <cstring>
#include <iostream>
#include <fstream>

/** This class is to export a tree in csv file
  *
  */
template <class OctreeClass, class ContainerClass , class ParticleClass>
class FTreeCsvSaver {
    const bool includeHeader;   //< To include a line of header
    int nbFrames;               //< The current frame
    char basefile[512];         //< The base file name like "~/OUT/simulation%d.csv"

public:
    /** Constructor
      * @param inBasefile is the output file name, you must put %d in it
      */
    FTreeCsvSaver(const char inBasefile[], const bool inIncludeHeader = false)
        : includeHeader(inIncludeHeader), nbFrames(0) {
        strcpy(basefile, inBasefile);
    }

    /** Virtual destructor
      */
    virtual ~FTreeCsvSaver(){
    }

    /** to know how many frame has been saved
      */
    int getNbFrames() const {
        return nbFrames;
    }

    /** export a tree
      */
    void exportTree(OctreeClass*const tree){
        char currentFilename[512];
        sprintf(currentFilename, basefile, nbFrames++);

        std::ofstream file(currentFilename, std::ofstream::out );
        if(includeHeader){
            file << "x, y, z, value\n";
        }

        typename OctreeClass::Iterator octreeIterator(tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
            //const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());

            while( iter.hasNotFinished() ){
                file << iter.data().getPosition().getX() << "," << iter.data().getPosition().getY() << "," <<
                        iter.data().getPosition().getZ() << "," << getValue(&iter.data()) << "\n";
                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());

        file.close();
    }

    /** Inherit from this class and customize this function if you need it
      */
    virtual FReal getValue(ParticleClass*const part){
        return FReal(0);
    }
};


#endif // FTREECSVSAVER_HPP
