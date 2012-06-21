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
			* @param inIncludeHeader tells if header must be included
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

            while( iter.hasNotFinished() ){
                file << iter.data().getPosition().getX() << "," << iter.data().getPosition().getY() << "," <<
                        iter.data().getPosition().getZ() << "," << getValue(&iter.data()) << "\n";
                iter.gotoNext();
            }

            const bool isUsingTsm = (octreeIterator.getCurrentListTargets() != octreeIterator.getCurrentListSrc());
            if( isUsingTsm ){
                typename ContainerClass::BasicIterator iterSources(*octreeIterator.getCurrentListSrc());

                while( iterSources.hasNotFinished() ){
                    file << iterSources.data().getPosition().getX() << "," << iterSources.data().getPosition().getY() << "," <<
                            iterSources.data().getPosition().getZ() << "," << getValue(&iterSources.data()) << "\n";
                    iterSources.gotoNext();
                }
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
