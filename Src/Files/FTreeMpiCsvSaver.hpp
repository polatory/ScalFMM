#ifndef FTREEMPICSVSAVER_HPP
#define FTREEMPICSVSAVER_HPP


#include "../Utils/FGlobal.hpp"
#include "../Utils/FMpi.hpp"

#include <cstring>
#include <iostream>
#include <fstream>

/** This class is to export a tree in csv file
  *
  */
template <class OctreeClass, class ContainerClass , class ParticleClass>
class FTreeMpiCsvSaver {
    FMpi::FComm comm;           //< Communicator
    const bool includeHeader;   //< To include a line of header
    int nbFrames;               //< The current frame
    char basefile[512];         //< The base file name like "~/OUT/simulation%d.csv"

public:
    /** Constructor
      * @param inBasefile is the output file name, you must put %d in it
      */
    FTreeMpiCsvSaver(const char inBasefile[], const FMpi::FComm& communicator, const bool inIncludeHeader = false)
        : comm(communicator), includeHeader(inIncludeHeader), nbFrames(0) {
        strcpy(basefile, inBasefile);
    }

    /** Virtual destructor
      */
    virtual ~FTreeMpiCsvSaver(){
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


#endif // FTREEMPICSVSAVER_HPP
