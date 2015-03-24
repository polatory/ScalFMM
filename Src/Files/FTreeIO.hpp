// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FTREEIO_HPP
#define FTREEIO_HPP

#include <iostream>
#include <fstream>

#include "../Containers/FBufferReader.hpp"
#include "../Containers/FBufferWriter.hpp"


/** This class proposes static methods to save and load
  * a tree.
  * It used binary format (so FReal must be the same!)
  *
  * The format is :
  * [nb leaves]
  * [nb particles in leaf][particles.....] ...
  * [nb cells at level]
  * [morton index][cell] ....
  * ...
  */
template <class FReal>
class FTreeIO{

    enum TargetSourceDiff{
        TsmUndef,
        TsmUsed,
        TsmUnused
    };

public:
    /** To save in memory */
    template <class OctreeClass, class CellClass, class LeafClass, class ContainerClass >
    static bool Save(const char filename[], OctreeClass& tree){
        std::ofstream file(filename, std::ofstream::binary | std::ofstream::out );
        FBufferWriter buffer;

        if(!file.good()){
            return false;
        }

        const size_t sizeof_freal = sizeof(FReal);
        file.write((const char*)&sizeof_freal, sizeof(size_t));

        const int height = tree.getHeight();
        file.write((const char*)&height, sizeof(int));

        const int subHeight = tree.getSubHeight();
        file.write((const char*)&subHeight, sizeof(int));

        const FReal width = tree.getBoxWidth();
        file.write((const char*)&width, sizeof(FReal));

        file.write((const char*)&tree.getBoxCenter(), sizeof(FPoint<FReal>));

        {
            typename OctreeClass::Iterator octreeIterator(&tree);

            int nbLeaf = 0;
            const std::ofstream::pos_type posNbLeaf = file.tellp();
            file.write((const char*)&nbLeaf,sizeof(int));

            octreeIterator.gotoBottomLeft();
            const bool useTargetSource = (octreeIterator.getCurrentListSrc() != octreeIterator.getCurrentListTargets());
            if( useTargetSource ){
                TargetSourceDiff tsm = TsmUsed;
                file.write((const char*)&tsm,sizeof(tsm));

                do{
                    const MortonIndex mindex = octreeIterator.getCurrentGlobalIndex();
                    file.write((const char*)&mindex,sizeof(mindex));

                    buffer.reset();
                    octreeIterator.getCurrentListSrc()->save(buffer);
                    octreeIterator.getCurrentListTargets()->save(buffer);

                    const int sizeOfLeaf = buffer.getSize();
                    file.write((const char*) &sizeOfLeaf, sizeof(int));
                    file.write(buffer.data(), buffer.getSize());

                    ++nbLeaf;
                } while(octreeIterator.moveRight());
            }
            else{
                TargetSourceDiff tsm = TsmUnused;
                file.write((const char*)&tsm,sizeof(tsm));
                do{
                    const MortonIndex mindex = octreeIterator.getCurrentGlobalIndex();
                    file.write((const char*)&mindex,sizeof(mindex));

                    buffer.reset();
                    octreeIterator.getCurrentListSrc()->save(buffer);

                    const int sizeOfLeaf = buffer.getSize();
                    file.write((const char*) &sizeOfLeaf, sizeof(int));
                    file.write(buffer.data(), buffer.getSize());

                    ++nbLeaf;
                } while(octreeIterator.moveRight());
            }

            const std::ofstream::pos_type currentPos = file.tellp();
            file.seekp(posNbLeaf);
            file.write((const char*)&nbLeaf,sizeof(int));
            file.seekp(currentPos);
        }

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = tree.getHeight() - 1 ; idxLevel > 1 ; --idxLevel ){
            int nbCells = 0;
            const std::ofstream::pos_type posNbCells = file.tellp();
            file.write((const char*)&nbCells,sizeof(int));

            do{
                const MortonIndex mindex = octreeIterator.getCurrentGlobalIndex();
                file.write((const char*)&mindex,sizeof(MortonIndex));

                buffer.reset();
                octreeIterator.getCurrentCell()->save(buffer);
                const int sizeOfCell = buffer.getSize();
                file.write((const char*) &sizeOfCell, sizeof(int));
                file.write(buffer.data(), buffer.getSize());

                ++nbCells;
            } while(octreeIterator.moveRight());

            const std::ofstream::pos_type currentPos = file.tellp();
            file.seekp(posNbCells);
            file.write((const char*)&nbCells,sizeof(int));
            file.seekp(currentPos);

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;// equal octreeIterator.moveUp(); octreeIterator.gotoLeft();
        }


        file.flush();
        file.close();

        return true;
    }


    /** To load from memory */
    template <class OctreeClass, class CellClass, class LeafClass, class ContainerClass  >
    static bool Load(const char filename[], OctreeClass& tree){
        std::ifstream file(filename, std::ifstream::binary | std::ifstream::in );
        FBufferReader buffer;

        if(!file.good()){
            return false;
        }

        size_t sizeof_freal = 0;
        file.read((char*)&sizeof_freal, sizeof(size_t));
        if( sizeof_freal != sizeof(FReal)){
            std::cerr << "[Loader] Error Freal do not coincide with file type:\n";
            std::cerr << "[Loader] In file : " << sizeof_freal << " Real : " << sizeof(FReal) << "\n";
            return false;
        }

        int treeHeight = 0;
        file.read((char*)&treeHeight, sizeof(int));
        int treeSubHeight = 0;
        file.read((char*)&treeSubHeight, sizeof(int));
        FReal boxWidth = 0;
        file.read((char*)&boxWidth, sizeof(FReal));
        FPoint<FReal> center;
        file.read((char*)&center, sizeof(FPoint<FReal>));

        tree.~OctreeClass();
        new (&tree) OctreeClass(treeHeight,treeSubHeight,boxWidth,center);

        {
            int nbLeaf = 0;
            file.read((char*)&nbLeaf, sizeof(int));

            TargetSourceDiff tsm = TsmUndef;
            file.read((char*)&tsm,sizeof(tsm));

            for(int idxLeaf = 0 ; idxLeaf < nbLeaf ; ++idxLeaf){
                MortonIndex mindex = 0;
                file.read((char*)&mindex, sizeof(mindex));

                int sizeOfLeaf = 0;
                file.read((char*)&sizeOfLeaf, sizeof(int));

                buffer.reserve(sizeOfLeaf);
                file.read((char*)buffer.data(), sizeOfLeaf);

                LeafClass*const leaf = tree.createLeaf(mindex);
                leaf->getSrc()->restore(buffer);

                if( tsm == TsmUsed ){
                    leaf->getTargets()->restore(buffer);
                }
            }
        }

        // Start from leal level - 1
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

        // for each levels
        for(int idxLevel = tree.getHeight() - 1 ; idxLevel > 1 ; --idxLevel ){
            int nbCells = 0;
            file.read((char*)&nbCells, sizeof(int));

            do{
                MortonIndex mindex;
                file.read((char*)&mindex, sizeof(MortonIndex));
                if(mindex != octreeIterator.getCurrentGlobalIndex()){
                    std::cerr << "[Loader] Error indexes are different\n";
                    return false;
                }

                int sizeOfCell = 0;
                file.read((char*)&sizeOfCell, sizeof(int));

                buffer.reserve(sizeOfCell);
                file.read((char*)buffer.data(), sizeOfCell);

                octreeIterator.getCurrentCell()->restore(buffer);

                --nbCells;
            } while(octreeIterator.moveRight());

            if(nbCells != 0){
                std::cerr << "[Loader] Wrong number of cells at level " << idxLevel << "\n";
                return false;
            }

            avoidGotoLeftIterator.moveUp();
            octreeIterator = avoidGotoLeftIterator;
        }

        file.close();

        return true;
    }

};


#endif // FTREEIO_HPP
