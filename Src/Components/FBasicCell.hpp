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
#ifndef FBASICCELL_HPP
#define FBASICCELL_HPP


#include "../Extensions/FExtendMortonIndex.hpp"
#include "../Extensions/FExtendCoordinate.hpp"



/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicCell
*
* This class defines a basic cell used for examples. It extends
* the mininum, only what is needed by FOctree and FFmmAlgorithm
* to make the things working.
* By using this extension it will implement the FAbstractCell without
* inheriting from it.
*
*
*/
class FBasicCell : public FExtendMortonIndex, public FExtendCoordinate {
public:
    /** Default destructor */
    virtual ~FBasicCell(){
    }

    /** Save the current cell in a buffer */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        FExtendMortonIndex::save(buffer);
        FExtendCoordinate::save(buffer);
    }
    /** Restore the current cell from a buffer */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        FExtendMortonIndex::restore(buffer);
        FExtendCoordinate::restore(buffer);
    }

    /** Do nothing */
    void resetToInitialState(){
    }
};


#endif //FBASICCELL_HPP


