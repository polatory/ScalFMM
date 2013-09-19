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
#ifndef FTESTCELL_HPP
#define FTESTCELL_HPP


#include "FBasicCell.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicCell*
* @brief This class is used in the FTestKernels, please look at this class to know how to customize a cell.
*
* This cell simply store the data when up/down.
* It also shows how to be restored and saved, etc.
*/
class FTestCell : public FBasicCell  {
protected:
    // To store data during upward and downward pass
    long long int dataUp, dataDown;
public:
    FTestCell(): dataUp(0) , dataDown(0){
    }
    /** Default destructor */
    virtual ~FTestCell(){
    }
    /** When doing the upward pass */
    long long int getDataUp() const {
        return this->dataUp;
    }
    /** When doing the upward pass */
    void setDataUp(const long long int inData){
        this->dataUp = inData;
    }
    /** When doing the downard pass */
    long long int getDataDown() const {
        return this->dataDown;
    }
    /** When doing the downard pass */
    void setDataDown(const long long int inData){
        this->dataDown = inData;
    }

    /////////////////////////////////////////////////

    /** Save the current cell in a buffer */
    void save(FBufferWriter& buffer) const{
        FBasicCell::save(buffer);
        buffer << dataDown << dataUp;
    }

    /** Restore the current cell from a buffer */
    void restore(FBufferReader& buffer){
        FBasicCell::restore(buffer);
        buffer >> dataDown >> dataUp;
    }

    /////////////////////////////////////////////////

    /** Serialize only up data in a buffer */
    void serializeUp(FBufferWriter& buffer) const {
        buffer << this->dataUp;
    }
    /** Deserialize only up data in a buffer */
    void deserializeUp(FBufferReader& buffer){
        buffer >> this->dataUp;
    }

    /** Serialize only down data in a buffer */
    void serializeDown(FBufferWriter& buffer) const {
        buffer << this->dataDown;
    }
    /** Deserialize only up data in a buffer */
    void deserializeDown(FBufferReader& buffer){
        buffer >> this->dataDown;
    }
};


#endif //FTESTCELL_HPP


