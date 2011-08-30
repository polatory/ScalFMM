#ifndef FTESTCELL_HPP
#define FTESTCELL_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "FBasicCell.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FBasicCell
* Please read the license
*
* This class is used in the FTestKernels, please
* look at this class to know what it is.
*
* This cell simply store the data when up/down
*/
class FTestCell : public FBasicCell  {
protected:
    // To store data during upward and downward pass
    long dataUp, dataDown;
public:
    FTestCell(): dataUp(0) , dataDown(0){
    }
    /** Default destructor */
    virtual ~FTestCell(){
    }
    /** When doing the upward pass */
    long getDataUp() const {
        return this->dataUp;
    }
    /** When doing the upward pass */
    void setDataUp(const long inData){
        this->dataUp = inData;
    }
    /** When doing the downard pass */
    long getDataDown() const {
        return this->dataDown;
    }
    /** When doing the downard pass */
    void setDataDown(const long inData){
        this->dataDown = inData;
    }
};


#endif //FTESTCELL_HPP

// [--LICENSE--]
