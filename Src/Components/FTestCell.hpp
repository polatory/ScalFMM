#ifndef FTESTCELL_HPP
#define FTESTCELL_HPP
// [--License--]

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
};


#endif //FTESTCELL_HPP

// [--END--]
