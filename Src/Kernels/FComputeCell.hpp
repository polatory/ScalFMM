#ifndef FCOMPUTECELL_HPP
#define FCOMPUTECELL_HPP
// [--License--]


#include "../Utils/FComplexe.hpp"
#include "../Utils/FMemUtils.hpp"

#include "../Components/FBasicCell.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FComputeCell
* Please read the license.
*
*/
class FComputeCell : public FBasicCell {
protected:
    static int DevP;
    static int ExpP;

    FComplexe* multipole_exp; //< For multipole extenssion
    FComplexe* local_exp;     //< For local extenssion

public:
    static void Init(const int inDevP){
        DevP = inDevP;
        ExpP = int((DevP+1) * (DevP+2) * 0.5);
    }

    static int GetP(){
        return DevP;
    }

    static int GetExp(){
        return ExpP;
    }


    /** Default constructor */
    FComputeCell()
        : multipole_exp(0), local_exp(0){
        multipole_exp = new FComplexe[ExpP];
        local_exp = new FComplexe[ExpP];
    }

    /** Constructor */
    FComputeCell(const FComputeCell& other)
        : multipole_exp(0), local_exp(0){
        multipole_exp = new FComplexe[ExpP];
        local_exp = new FComplexe[ExpP];
        (*this) = other;
    }

    /** Default destructor */
    virtual ~FComputeCell(){
        delete[] multipole_exp;
        delete[] local_exp;
    }

    /** Copy constructor */
    FComputeCell& operator=(const FComputeCell& other) {
        FMemUtils::copyall(multipole_exp, other.multipole_exp, ExpP);
        FMemUtils::copyall(local_exp, other.local_exp, ExpP);
        return *this;
    }

    /** Get Multipole */
    const FComplexe* getMultipole() const {
        return multipole_exp;
    }
    /** Get Local */
    const FComplexe* getLocal() const {
        return local_exp;
    }

    /** Get Multipole */
    FComplexe* getMultipole() {
        return multipole_exp;
    }
    /** Get Local */
    FComplexe* getLocal() {
        return local_exp;
    }
};


int FComputeCell::DevP(-1);
int FComputeCell::ExpP(-1);


#endif //FCOMPUTECELL_HPP

// [--END--]
