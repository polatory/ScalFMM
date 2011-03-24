#ifndef FEXTENDPOTENTIAL_HPP
#define FEXTENDPOTENTIAL_HPP
// /!\ Please, you must read the license at the bottom of this page

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendPotential
* Please read the license
* This class is an extenssion.
* It proposes a Potential (double).
*/
class FExtendPotential {
protected:
    double potential;   //< The potential extended

public:
    /** Default constructor */
    FExtendPotential() : potential(0) {
    }

    /** Copy constructor */
    FExtendPotential(const FExtendPotential& other) : potential(other.potential) {
    }

    /** Destructor */
    virtual ~FExtendPotential(){
    }

    /** Copy operator */
    FExtendPotential& operator=(const FExtendPotential& other) {
        this->potential = other.potential;
        return *this;
    }

    /** To get the potential */
    double getPotential() const {
        return this->potential;
    }

    /** To set the potential */
    void setPotential(const double inPotential) {
        this->potential = inPotential;
    }

};


#endif //FEXTENDPOTENTIAL_HPP

// [--LICENSE--]
