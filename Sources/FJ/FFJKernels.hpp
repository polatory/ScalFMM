#ifndef FFJKERNELS_HPP
#define FFJKERNELS_HPP
// /!\ Please, you must read the license at the bottom of this page

#include "../Core/FAbstractKernels.hpp"
#include "../Containers/FList.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FFJKernels
* @brief
* Please read the license
*
*/
template< class ParticuleClass, class CellClass>
class FFJKernels : public FAbstractKernels<ParticuleClass,CellClass> {
public:
    /** Default destructor */
    virtual ~FFJKernels(){
    }

    /** When init the kernel */
    void init(){
    }

    Complex[] getSCoeff(Complex xi, Complex xstar) {
        Complex[] ans = new Complex[p];
        ans[0] = new Complex(1d);
        for (int i=1; i<p; i++) {
            ans[i] = xi.subtract(xstar).pow(i).divide(i).negate();
        }
        return ans;
    }

    int getIndex(Point[] z, Point p) {
        int ans = -1;
        for (int i=0; i<z.length; i++)
            if (z[i].equals(p))
                return i;
        return ans;
    }

    /** Print the number of particules */
    void P2M(CellClass* const pole, FList<ParticuleClass*>* const particules) {
        for(typename FList<ParticuleClass*>::BasicIterator iterParticule(*inParticules);
                                iterParticule.isValide() ; iterParticule.progress()){
            //Point thisX = xPoints[j];
            Complex[] B = getSCoeff(iterParticule.value()->getPosition(), inPole->getPosition());
            // compute
            const double thisU = u[getIndex(x, thisX)];
            for( int k = 0 ; k < B.length; ++k ){
                pole->setMultipole(k, B[k].multiply(thisU));//add?
            }
        }
    }

    /** Print the morton index */
    void M2M(CellClass* const pole, CellClass** const child, const int inLevel) {
    }

    /** Print the morton index */
    void M2L(CellClass* const pole, CellClass** const distantNeighbors, const int size, const int inLevel) {
    }

    /** Print the morton index */
    void L2L(CellClass* const local, CellClass** const child, const int inLevel) {
    }

    /** Print the number of particules */
    void L2P(CellClass* const local, FList<ParticuleClass*>* const particules){
    }

    /** Print the number of particules */
    void P2P(FList<ParticuleClass*>* const currentBox, FList<ParticuleClass*>** directNeighbors, const int size) {
    }
};


#endif //FFJKERNELS_HPP

// [--LICENSE--]
