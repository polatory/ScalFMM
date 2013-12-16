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
#ifndef FABSTRACTSPHERICALKERNEL_HPP
#define FABSTRACTSPHERICALKERNEL_HPP

#include "../../Components/FAbstractKernels.hpp"

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FTrace.hpp"
#include "../../Utils/FMemUtils.hpp"
#include "../../Utils/FSmartPointer.hpp"
#include "../../Utils/FPoint.hpp"
#include "../../Utils/FAssert.hpp"

#include "../../Containers/FTreeCoordinate.hpp"

#include "../P2P/FP2P.hpp"

#include "FHarmonic.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* This is the abstract spherical harmonic kernel
*/
template< class CellClass, class ContainerClass>
class FAbstractSphericalKernel : public FAbstractKernels<CellClass,ContainerClass> {
protected:
    const int   devP;           //< The P
    const FReal boxWidth;       //< the box width at leaf level
    const int   treeHeight;     //< The height of the tree

    const FReal widthAtLeafLevel;       //< the width of a box at leaf level
    const FReal widthAtLeafLevelDiv2;   //< the width of a box at leaf level divided by 2
    const FPoint boxCorner;        //< the corner of the box system

    FHarmonic harmonic; //< The harmonic computation class

    // For normal computation
    FSmartPointer<FComplexe*> preL2LTransitions; //< The pre-computation for the L2L based on the level
    FSmartPointer<FComplexe*> preM2MTransitions; //< The pre-computation for the M2M based on the level


    /** Alloc and init pre-vectors*/
    void allocAndInit(){
        preL2LTransitions = new FComplexe*[treeHeight ];
        memset(preL2LTransitions.getPtr(), 0, (treeHeight) * sizeof(FComplexe*));
        preM2MTransitions = new FComplexe*[treeHeight];
        memset(preM2MTransitions.getPtr(), 0, (treeHeight) * sizeof(FComplexe*));

        FReal treeWidthAtLevel = (boxWidth)/2;
        for(int idxLevel = 0 ; idxLevel < treeHeight - 1 ; ++idxLevel ){
            preL2LTransitions[idxLevel] = new FComplexe[ 8 * harmonic.getExpSize()];
            preM2MTransitions[idxLevel] = new FComplexe[ 8 * harmonic.getExpSize()];

            const FPoint father(treeWidthAtLevel,treeWidthAtLevel,treeWidthAtLevel);
            treeWidthAtLevel /= 2;

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild ){
                FTreeCoordinate childBox;
                childBox.setPositionFromMorton(idxChild,1);

                const FPoint M2MVector (
                        father.getX() - (treeWidthAtLevel * FReal(1 + (childBox.getX() * 2))),
                        father.getY() - (treeWidthAtLevel * FReal(1 + (childBox.getY() * 2))),
                        father.getZ() - (treeWidthAtLevel * FReal(1 + (childBox.getZ() * 2)))
                        );

                harmonic.computeInner(FSpherical(M2MVector));
                FMemUtils::copyall<FComplexe>(&preM2MTransitions[idxLevel][harmonic.getExpSize() * idxChild], harmonic.result(), harmonic.getExpSize());

                const FPoint L2LVector (
                        (treeWidthAtLevel * FReal(1 + (childBox.getX() * 2))) - father.getX(),
                        (treeWidthAtLevel * FReal(1 + (childBox.getY() * 2))) - father.getY(),
                        (treeWidthAtLevel * FReal(1 + (childBox.getZ() * 2))) - father.getZ()
                        );

                harmonic.computeInner(FSpherical(L2LVector));
                FMemUtils::copyall<FComplexe>(&preL2LTransitions[idxLevel][harmonic.getExpSize() * idxChild], harmonic.result(), harmonic.getExpSize());
           }
        }
    }

    /** Get a leaf real position from its tree coordinate */
    FPoint getLeafCenter(const FTreeCoordinate coordinate) const {
        return FPoint(
                    FReal(coordinate.getX()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getX(),
                    FReal(coordinate.getY()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getX(),
                    FReal(coordinate.getZ()) * widthAtLeafLevel + widthAtLeafLevelDiv2 + boxCorner.getX());
    }


public:
    /** Kernel constructor */
    FAbstractSphericalKernel(const int inDevP, const int inTreeHeight, const FReal inBoxWidth, const FPoint& inBoxCenter)
        : devP(inDevP),
          boxWidth(inBoxWidth),
          treeHeight(inTreeHeight),
          widthAtLeafLevel(inBoxWidth/FReal(1 << (inTreeHeight-1))),
          widthAtLeafLevelDiv2(widthAtLeafLevel/2),
          boxCorner(inBoxCenter.getX()-(inBoxWidth/2),inBoxCenter.getY()-(inBoxWidth/2),inBoxCenter.getZ()-(inBoxWidth/2)),
          harmonic(inDevP),
          preL2LTransitions(0),
          preM2MTransitions(0) {

        allocAndInit();
    }

    /** Copy constructor */
    FAbstractSphericalKernel(const FAbstractSphericalKernel& other)
        : devP(other.devP),
          boxWidth(other.boxWidth),
          treeHeight(other.treeHeight),
          widthAtLeafLevel(other.widthAtLeafLevel),
          widthAtLeafLevelDiv2(other.widthAtLeafLevelDiv2),
          boxCorner(other.boxCorner),
          harmonic(other.devP),
          preL2LTransitions(other.preL2LTransitions),
          preM2MTransitions(other.preM2MTransitions) {

    }

    /** Default destructor */
    virtual ~FAbstractSphericalKernel(){
        if(preL2LTransitions.isLast()){
            FMemUtils::DeleteAll(preL2LTransitions.getPtr(), treeHeight);
        }
        if(preM2MTransitions.isLast()){
            FMemUtils::DeleteAll(preM2MTransitions.getPtr(), treeHeight);
        }
    }

    /** P2M with a cell and all its particles */
    void P2M(CellClass* const inPole, const ContainerClass* const inParticles) {
        FComplexe* FRestrict const cellMultiPole = inPole->getMultipole();
        // Copying the position is faster than using cell position
        const FPoint polePosition = getLeafCenter(inPole->getCoordinate());
        // For all particles in the leaf box
        const FReal*const physicalValues = inParticles->getPhysicalValues();
        const FReal*const positionsX = inParticles->getPositions()[0];
        const FReal*const positionsY = inParticles->getPositions()[1];
        const FReal*const positionsZ = inParticles->getPositions()[2];
        for(int idxPart = 0 ; idxPart < inParticles->getNbParticles() ; ++idxPart){
            // P2M
            particleToMultiPole(cellMultiPole, polePosition,
                                FPoint(positionsX[idxPart],positionsY[idxPart],positionsZ[idxPart]),
                                physicalValues[idxPart]);
        }
    }

    /** M2M with a cell and all its child */
    void M2M(CellClass* const FRestrict inPole, const CellClass *const FRestrict *const FRestrict inChild, const int inLevel) {
        FComplexe* FRestrict const multipole_exp_target = inPole->getMultipole();
        // iter on each child and process M2M
        const FComplexe* FRestrict const preM2MTransitionsAtLevel = preM2MTransitions[inLevel];
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            if(inChild[idxChild]){
                multipoleToMultipole(multipole_exp_target, inChild[idxChild]->getMultipole(), &preM2MTransitionsAtLevel[idxChild * harmonic.getExpSize()]);
            }
        }
    }

    /** M2L with a cell and all the existing neighbors */
    virtual void M2L(CellClass* const FRestrict pole, const CellClass* distantNeighbors[343],
             const int size, const int inLevel) = 0;

    /** L2L with a cell and all its child */
    void L2L(const CellClass* const FRestrict pole, CellClass* FRestrict *const FRestrict child, const int inLevel) {
        // iter on each child and process L2L
        const FComplexe* FRestrict const preL2LTransitionsAtLevel = preL2LTransitions[inLevel];
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            if(child[idxChild]){
                localToLocal(child[idxChild]->getLocal(), pole->getLocal(), &preL2LTransitionsAtLevel[idxChild * harmonic.getExpSize()]);
            }
        }
    }

    /** L2P with a cell and all its particles */
    void L2P(const CellClass* const local, ContainerClass* const inParticles){
        const FComplexe* const cellLocal = local->getLocal();
        // Copying the position is faster than using cell position
        const FPoint localPosition = getLeafCenter(local->getCoordinate());
        // For all particles in the leaf box
        const FReal*const physicalValues = inParticles->getPhysicalValues();
        const FReal*const positionsX = inParticles->getPositions()[0];
        const FReal*const positionsY = inParticles->getPositions()[1];
        const FReal*const positionsZ = inParticles->getPositions()[2];
        FReal*const potentials = inParticles->getPotentials();
        FReal*const forcesX = inParticles->getForcesX();
        FReal*const forcesY = inParticles->getForcesY();
        FReal*const forcesZ = inParticles->getForcesZ();
        for(int idxPart = 0 ; idxPart < inParticles->getNbParticles() ; ++idxPart){
            // L2P
            localToParticle(localPosition, cellLocal,
                            FPoint(positionsX[idxPart],positionsY[idxPart],positionsZ[idxPart]),
                            physicalValues[idxPart], &potentials[idxPart],
                            &forcesX[idxPart],&forcesY[idxPart],&forcesZ[idxPart]);
        }
    }

    /** This P2P has to be used when target != sources
      * It will proceed an direct interation no mutual
      *
      * It takes all the target particles from the current leaf,
      * then it computes the sources/targets interaction in this leaf,
      * then it computes the sources/targets inteactions between this leaf and the
      * neighbors.
      */
    void P2P(const FTreeCoordinate& inLeafPosition,
                  ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
                  ContainerClass* const directNeighborsParticles[27], const int /*size*/){
        FP2P::FullMutual(targets,directNeighborsParticles,14);
    }

    /** This P2P has to be used when target != sources
      * It will proceed an direct interation no mutual
      *
      * It takes all the target particles from the current leaf,
      * then it computes the sources/targets interaction in this leaf,
      * then it computes the sources/targets inteactions between this leaf and the
      * neighbors.
      */
    void P2PRemote(const FTreeCoordinate& ,
                  ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict ,
                  ContainerClass* const directNeighborsParticles[27], const int /*size*/){
        FP2P::FullRemote(targets,directNeighborsParticles,27);
    }

private:


    ///////////////////////////////////////////////////////////////////////////////
    //                                  Computation
    ///////////////////////////////////////////////////////////////////////////////


    /** P2M computation
    * expansion_P2M_add
    * Multipole expansion with m charges q_i in Q_i=(rho_i, alpha_i, beta_i)
    *whose relative coordinates according to *p_center are:
    *Q_i - *p_center = (rho'_i, alpha'_i, beta'_i);
    *
    *For j=0..P, k=-j..j, we have:
    *
    *M_j^k = (-1)^j { sum{i=1..m} q_i Inner_j^k(rho'_i, alpha'_i, beta'_i) }
    *
    *However the extern loop is over the bodies (i=1..m) in our code and as an
    *intern loop we have: j=0..P, k=-j..j
    *
    *and the potential is then given by:
    *
    * Phi(x) = sum_{n=0}^{+} sum_{m=-n}^{n} M_n^m O_n^{-m} (x - *p_center)
    *
    */
    void particleToMultiPole(FComplexe* const cellMultiPole, const FPoint& inPolePosition ,
                             const FPoint& particlePosition, const FReal particlePhysicalValue){

        // Inner of Qi - Z0 => harmonic.result
        harmonic.computeInner( FSpherical(particlePosition - inPolePosition) );

        FReal minus_one_pow_j = 1.0;    // (-1)^j => be in turn 1 and -1
        const FReal qParticle = particlePhysicalValue; // q in the formula
        int index_j_k = 0; // p_exp_term & p_Y_term

        // J from 0 to P
        for(int j = 0 ; j <= devP ; ++j){
            // k from 0 to J
            for(int k = 0 ; k <= j ; ++k, ++index_j_k){
                harmonic.result(index_j_k).mulRealAndImag( qParticle * minus_one_pow_j );
                cellMultiPole[index_j_k] += harmonic.result(index_j_k);
            }
            // (-1)^J => -1 becomes 1 or 1 becomes -1
            minus_one_pow_j = -minus_one_pow_j;
        }
    }

    /* M2M
    *We compute the translation of multipole_exp_src from *p_center_of_exp_src to
    *p_center_of_exp_target, and add the result to multipole_exp_target.
    *
    * O_n^l (with n=0..P, l=-n..n) being the former multipole expansion terms
    * (whose center is *p_center_of_multipole_exp_src) we have for the new multipole
    * expansion terms (whose center is *p_center_of_multipole_exp_target):

    * M_j^k = sum{n=0..j}
    * sum{l=-n..n, |k-l|<=j-n}
    * O_n^l Inner_{j-n}^{k-l}(rho, alpha, beta)
    *
    * where (rho, alpha, beta) are the spherical coordinates of the vector :
    * p_center_of_multipole_exp_target - *p_center_of_multipole_exp_src
    *
    * Warning: if j-n < |k-l| we do nothing.
     */
    void multipoleToMultipole(FComplexe* const FRestrict multipole_exp_target,
                              const FComplexe* const FRestrict multipole_exp_src,
                              const FComplexe* const FRestrict M2M_Inner_transfer){

        // n from 0 to P
        for(int n = 0 ; n <= devP ; ++n ){
            // l<0 // (-1)^l
            FReal pow_of_minus_1_for_l = ( n & 1 ? FReal(-1.0) : FReal(1.0) );

            // O_n^l : here points on the source multipole expansion term of degree n and order |l|
            const int index_n = harmonic.getPreExpRedirJ(n);

            // l from -n to <0
            for(int l = -n ; l < 0 ; ++l){
                const FComplexe M_n__n_l = multipole_exp_src[index_n -l];

                // j from n to P
                for(int j = n ; j <= devP ; ++j ){
                    // M_j^k
                    const int index_j = harmonic.getPreExpRedirJ(j);
                    // Inner_{j-n}^{k-l} : here points on the M2M transfer function/expansion term of degree n-j and order |k-l|
                    const int index_j_n = harmonic.getPreExpRedirJ(j-n); /* k==0 */

                    // since n-j+l<0
                    for(int k = 0 ; k <= (j-n+l) ; ++k ){ // l<0 && k>=0 => k-l>0
                        const FComplexe I_j_n__k_l = M2M_Inner_transfer[index_j_n + k - l];

                        multipole_exp_target[index_j + k].incReal( pow_of_minus_1_for_l *
                                                    ((M_n__n_l.getReal() * I_j_n__k_l.getReal()) +
                                                     (M_n__n_l.getImag() * I_j_n__k_l.getImag())));
                        multipole_exp_target[index_j + k].incImag( pow_of_minus_1_for_l *
                                                    ((M_n__n_l.getReal() * I_j_n__k_l.getImag()) -
                                                     (M_n__n_l.getImag() * I_j_n__k_l.getReal())));

                     } // for k
                } // for j

                pow_of_minus_1_for_l = -pow_of_minus_1_for_l;
            } // for l

            // l from 0 to n
            for(int l = 0 ; l <= n ; ++l){
                const FComplexe M_n__n_l = multipole_exp_src[index_n + l];

                // j from n to P
                for( int j = n ; j <= devP ; ++j ){
                    const int first_k = FMath::Max(0,n-j+l);
                    // (-1)^k
                    FReal pow_of_minus_1_for_k = static_cast<FReal>( first_k&1 ? -1.0 : 1.0 );
                    // M_j^k
                    const int index_j = harmonic.getPreExpRedirJ(j);
                    // Inner_{j-n}^{k-l} : here points on the M2M transfer function/expansion term of degree n-j and order |k-l|
                    const int index_j_n = harmonic.getPreExpRedirJ(j-n);

                    int k = first_k;
                    for(; k <= (j-n+l) && k < l ; ++k){ /* l>=0 && k-l<0 */
                        const FComplexe I_j_n__l_k = M2M_Inner_transfer[index_j_n + l - k];

                        multipole_exp_target[index_j + k].incReal( pow_of_minus_1_for_k * pow_of_minus_1_for_l *
                                                    ((M_n__n_l.getReal() * I_j_n__l_k.getReal()) +
                                                     (M_n__n_l.getImag() * I_j_n__l_k.getImag())));
                        multipole_exp_target[index_j + k].incImag(pow_of_minus_1_for_k * pow_of_minus_1_for_l *
                                                   ((M_n__n_l.getImag() * I_j_n__l_k.getReal()) -
                                                    (M_n__n_l.getReal() * I_j_n__l_k.getImag())));

                        pow_of_minus_1_for_k = -pow_of_minus_1_for_k;
                    } // for k

                    for(/* k = l */; k <= (j - n + l) ; ++k){ // l>=0 && k-l>=0
                        const FComplexe I_j_n__k_l = M2M_Inner_transfer[index_j_n + k - l];

                        multipole_exp_target[index_j + k].incReal(
                                (M_n__n_l.getReal() * I_j_n__k_l.getReal()) -
                                (M_n__n_l.getImag() * I_j_n__k_l.getImag()));
                        multipole_exp_target[index_j + k].incImag(
                                (M_n__n_l.getImag() * I_j_n__k_l.getReal()) +
                                (M_n__n_l.getReal() * I_j_n__k_l.getImag()));

                    } // for k
                } // for j

                pow_of_minus_1_for_l = -pow_of_minus_1_for_l;
            } // for l
        } // for n
    }


    /** L2L
      *We compute the shift of local_exp_src from *p_center_of_exp_src to
      *p_center_of_exp_target, and set the result to local_exp_target.
      *
      *O_n^l (with n=0..P, l=-n..n) being the former local expansion terms
      *(whose center is *p_center_of_exp_src) we have for the new local
      *expansion terms (whose center is *p_center_of_exp_target):
      *
      *L_j^k = sum{n=j..P}
      *sum{l=-n..n}
      *O_n^l Inner_{n-j}^{l-k}(rho, alpha, beta)
      *
      *where (rho, alpha, beta) are the spherical coordinates of the vector :
      *p_center_of_exp_target - *p_center_of_exp_src
      *
      *Warning: if |l-k| > n-j, we do nothing.
      */
    void localToLocal(FComplexe* const FRestrict local_exp_target, const FComplexe* const FRestrict local_exp_src,
                      const FComplexe* const FRestrict L2L_tranfer){
        // L_j^k
        int index_j_k = 0;

        for (int j = 0 ; j <= devP ; ++j ){
            // (-1)^k
            FReal pow_of_minus_1_for_k = 1.0;

            for (int k = 0 ; k <= j ; ++k, ++index_j_k ){
                FComplexe L_j_k = local_exp_target[index_j_k];

                for (int n=j; n <= devP;++n){
                    // O_n^l : here points on the source multipole expansion term of degree n and order |l|
                    const int index_n = harmonic.getPreExpRedirJ(n);

                    int l = n - j + k;
                    // Inner_{n-j}^{l-k} : here points on the L2L transfer function/expansion term of degree n-j and order |l-k|
                    const int index_n_j = harmonic.getPreExpRedirJ(n-j);

                    for(/*l = n - j + k*/ ; l-k > 0 ;  --l){ /* l>0 && l-k>0 */
                        const FComplexe L_j_l = local_exp_src[index_n + l];
                        const FComplexe I_l_j__l_k = L2L_tranfer[index_n_j  + l - k];

                        L_j_k.incReal( (L_j_l.getReal() * I_l_j__l_k.getReal()) -
                                                    (L_j_l.getImag() * I_l_j__l_k.getImag()));
                        L_j_k.incImag( (L_j_l.getImag() * I_l_j__l_k.getReal()) +
                                                    (L_j_l.getReal() * I_l_j__l_k.getImag()));

                    }

                    // (-1)^l
                    FReal pow_of_minus_1_for_l = ((l&1) ? FReal(-1.0) : FReal(1.0));
                    for(/*l = k*/; l>0 && l>=j-n+k; --l){ /* l>0 && l-k<=0 */
                        const FComplexe L_j_l = local_exp_src[index_n + l];
                        const FComplexe I_l_j__l_k = L2L_tranfer[index_n_j  - l + k];

                        L_j_k.incReal( pow_of_minus_1_for_l * pow_of_minus_1_for_k *
                                                    ((L_j_l.getReal() * I_l_j__l_k.getReal()) +
                                                     (L_j_l.getImag() * I_l_j__l_k.getImag())));
                        L_j_k.incImag( pow_of_minus_1_for_l * pow_of_minus_1_for_k *
                                                    ((L_j_l.getImag() * I_l_j__l_k.getReal()) -
                                                     (L_j_l.getReal() * I_l_j__l_k.getImag())));

                        pow_of_minus_1_for_l = -pow_of_minus_1_for_l;
                     }

                    // l<=0 && l-k<=0
                    for(/*l = 0 ou l = j-n+k-1*/; l>=j-n+k; --l){
                        const FComplexe L_j_l = local_exp_src[index_n - l];
                        const FComplexe I_l_j__l_k = L2L_tranfer[index_n_j  - l + k];

                        L_j_k.incReal( pow_of_minus_1_for_k *
                                                    ((L_j_l.getReal() * I_l_j__l_k.getReal()) -
                                                     (L_j_l.getImag() * I_l_j__l_k.getImag())));
                        L_j_k.decImag( pow_of_minus_1_for_k *
                                                    ((L_j_l.getImag() * I_l_j__l_k.getReal()) +
                                                     (L_j_l.getReal() * I_l_j__l_k.getImag())));


                    }
                }//n

                local_exp_target[index_j_k] = L_j_k;

                pow_of_minus_1_for_k = -pow_of_minus_1_for_k;
            }//k
        }//j
    }

    /** L2P
      */
    void localToParticle(const FPoint& local_position,const FComplexe*const local_exp,
                         const FPoint& particlePosition,
                         const FReal physicalValue, FReal*const potential,
                         FReal*const forcesX,FReal*const forcesY,FReal*const forcesZ){
        //--------------- Forces ----------------//

        FReal force_vector_in_local_base_x = 0.0;
        FReal force_vector_in_local_base_y = 0.0;
        FReal force_vector_in_local_base_z = 0.0;

        const FSpherical spherical(particlePosition - local_position);
        harmonic.computeInnerTheta( spherical );

        int index_j_k = 1;

        for (int j = 1 ; j <= devP ; ++j ){
            {
                // k=0:
                // F_r:
                const FReal exp_term_aux_real = ( (harmonic.result(index_j_k).getReal() * local_exp[index_j_k].getReal()) - (harmonic.result(index_j_k).getImag() * local_exp[index_j_k].getImag()) );
                //const FReal exp_term_aux_imag = ( (harmonic.result(index_j_k).getReal() * local_exp[index_j_k].getImag()) + harmonic.result(index_j_k).getImag() * local_exp[index_j_k].getReal()) );
                force_vector_in_local_base_x = ( force_vector_in_local_base_x  + FReal(j) * exp_term_aux_real );
            }
            {
                // F_phi: k=0 => nothing to do for F_phi
                // F_theta:
                const FReal exp_term_aux_real = ( (harmonic.resultThetaDerivated(index_j_k).getReal() * local_exp[index_j_k].getReal()) - (harmonic.resultThetaDerivated(index_j_k).getImag() * local_exp[index_j_k].getImag()) );
                //const FReal exp_term_aux_imag = ( (harmonic.resultThetaDerivated(index_j_k).getReal() * local_exp[index_j_k].getImag()) + (harmonic.resultThetaDerivated(index_j_k).getImag() * local_exp[index_j_k].getReal()) );
                force_vector_in_local_base_y = ( force_vector_in_local_base_y + exp_term_aux_real );
            }

            ++index_j_k;

            // k>0:
            for (int k=1; k<=j ;++k, ++index_j_k){
                {
                    // F_r:
                    const FReal exp_term_aux_real = ( (harmonic.result(index_j_k).getReal() * local_exp[index_j_k].getReal()) - (harmonic.result(index_j_k).getImag() * local_exp[index_j_k].getImag()) );
                    const FReal exp_term_aux_imag = ( (harmonic.result(index_j_k).getReal() * local_exp[index_j_k].getImag()) + (harmonic.result(index_j_k).getImag() * local_exp[index_j_k].getReal()) );
                    force_vector_in_local_base_x = (force_vector_in_local_base_x  + FReal(2 * j) * exp_term_aux_real );
                    // F_phi:
                    force_vector_in_local_base_z = ( force_vector_in_local_base_z - FReal(2 * k) * exp_term_aux_imag);
                }
                {
                    // F_theta:
                    const FReal exp_term_aux_real = ( (harmonic.resultThetaDerivated(index_j_k).getReal() * local_exp[index_j_k].getReal()) - (harmonic.resultThetaDerivated(index_j_k).getImag() * local_exp[index_j_k].getImag()) );
                    //const FReal exp_term_aux_imag = ( (harmonic.resultThetaDerivated(index_j_k).getReal() * local_exp[index_j_k].getImag()) + (harmonic.resultThetaDerivated(index_j_k).getImag() * local_exp[index_j_k].getReal()) );
                    force_vector_in_local_base_y = (force_vector_in_local_base_y + FReal(2.0) * exp_term_aux_real );
                }

            }

        }
        // We want: - gradient(POTENTIAL_SIGN potential).
        // The -(- 1.0) computing is not the most efficient programming ...
        const FReal signe = 1.0;
        //if( FMath::Epsilon < spherical.getR()){
        // POURQUOI diviser par spherical.getR()
            force_vector_in_local_base_x = ( force_vector_in_local_base_x  * signe / spherical.getR()); // ????
        	// force_vector_in_local_base_x = ( force_vector_in_local_base_x ); // Better
            //
            force_vector_in_local_base_y = ( force_vector_in_local_base_y * signe / spherical.getR());
            force_vector_in_local_base_z = ( force_vector_in_local_base_z * signe / (spherical.getR() * spherical.getSinTheta()));
        //}
        /////////////////////////////////////////////////////////////////////
        //	???????????????????
        //spherical_position_Set_ph
        FReal ph = FMath::Fmod(spherical.getPhi(), FReal(2)*FMath::FPi);
        if (ph > M_PI) ph -= FReal(2) * FMath::FPi;
        if (ph < -M_PI + FMath::Epsilon)  ph += FReal(2) * FMath::FPi;

        //spherical_position_Set_th
        FReal th = FMath::Fmod(spherical.getTheta(), FReal(2) * FMath::FPi);
        if (th < 0.0) th += 2*FMath::FPi;
        if (th > FMath::FPi){
            th = 2*FMath::FPi - th;
            //spherical_position_Set_ph(p, spherical_position_Get_ph(p) + M_PI);
            ph = FMath::Fmod(ph + FMath::FPi, 2*FMath::FPi);
            if (ph > M_PI) ph -= 2*FMath::FPi;
            if (ph < -M_PI + FMath::Epsilon)  ph += 2 * FMath::FPi;
        }
        //spherical_position_Set_r
        //FReal rh = spherical.r;
        FAssertLF(spherical.getR() >= 0 , "R should be < 0!");

        const FReal cos_theta   = FMath::Cos(th);
        const FReal cos_phi     = FMath::Cos(ph);
        const FReal sin_theta   = FMath::Sin(th);
        const FReal sin_phi     = FMath::Sin(ph);
//
        //
        // Formulae below are OK
        //
        FReal force_vector_tmp_x = (
                cos_phi * sin_theta * force_vector_in_local_base_x  +
                cos_phi * cos_theta * force_vector_in_local_base_y +
                (-sin_phi) * force_vector_in_local_base_z);

        FReal force_vector_tmp_y = (
                sin_phi * sin_theta * force_vector_in_local_base_x  +
                sin_phi * cos_theta * force_vector_in_local_base_y +
                cos_phi * force_vector_in_local_base_z);

        FReal force_vector_tmp_z = (
                cos_theta * force_vector_in_local_base_x +
                (-sin_theta) * force_vector_in_local_base_y);

        force_vector_tmp_x *= physicalValue;
        force_vector_tmp_y *= physicalValue;
        force_vector_tmp_z *= physicalValue;

        (*forcesX) += force_vector_tmp_x;
        (*forcesY) += force_vector_tmp_y;
        (*forcesZ) += force_vector_tmp_z;
        //
        //--------------- Potential ----------------//
        //
        FReal result = 0.0;
        index_j_k    = 0;

        for(int j = 0 ; j<= devP ; ++j){
            // k=0
            harmonic.result(index_j_k) *= local_exp[index_j_k];
            result += harmonic.result(index_j_k).getReal();

            ++index_j_k;            

            // k>0
            for (int k = 1 ; k <= j ; ++k, ++index_j_k){
                harmonic.result(index_j_k) *= local_exp[index_j_k];
                result += 2 * harmonic.result(index_j_k).getReal();
            }
        }

        (*potential) += (result /* * physicalValue*/);
    }

public:
    /** Update a velocity of a particle
      *
      */
    void computeVelocity(ContainerClass*const FRestrict inParticles, const FReal DT){
        const FReal*const physicalValues = inParticles->getPhysicalValues();
        FReal*const forcesX = inParticles->getForcesX();
        FReal*const forcesY = inParticles->getForcesY();
        FReal*const forcesZ = inParticles->getForcesZ();
        FVector<FPoint>& velocities = inParticles->getVelocities();

        for(int idxPart = 0 ; idxPart < inParticles->getNbParticles() ; ++idxPart){
            const FReal physicalValue = physicalValues[idxPart];
            // Coef = 1/m * time/2
            const FReal coef = (FReal(1.0)/physicalValue) * (DT/FReal(2.0));

            // velocity = velocity + forces * coef
            FPoint forces_coef(forcesX[idxPart], forcesY[idxPart], forcesZ[idxPart]);
            forces_coef *= coef;
            velocities[idxPart] += (forces_coef);
        }
    }

    /** Update a position of a particle
      *
      */
    void updatePosition(ContainerClass*const FRestrict inParticles, const FReal DT){
        FReal*const positionsX = inParticles->getWPositions()[0];
        FReal*const positionsY = inParticles->getWPositions()[1];
        FReal*const positionsZ = inParticles->getWPositions()[2];
        FVector<FPoint>& velocities = inParticles->getVelocities();

        for(int idxPart = 0 ; idxPart < inParticles->getNbParticles() ; ++idxPart){
            FPoint velocity_dt( velocities[idxPart] );
            velocity_dt *= DT;
            positionsX[idxPart] += velocity_dt.getX();
            positionsY[idxPart] += velocity_dt.getY();
            positionsZ[idxPart] += velocity_dt.getZ();
        }
    }
};


#endif //FABSTRACTSPHERICALKERNEL_HPP


