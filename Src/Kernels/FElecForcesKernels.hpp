#ifndef FELECBASICKERNELS_HPP
#define FELECBASICKERNELS_HPP
// [--License--]

#include "../Components/FAbstractKernels.hpp"

#include "../Utils/FGlobal.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FMemUtils.hpp"

#include "../Containers/FTreeCoordinate.hpp"

#include "FHarmonic.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @brief
* Please read the license
*
*/
template< class ParticleClass, class CellClass, class ContainerClass>
class FElecForcesKernels : public FAbstractKernels<ParticleClass,CellClass,ContainerClass> {
    int devP;           //< The P
    int devM2lP;        //< A secondary P
    int treeHeight;     //< The height of the tree
    FHarmonic harmonic; //< The harmonic computation class

    FComplexe* preL2LTransitions; //< The pre-computation for the L2L based on the level
    FComplexe* preM2MTransitions; //< The pre-computation for the M2M based on the level

    FComplexe* preM2LTransitions; //< The pre-computation for the M2L based on the level and the 189 possibilities

    /** To access te preL2L/preM2M right vector */
    int indexTransition(const int level, const int child){
        return level * 8 * harmonic.getExpSize() + child * harmonic.getExpSize();
    }

    /** To access te preM2L right vector */
    int indexM2LTransition(const int idxLevel,const int idxX,const int idxY,const int idxZ){
        return (idxLevel * 7 * 7 * 7 * devM2lP) + (((((idxX+3) * 7) + (idxY+3)) * 7 ) + (idxZ+3)) * devM2lP;
    }

    /** Alloc and init pre-vectors*/
    void allocAndInit(const FReal boxWidth){
        preL2LTransitions = new FComplexe[treeHeight * 8 * harmonic.getExpSize()];
        preM2MTransitions = new FComplexe[treeHeight * 8 * harmonic.getExpSize()];

        FReal treeWidthAtLevel = boxWidth/2;
        for(int idxLevel = 0 ; idxLevel < treeHeight - 1 ; ++idxLevel ){
            const F3DPosition father(treeWidthAtLevel,treeWidthAtLevel,treeWidthAtLevel);
            treeWidthAtLevel /= 2;

            for(int idxChild = 0 ; idxChild < 8 ; ++idxChild ){
                FTreeCoordinate childBox;
                childBox.setPositionFromMorton(idxChild,1);

                const F3DPosition M2MVector (
                        father.getX() - (treeWidthAtLevel * FReal(1 + (childBox.getX() * 2))),
                        father.getY() - (treeWidthAtLevel * FReal(1 + (childBox.getY() * 2))),
                        father.getZ() - (treeWidthAtLevel * FReal(1 + (childBox.getZ() * 2)))
                        );

                harmonic.computeInner(FSpherical(M2MVector));
                FMemUtils::copyall<FComplexe>(&preM2MTransitions[indexTransition(idxLevel,idxChild)], harmonic.result(), harmonic.getExpSize());

                const F3DPosition L2LVector (
                        (treeWidthAtLevel * FReal(1 + (childBox.getX() * 2))) - father.getX(),
                        (treeWidthAtLevel * FReal(1 + (childBox.getY() * 2))) - father.getY(),
                        (treeWidthAtLevel * FReal(1 + (childBox.getZ() * 2))) - father.getZ()
                        );

                harmonic.computeInner(FSpherical(L2LVector));
                FMemUtils::copyall<FComplexe>(&preL2LTransitions[indexTransition(idxLevel,idxChild)], harmonic.result(), harmonic.getExpSize());
           }
        }

        // M2L transfer, there is a maximum of 3 neighbors in each direction,
        // so 6 in each dimension
        treeWidthAtLevel = boxWidth/2;
        preM2LTransitions = new FComplexe[treeHeight * (7 * 7 * 7) * devM2lP];
        for(int idxLevel = 0 ; idxLevel < treeHeight ; ++idxLevel ){
            for(int idxX = -3 ; idxX <= 3 ; ++idxX ){
                for(int idxY = -3 ; idxY <= 3 ; ++idxY ){
                    for(int idxZ = -3 ; idxZ <= 3 ; ++idxZ ){
                        if(idxX || idxY || idxZ){
                            const F3DPosition relativePos( FReal(idxX) * treeWidthAtLevel , FReal(idxY) * treeWidthAtLevel , FReal(idxZ) * treeWidthAtLevel );
                            harmonic.computeOuter(FSpherical(relativePos));
                            FMemUtils::copyall<FComplexe>(&preM2LTransitions[indexM2LTransition(idxLevel,idxX,idxY,idxZ)], harmonic.result(), harmonic.getExpSize());
                        }
                    }
                }
            }
            treeWidthAtLevel /= 2;
        }
    }


public:
    /** Kernel constructor */
    FElecForcesKernels(const int inDevP, const int inTreeHeight, const FReal boxWidth)
        : devP(inDevP), devM2lP(int(((inDevP*2)+1) * ((inDevP*2)+2) * 0.5)), treeHeight(inTreeHeight), harmonic(inDevP),
          preL2LTransitions(0), preM2MTransitions(0), preM2LTransitions(0) {

        allocAndInit(boxWidth);
    }

    /** Default destructor */
    ~FElecForcesKernels(){
        delete[] preL2LTransitions;
        delete[] preM2MTransitions;
        delete[] preM2LTransitions;
    }

    /** P2M with a cell and all its particles */
    void P2M(CellClass* const inPole, const ContainerClass* const inParticles) {
        FComplexe* const cellMultiPole = inPole->getMultipole();
        // Copying the position is faster than using cell position
        const F3DPosition polePosition = inPole->getPosition();
        // For all particles in the leaf box
        typename ContainerClass::ConstBasicIterator iterParticle(*inParticles);
        while( iterParticle.hasNotFinished()){
            // P2M
            particleToMultiPole(cellMultiPole, polePosition, iterParticle.data());
            iterParticle.gotoNext();
        }
    }

    /** M2M with a cell and all its child */
    void M2M(CellClass* const FRestrict inPole, const CellClass *const FRestrict *const FRestrict inChild, const int inLevel) {
        FComplexe* const multipole_exp_target = inPole->getMultipole();
        // iter on each child and process M2M
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            if(inChild[idxChild]){
                multipoleToMultipole(multipole_exp_target, inChild[idxChild]->getMultipole(), &preM2MTransitions[indexTransition(inLevel,idxChild)]);
            }
        }
    }

    /** M2L with a cell and all the existing neighbors */
    void M2L(CellClass* const FRestrict pole, const CellClass* distantNeighbors[189],
             const int size, const int inLevel) {
        const FTreeCoordinate& coordCenter = pole->getCoordinate();
        // For all neighbors compute M2L
        for(int idxSize = 0 ; idxSize < size ; ++idxSize){
            const FTreeCoordinate& coordNeighbors = distantNeighbors[idxSize]->getCoordinate();
            const FComplexe* const transitionVector = &preM2LTransitions[indexM2LTransition(inLevel,
                                                  (coordNeighbors.getX() - coordCenter.getX()),
                                                  (coordNeighbors.getY() - coordCenter.getY()),
                                                  (coordNeighbors.getZ() - coordCenter.getZ()))];

            multipoleToLocal(pole->getLocal(), distantNeighbors[idxSize]->getMultipole(), transitionVector);
        }
    }

    /** L2L with a cell and all its child */
    void L2L(const CellClass* const FRestrict pole, CellClass* FRestrict *const FRestrict child, const int inLevel) {
        // iter on each child and process L2L
        for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
            if(child[idxChild]){
                localToLocal(child[idxChild]->getLocal(), pole->getLocal(), &preL2LTransitions[indexTransition(inLevel,idxChild)]);
            }
        }
    }

    /** L2P with a cell and all its particles */
    void L2P(const CellClass* const local, ContainerClass* const particles){
        const FComplexe* const cellLocal = local->getLocal();
        // Copying the position is faster than using cell position
        const F3DPosition localPosition = local->getPosition();
        // For all particles in the leaf box
        typename ContainerClass::BasicIterator iterTarget(*particles);
        while( iterTarget.hasNotFinished() ){
            // L2P
            localToParticle(&iterTarget.data(), localPosition, cellLocal);
            iterTarget.gotoNext();
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
    void P2P(ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict sources,
             const ContainerClass* const directNeighbors[26], const int size) {

        { // Compute interaction in this leaf
            typename ContainerClass::BasicIterator iterTarget(*targets);
            while( iterTarget.hasNotFinished() ){
                // We copy the target particle to work with a particle in the heap
                ParticleClass target( iterTarget.data() );

                // For all the source particles in the same leaf
                typename ContainerClass::ConstBasicIterator iterSameBox(*sources);
                while( iterSameBox.hasNotFinished() ){
                    //(&iterSameBox.data() != &iterTarget.data())
                    directInteraction(&target, iterSameBox.data());
                    iterSameBox.gotoNext();
                }
                // Set data and progress
                iterTarget.setData(target);
                iterTarget.gotoNext();
            }
        }
        { // Compute interactions with other leaves
            // For all the neigbors leaves
            for(int idxDirectNeighbors = 0 ; idxDirectNeighbors < size ; ++idxDirectNeighbors){
                // For all particles in current leaf
                typename ContainerClass::BasicIterator iterTarget(*targets);
                while( iterTarget.hasNotFinished() ){
                    ParticleClass target( iterTarget.data() );
                    // For all the particles in the other leaf
                    typename ContainerClass::ConstBasicIterator iterSource(*directNeighbors[idxDirectNeighbors]);
                    while( iterSource.hasNotFinished() ){
                        directInteraction(&target, iterSource.data());
                        iterSource.gotoNext();
                    }
                    // Set data and progress
                    iterTarget.setData(target);
                    iterTarget.gotoNext();
                }
            }
        }
    }

    /** This P2P has to be used when target == sources
      * It will proceed a direct interation >> mutual
      *
      * It takes all the particles from the current leaf,
      * then it computes the interactions in this leaf,
      * then it computes the  inteactions between this leaf and the
      * neighbors.
      */
    void P2P(const MortonIndex inCurrentIndex,
             ContainerClass* const FRestrict targets, const ContainerClass* const FRestrict /*sources*/,
             ContainerClass* const directNeighbors[26], const MortonIndex inNeighborsIndex[26], const int size){
        { // Compute interaction in this leaf
            typename ContainerClass::BasicIterator iterTarget(*targets);
            while( iterTarget.hasNotFinished() ){
                // We copy the target particle to work with a particle in the heap
                ParticleClass target( iterTarget.data() );

                // For all particles after the current one
                typename ContainerClass::BasicIterator iterSameBox = iterTarget;
                iterSameBox.gotoNext();
                while( iterSameBox.hasNotFinished() ){
                    directInteractionMutual(&target, &iterSameBox.data());
                    iterSameBox.gotoNext();
                }
                // Set data and progress
                iterTarget.setData(target);
                iterTarget.gotoNext();
            }
        }
        { // Compute interactions with other leaves
            // For all the neigbors leaves
            for(int idxDirectNeighbors = 0 ; idxDirectNeighbors < size ; ++idxDirectNeighbors){
                if(inCurrentIndex < inNeighborsIndex[idxDirectNeighbors] ){
                    // For all particles in current leaf
                    typename ContainerClass::BasicIterator iterTarget(*targets);
                    while( iterTarget.hasNotFinished() ){
                        ParticleClass target( iterTarget.data() );
                        // For all the particles in the other leaf
                        typename ContainerClass::BasicIterator iterSource(*directNeighbors[idxDirectNeighbors]);
                        while( iterSource.hasNotFinished() ){
                            directInteractionMutual(&target, &iterSource.data());
                            iterSource.gotoNext();
                        }
                        // Set data and progress
                        iterTarget.setData(target);
                        iterTarget.gotoNext();
                    }
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////
    //                                  Periodic
    ///////////////////////////////////////////////////////////////////////////////

    /** Before Downward */
    void M2L(CellClass* const FRestrict , const CellClass* [189], FTreeCoordinate [189], const int , const int ) {
    }


    /** After Downward */
    void P2P(const MortonIndex ,
             ContainerClass* const FRestrict , const ContainerClass* const FRestrict ,
             ContainerClass* const [26], const FTreeCoordinate [26], const int ) {
    }

    ///////////////////////////////////////////////////////////////////////////////
    //                                  Computation
    ///////////////////////////////////////////////////////////////////////////////
private:

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
    void particleToMultiPole(FComplexe*const cellMultiPole, const F3DPosition& inPolePosition ,
                             const ParticleClass& particle){
        harmonic.computeInner( FSpherical(particle.getPosition() - inPolePosition) );

        FReal minus_one_pow_j = 1.0;//(-1)^j
        const FReal valueParticle = particle.getPhysicalValue();
        int p_exp_term = 0; // p_Y_term

        for(int jIdx = 0 ; jIdx <= devP ; ++jIdx){
            for(int kIdx = 0 ; kIdx <= jIdx ; ++kIdx, ++p_exp_term){
                harmonic.result(p_exp_term).mulRealAndImag( valueParticle * minus_one_pow_j );
                cellMultiPole[p_exp_term] += harmonic.result(p_exp_term);
            }

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
                              const FComplexe* const FRestrict M2M_transfer){

        for(int n = 0 ; n <= FMB_Info_P ; ++n ){
            // l<0 // (-1)^l
            FReal pow_of_minus_1_for_l = static_cast<FReal>( n % 2 ? -1.0 : 1.0);

            // O_n^l : here points on the source multipole expansion term of degree n and order |l|
            const FComplexe* p_src_exp_term = multipole_exp_src + harmonic.getPreExpRedirJ(n) + n;

            int l = -n;
            for(; l<0 ; ++l, --p_src_exp_term, pow_of_minus_1_for_l = -pow_of_minus_1_for_l){

                for(int j = n ; j<= FMB_Info_P ; ++j ){
                    // M_j^k
                    FComplexe *p_target_exp_term = multipole_exp_target + harmonic.getPreExpRedirJ(j);
                    // Inner_{j-n}^{k-l} : here points on the M2M transfer function/expansion term of degree n-j and order |k-l|
                    const FComplexe *p_Inner_term= M2M_transfer + harmonic.getPreExpRedirJ(j-n)-l /* k==0 */;

                    // since n-j+l<0
                    for(int k=0 ; k <= (j-n+l) ; ++k, ++p_target_exp_term, ++p_Inner_term){ // l<0 && k>=0 => k-l>0
                        p_target_exp_term->incReal( pow_of_minus_1_for_l *
                                                    ((p_src_exp_term->getReal() * p_Inner_term->getReal()) +
                                                     (p_src_exp_term->getImag() * p_Inner_term->getImag())));
                        p_target_exp_term->incImag( pow_of_minus_1_for_l *
                                                    ((p_src_exp_term->getReal() * p_Inner_term->getImag()) -
                                                     (p_src_exp_term->getImag() * p_Inner_term->getReal())));

                     } // for k
                } // for j
            } // for l

            // l>=0
            for(; l <= n ; ++l, ++p_src_exp_term, pow_of_minus_1_for_l = -pow_of_minus_1_for_l){

                for( int j=n ; j <= FMB_Info_P ; ++j ){
                    // (-1)^k
                    FReal pow_of_minus_1_for_k = static_cast<FReal>( FMath::Max(0,n-j+l) %2 ? -1.0 : 1.0 );
                    // M_j^k
                    FComplexe *p_target_exp_term = multipole_exp_target + harmonic.getPreExpRedirJ(j) + FMath::Max(0,n-j+l);
                    // Inner_{j-n}^{k-l} : here points on the M2M transfer function/expansion term of degree n-j and order |k-l|
                    const FComplexe *p_Inner_term = M2M_transfer + harmonic.getPreExpRedirJ(j-n) + l - FMath::Max(0,n-j+l);// -(k-l)

                    int k = FMath::Max(0,n-j+l);
                    for(; k <= (j-n+l) && (k-l) < 0 ; ++k, ++p_target_exp_term, --p_Inner_term, pow_of_minus_1_for_k = -pow_of_minus_1_for_k){ /* l>=0 && k-l<0 */
                        p_target_exp_term->incReal( pow_of_minus_1_for_k * pow_of_minus_1_for_l *
                                                    ((p_src_exp_term->getReal() * p_Inner_term->getReal()) +
                                                     (p_src_exp_term->getImag() * p_Inner_term->getImag())));
                        p_target_exp_term->incImag(pow_of_minus_1_for_k * pow_of_minus_1_for_l *
                                                   ((p_src_exp_term->getImag() * p_Inner_term->getReal()) -
                                                    (p_src_exp_term->getReal() * p_Inner_term->getImag())));
                    } // for k

                    for(; k <= (j - n + l) ; ++k, ++p_target_exp_term, ++p_Inner_term){ // l>=0 && k-l>=0
                        p_target_exp_term->incReal(
                                (p_src_exp_term->getReal() * p_Inner_term->getReal()) -
                                (p_src_exp_term->getImag() * p_Inner_term->getImag()));
                        p_target_exp_term->incImag(
                                (p_src_exp_term->getImag() * p_Inner_term->getReal()) +
                                (p_src_exp_term->getReal() * p_Inner_term->getImag()));
                    } // for k
                } // for j
            } // for l
        } // for n
    }


    /** M2L
    *We compute the conversion of multipole_exp_src in *p_center_of_exp_src to
    *a local expansion in *p_center_of_exp_target, and add the result to local_exp_target.
    *
    *O_n^l (with n=0..P, l=-n..n) being the former multipole expansion terms
    *(whose center is *p_center_of_multipole_exp_src) we have for the new local
    *expansion terms (whose center is *p_center_of_local_exp_target):
    *
    *L_j^k = sum{n=0..+}
    *sum{l=-n..n}
    *O_n^l Outer_{j+n}^{-k-l}(rho, alpha, beta)
    *
    *where (rho, alpha, beta) are the spherical coordinates of the vector :
    *p_center_of_local_exp_src - *p_center_of_multipole_exp_target
    *
    *Remark: here we have always j+n >= |-k-l|
      *
      */
    void multipoleToLocal(FComplexe*const FRestrict local_exp, const FComplexe* const FRestrict multipole_exp_src,
                          const FComplexe* const FRestrict M2L_transfer){
        FComplexe* p_target_exp_term = local_exp;

        // L_j^k
        int start_for_j = 0;

        //    HPMSTART(51, "M2L computation (loops)");
        for (int j = start_for_j ; j <= FMB_Info_P ; ++j){

            int stop_for_n = devP;
            //stop_for_n = FMB_Info_P - j;

            // (-1)^k
            FReal pow_of_minus_1_for_k = 1.0;
            for (int k = 0 ; k <= j ; ++k, pow_of_minus_1_for_k = -pow_of_minus_1_for_k, ++p_target_exp_term){

                // (-1)^n
                FReal pow_of_minus_1_for_n = 1.0;
                for (int n = 0 ; n <= stop_for_n ; ++n, pow_of_minus_1_for_n = -pow_of_minus_1_for_n){

                    // O_n^l : here points on the source multipole expansion term of degree n and order |l|
                    const FComplexe *p_src_exp_term = multipole_exp_src + harmonic.getPreExpRedirJ(n) + n;
                    // Outer_{j+n}^{-k-l} : here points on the M2L transfer function/expansion term of degree j+n and order |-k-l|
                    const FComplexe *p_Outer_term = M2L_transfer + harmonic.getPreExpRedirJ(n+j) + k+n;
                    FReal pow_of_minus_1_for_l = pow_of_minus_1_for_n; // (-1)^l
                    // We start with l=n (and not l=-n) so that we always set p_Outer_term to a correct value in the first loop.
                    int l=n;
                    for ( ; l>0 ; --l, pow_of_minus_1_for_l = -pow_of_minus_1_for_l, --p_src_exp_term, --p_Outer_term){ // we have -k-l<0 and l>0
                        p_target_exp_term->incReal( pow_of_minus_1_for_l * pow_of_minus_1_for_k *
                                                    ((p_src_exp_term->getReal() * p_Outer_term->getReal()) +
                                                     (p_src_exp_term->getImag() * p_Outer_term->getImag())));
                        p_target_exp_term->incImag( pow_of_minus_1_for_l * pow_of_minus_1_for_k *
                                                    ((p_src_exp_term->getImag() * p_Outer_term->getReal()) -
                                                     (p_src_exp_term->getReal() * p_Outer_term->getImag())));
                    }

                    for (; l>=-n && -k-l<0 ; --l, pow_of_minus_1_for_l = -pow_of_minus_1_for_l, ++p_src_exp_term, --p_Outer_term){ // we have -k-l<0 and l<=0
                        p_target_exp_term->incReal( pow_of_minus_1_for_k *
                                                    ((p_src_exp_term->getReal() * p_Outer_term->getReal()) -
                                                     (p_src_exp_term->getImag() * p_Outer_term->getImag())));
                        p_target_exp_term->decImag(  pow_of_minus_1_for_k *
                                                     ((p_src_exp_term->getImag() * p_Outer_term->getReal()) +
                                                      (p_src_exp_term->getReal() * p_Outer_term->getImag())));
                    }

                    for (; l>=-n; --l, pow_of_minus_1_for_l = -pow_of_minus_1_for_l, ++p_src_exp_term, ++p_Outer_term){ // we have -k-l>=0 and l<=0
                        p_target_exp_term->incReal( pow_of_minus_1_for_l *
                                                    ((p_src_exp_term->getReal() * p_Outer_term->getReal()) +
                                                     (p_src_exp_term->getImag() * p_Outer_term->getImag())));
                        p_target_exp_term->incImag( pow_of_minus_1_for_l *
                                                    ((p_src_exp_term->getReal() * p_Outer_term->getImag()) -
                                                     (p_src_exp_term->getImag() * p_Outer_term->getReal())));
                    }
                }
            }
        }
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
        FComplexe* p_target_exp_term = local_exp_target;
        for (int j=0 ; j<= FMB_Info_P ; ++j){
            // (-1)^k
            FReal pow_of_minus_1_for_k = 1.0;
            for (int k=0 ; k <= j ; ++k, pow_of_minus_1_for_k = -pow_of_minus_1_for_k, ++p_target_exp_term){
                for (int n=j; n<=FMB_Info_P;++n){
                    // O_n^l : here points on the source multipole expansion term of degree n and order |l|
                    const FComplexe* p_src_exp_term = local_exp_src + harmonic.getPreExpRedirJ(n) + n-j+k;

                    int l = n-j+k;
                    // Inner_{n-j}^{l-k} : here points on the L2L transfer function/expansion term of degree n-j and order |l-k|
                    const FComplexe* p_Inner_term = L2L_tranfer + harmonic.getPreExpRedirJ(n-j) + l-k;

                    for ( ; l-k>0;  --l, --p_src_exp_term, --p_Inner_term){ /* l>0 && l-k>0 */
                        p_target_exp_term->incReal( (p_src_exp_term->getReal() * p_Inner_term->getReal()) -
                                                    (p_src_exp_term->getImag() * p_Inner_term->getImag()));
                        p_target_exp_term->incImag( (p_src_exp_term->getImag() * p_Inner_term->getReal()) +
                                                    (p_src_exp_term->getReal() * p_Inner_term->getImag()));
                    }

                    // (-1)^l
                    FReal pow_of_minus_1_for_l = static_cast<FReal>((l%2) ? -1.0 : 1.0);
                    for (; l>0 && l>=j-n+k; --l, pow_of_minus_1_for_l = -pow_of_minus_1_for_l, --p_src_exp_term, ++p_Inner_term){ /* l>0 && l-k<=0 */
                        p_target_exp_term->incReal( pow_of_minus_1_for_l * pow_of_minus_1_for_k *
                                                    ((p_src_exp_term->getReal() * p_Inner_term->getReal()) +
                                                     (p_src_exp_term->getImag() * p_Inner_term->getImag())));
                        p_target_exp_term->incImag( pow_of_minus_1_for_l * pow_of_minus_1_for_k *
                                                    ((p_src_exp_term->getImag() * p_Inner_term->getReal()) -
                                                     (p_src_exp_term->getReal() * p_Inner_term->getImag())));
                     }

                    // l<=0 && l-k<=0
                    for (; l>=j-n+k; --l, ++p_src_exp_term, ++p_Inner_term){
                        p_target_exp_term->incReal( pow_of_minus_1_for_k *
                                                    ((p_src_exp_term->getReal() * p_Inner_term->getReal()) -
                                                     (p_src_exp_term->getImag() * p_Inner_term->getImag())));
                        p_target_exp_term->decImag( pow_of_minus_1_for_k *
                                                    ((p_src_exp_term->getImag() * p_Inner_term->getReal()) +
                                                     (p_src_exp_term->getReal() * p_Inner_term->getImag())));
                    }
                }
            }
        }
    }

    /** L2P
      */
    void localToParticle(ParticleClass*const particle, const F3DPosition& local_position,
                         const FComplexe*const local_exp){
        FReal force_vector_in_local_base_x = 0;
        FReal force_vector_in_local_base_y = 0;
        FReal force_vector_in_local_base_z = 0;

        const FSpherical spherical(particle->getPosition() - local_position);
        harmonic.computeInnerTheta( spherical );

        // The maximum degree used here will be P.
        const FComplexe* p_Y_term = harmonic.result() + 1;
        const FComplexe* p_Y_theta_derivated_term = harmonic.resultThetaDerivated() + 1;
        const FComplexe* p_local_exp_term = local_exp + 1;

        for (int j = 1 ; j <= devP ; ++j ){
            FReal exp_term_aux_real = 0.0;
            FReal exp_term_aux_imag = 0.0;

            // k=0:
            // F_r:
            exp_term_aux_real = ( (p_Y_term->getReal() * p_local_exp_term->getReal()) - (p_Y_term->getImag() * p_local_exp_term->getImag()) );
            exp_term_aux_imag = ( (p_Y_term->getReal() * p_local_exp_term->getImag()) + (p_Y_term->getImag() * p_local_exp_term->getReal()) );

            force_vector_in_local_base_x = ( force_vector_in_local_base_x  + FReal(j) * exp_term_aux_real );
            // F_phi: k=0 => nothing to do for F_phi
            // F_theta:
            exp_term_aux_real = ( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getReal()) - (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getImag()) );
            exp_term_aux_imag = ( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getImag()) + (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getReal()) );

            force_vector_in_local_base_y = ( force_vector_in_local_base_y + exp_term_aux_real );

            ++p_local_exp_term;
            ++p_Y_term;
            ++p_Y_theta_derivated_term;


            // k>0:
            for (int k=1; k<=j ;++k, ++p_local_exp_term, ++p_Y_term, ++p_Y_theta_derivated_term){
                // F_r:

                exp_term_aux_real = ( (p_Y_term->getReal() * p_local_exp_term->getReal()) - (p_Y_term->getImag() * p_local_exp_term->getImag()) );
                exp_term_aux_imag = ( (p_Y_term->getReal() * p_local_exp_term->getImag()) + (p_Y_term->getImag() * p_local_exp_term->getReal()) );

                force_vector_in_local_base_x = (force_vector_in_local_base_x  + FReal(2 * j) * exp_term_aux_real );
                // F_phi:
                force_vector_in_local_base_z = ( force_vector_in_local_base_z - FReal(2 * k) * exp_term_aux_imag);
                // F_theta:

                exp_term_aux_real = ( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getReal()) - (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getImag()) );
                exp_term_aux_imag = ( (p_Y_theta_derivated_term->getReal() * p_local_exp_term->getImag()) + (p_Y_theta_derivated_term->getImag() * p_local_exp_term->getReal()) );

                force_vector_in_local_base_y = (force_vector_in_local_base_y + FReal(2.0) * exp_term_aux_real );

           }

        }
        // We want: - gradient(POTENTIAL_SIGN potential).
        // The -(- 1.0) computing is not the most efficient programming ...
        force_vector_in_local_base_x = ( force_vector_in_local_base_x  * FReal(-1.0) / spherical.getR());
        force_vector_in_local_base_y = ( force_vector_in_local_base_y * FReal(-1.0) / spherical.getR());
        force_vector_in_local_base_z = ( force_vector_in_local_base_z * FReal(-1.0) / (spherical.getR() * spherical.getSinTheta()));

        /////////////////////////////////////////////////////////////////////

        //spherical_position_Set_ph
        //FMB_INLINE COORDINATES_T angle_Convert_in_MinusPi_Pi(COORDINATES_T a){
        FReal ph = FMath::Fmod(spherical.getPhi(), FReal(2)*FMath::FPi);
        if (ph > M_PI) ph -= FReal(2) * FMath::FPi;
        if (ph < -M_PI + FMath::Epsilon)  ph += FReal(2) * FMath::Epsilon;

        //spherical_position_Set_th
        FReal th = FMath::Fmod(FMath::ACos(spherical.getCosTheta()), FReal(2) * FMath::FPi);
        if (th < 0.0) th += 2*FMath::FPi;
        if (th > FMath::FPi){
            th = 2*FMath::FPi - th;
            //spherical_position_Set_ph(p, spherical_position_Get_ph(p) + M_PI);
            ph = FMath::Fmod(ph + FMath::FPi, 2*FMath::FPi);
            if (ph > M_PI) ph -= 2*FMath::FPi;
            if (ph < -M_PI + FMath::Epsilon)  ph += 2 * FMath::Epsilon;
            th = FMath::Fmod(th, 2*FMath::FPi);
            if (th > M_PI) th -= 2*FMath::FPi;
            if (th < -M_PI + FMath::Epsilon)  th += 2 * FMath::Epsilon;
        }
        //spherical_position_Set_r
        //FReal rh = spherical.r;
        if (spherical.getR() < 0){
            //rh = -spherical.r;
            //spherical_position_Set_ph(p, M_PI - spherical_position_Get_th(p));
            ph = FMath::Fmod(FMath::FPi - th, 2*FMath::FPi);
            if (ph > M_PI) ph -= 2*FMath::FPi;
            if (ph < -M_PI + FMath::Epsilon)  ph += 2 * FMath::Epsilon;
            //spherical_position_Set_th(p, spherical_position_Get_th(p) + M_PI);
            th = FMath::Fmod(th + FMath::FPi, 2*FMath::FPi);
            if (th < 0.0) th += 2*FMath::FPi;
            if (th > FMath::FPi){
                th = 2*FMath::FPi - th;
                //spherical_position_Set_ph(p, spherical_position_Get_ph(p) + M_PI);
                ph = FMath::Fmod(ph + FMath::FPi, 2*FMath::FPi);
                if (ph > M_PI) ph -= 2*FMath::FPi;
                if (ph < -M_PI + FMath::Epsilon)  ph += 2 * FMath::Epsilon;
                th = FMath::Fmod(th, 2*FMath::FPi);
                if (th > M_PI) th -= 2*FMath::FPi;
                if (th < -M_PI + FMath::Epsilon)  th += 2 * FMath::Epsilon;
            }
        }

        const FReal cos_theta   = FMath::Cos(th);
        const FReal cos_phi     = FMath::Cos(ph);
        const FReal sin_theta   = FMath::Sin(th);
        const FReal sin_phi     = FMath::Sin(ph);

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

        const FReal physicalValue = particle->getPhysicalValue();
        force_vector_tmp_x *= physicalValue;
        force_vector_tmp_y *= physicalValue;
        force_vector_tmp_z *= physicalValue;

        particle->incForces( force_vector_tmp_x, force_vector_tmp_y, force_vector_tmp_z );


        FReal result = 0.0;
        p_Y_term = harmonic.result() + 1;;
        const FComplexe* local_exp_iter = local_exp;
        for(int j = 0 ; j<= FMB_Info_P ; ++j){
            // k=0
            (*p_Y_term) *= (*local_exp_iter);
            result += p_Y_term->getReal();
            ++p_Y_term;
            ++local_exp_iter;

            // k>0
            for (int k=1; k<=j ;++k, ++p_Y_term, ++local_exp_iter){
                (*p_Y_term) *= (*local_exp_iter);
                result += 2 * p_Y_term->getReal();
            }
        }

        particle->incPotential(result);
    }


    /** P2P mutual interaction
      * F = q * q' / r²
      */
    void directInteractionMutual(ParticleClass*const FRestrict target, ParticleClass*const FRestrict source){

        FReal dx = -(target->getPosition().getX() - source->getPosition().getX());
        FReal dy = -(target->getPosition().getY() - source->getPosition().getY());
        FReal dz = -(target->getPosition().getZ() - source->getPosition().getZ());

        FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
        FReal inv_distance = FMath::Sqrt(inv_square_distance);

        inv_square_distance *= inv_distance;
        inv_square_distance *= target->getPhysicalValue() * source->getPhysicalValue();

        dx *= inv_square_distance;
        dy *= inv_square_distance;
        dz *= inv_square_distance;

        target->incForces( dx, dy, dz);
        target->incPotential( inv_distance  * source->getPhysicalValue() );

        source->incForces( (-dx), (-dy), (-dz));
        source->incPotential( inv_distance * target->getPhysicalValue() );
    }

    /** P2P NO mutual interaction
      * F = q * q' / r²
      */
    void directInteraction(ParticleClass*const FRestrict target, const ParticleClass& source){

        FReal dx = -(target->getPosition().getX() - source.getPosition().getX());
        FReal dy = -(target->getPosition().getY() - source.getPosition().getY());
        FReal dz = -(target->getPosition().getZ() - source.getPosition().getZ());

        FReal inv_square_distance = FReal(1.0) / (dx*dx + dy*dy + dz*dz);
        FReal inv_distance = FMath::Sqrt(inv_square_distance);

        inv_square_distance *= inv_distance;
        inv_square_distance *= target->getPhysicalValue() * source.getPhysicalValue();

        dx *= inv_square_distance;
        dy *= inv_square_distance;
        dz *= inv_square_distance;

        target->incForces( dx, dy, dz);
        target->incPotential( inv_distance  * source.getPhysicalValue() );
    }
};


#endif //FELECBASICKERNELS_HPP

// [--END--]
