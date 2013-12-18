// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================

// ==== CMAKE =====
// @FUSE_BLAS
// ================


#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FBlas.hpp"

#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssert.hpp"
#include "../../Src/Utils/FPoint.hpp"


#include "../../Src/Kernels/Uniform/FUnifInterpolator.hpp"
#include "../../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Uniform/FUnifTensor.hpp"

// Check DFT
#include "../../Src/Utils/FDft.hpp"
#include "../../Src/Utils/FComplexe.hpp"


#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"




/**
 * In this file we show how to use octree
 */

int main(int, char **){
  typedef FP2PParticleContainer ContainerClass;
  typedef FSimpleLeaf<ContainerClass> LeafClass;
  typedef FInterpMatrixKernelR MatrixKernelClass;


  ///////////////////////What we do/////////////////////////////
  std::cout << "\nTask: Compute interactions between source particles in leaf Y and target\n";
  std::cout << " particles in leaf X. Compare the fast summation K ~ Lx K Ly' with the\n";
  std::cout << " direct computation.\n" << std::endl;
  //////////////////////////////////////////////////////////////

  MatrixKernelClass MatrixKernel;
  const FReal FRandMax = FReal(RAND_MAX);
  FTic time;


  // Leaf size
  FReal width = FReal(3.723);

  ////////////////////////////////////////////////////////////////////
  LeafClass X;
  FPoint cx(0., 0., 0.);
  const unsigned long M = 20000;
  std::cout << "Fill the leaf X of width " << width
            << " centered at cx=" << cx << " with M=" << M << " target particles" << std::endl;
  {
    for(unsigned long i=0; i<M; ++i){
      FReal x = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getX();
      FReal y = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getY();
      FReal z = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getZ();
      X.push(FPoint(x, y, z), FReal(rand())/FRandMax);
    }
  }


  ////////////////////////////////////////////////////////////////////
  LeafClass Y;
  //  FPoint cy(FReal(2.)*width, 0., 0.);
  FPoint cy(FReal(2.)*width, FReal(2.)*width, 0.);

  const unsigned long N = 20000;
  std::cout << "Fill the leaf Y of width " << width
            << " centered at cy=" << cy	<< " with N=" << N << " target particles" << std::endl;
  {
    for(unsigned long i=0; i<N; ++i){
      FReal x = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getX();
      FReal y = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getY();
      FReal z = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getZ();
      Y.push(FPoint(x, y, z), FReal(rand())/FRandMax);
    }
  }



  ////////////////////////////////////////////////////////////////////
  // approximative computation
  const unsigned int ORDER = 3;
  const unsigned int nnodes = TensorTraits<ORDER>::nnodes;
  typedef FUnifInterpolator<ORDER> InterpolatorClass;
  InterpolatorClass S;

  std::cout << "\nCompute interactions approximatively, interpolation order = " << ORDER << " ..." << std::endl;

  std::cout << "\nP2M ... " << std::flush;
  time.tic();
  // Anterpolate: W_n = \sum_j^N S(y_j,\bar y_n) * w_j
  FReal W[nnodes]; // multipole expansion
  S.applyP2M(cy, width, W, Y.getSrc()); // the multipole expansions are set to 0 in S.applyP2M
  std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

  std::cout << "M2L ... " << std::flush;
  time.tic();
  // Multipole to local: F_m = \sum_n^L K(\bar x_m, \bar y_n) * W_n
  FPoint rootsX[nnodes], rootsY[nnodes];
  FUnifTensor<ORDER>::setRoots(cx, width, rootsX);
  FUnifTensor<ORDER>::setRoots(cy, width, rootsY);

  FReal F[nnodes]; // local expansion
  for (unsigned int i=0; i<nnodes; ++i) {
    F[i] = FReal(0.);
    for (unsigned int j=0; j<nnodes; ++j){

      F[i] += MatrixKernel.evaluate(rootsX[i], rootsY[j]) * W[j];

    }
  }
  std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;


  ////////////////////////////////////////////////////////////////////////////
  // Store M2L in K and apply K
  FReal K[nnodes*nnodes]; // local expansion
  for (unsigned int i=0; i<nnodes; ++i) {
    for (unsigned int j=0; j<nnodes; ++j){
      K[i*nnodes+j] = MatrixKernel.evaluate(rootsX[i], rootsY[j]);
    }
  }
  std::cout<< "Apply M2L in usual sense: ";
  time.tic();

  for (unsigned int i=0; i<nnodes; ++i) {
    F[i] = FReal(0.);
    for (unsigned int j=0; j<nnodes; ++j){

      F[i] += K[i*nnodes+j] * W[j];

    }
  }

  time.tac();
  std::cout << "took " << time.elapsed() << "sec." << std::endl;

  // PB: Verify storage improvement works (indexing etc...)
  // 1) store circulant matrix
  const unsigned int rc = (2*ORDER-1)*(2*ORDER-1)*(2*ORDER-1);
  FReal C[rc];

  //  std::cout << "2*order-1=" << 2*ORDER-1 << std::endl;
  unsigned int li,lj;
  unsigned int mi,mj;
  unsigned int ni,nj;
  unsigned int idi, idj, ido=0;

  for(unsigned int l=0; l<2*ORDER-1; ++l)
    for(unsigned int m=0; m<2*ORDER-1; ++m)
      for(unsigned int n=0; n<2*ORDER-1; ++n){

        // compute position in nnodes vector of roots
        if(l<ORDER-1) lj=ORDER-1; else lj=0;
        if(m<ORDER-1) mj=ORDER-1; else mj=0;
        if(n<ORDER-1) nj=ORDER-1; else nj=0;
        li=(l-(ORDER-1))+lj; mi=(m-(ORDER-1))+mj; ni=(n-(ORDER-1))+nj;

        idi=li*ORDER*ORDER + mi*ORDER + ni;
        idj=lj*ORDER*ORDER + mj*ORDER + nj;

        C[ido]
          = MatrixKernel.evaluate(rootsX[idi], rootsY[idj]);
        ido++;
      }
//  // Display C (gathers every values of K that need to be stored,
//  // corresponds to the first line of the padded matrix (reverse order?))
//  std::cout<<"C="<<std::endl;
//  for (unsigned int n=0; n<rankpad; ++n)
//    std::cout<< C[ido] << ", ";
//  std::cout<<std::endl;

  typedef FUnifTensor<ORDER> TensorType;
  unsigned int node_diff[nnodes*nnodes];
  TensorType::setNodeIdsDiff(node_diff);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // K is a block Toeplitz matrix
  // i.e. a blockwise Toeplitz matrix where the block also have the Toeplitz structure.
  // e.g. for ORDER=3: K=[K_{1,1} K_{1,2} K_{1,3},  where K_{i,j}=[k11 k12 k13,
  //                      K_{2,1} K_{1,1} K_{1,2},                 k21 k11 k12,
  //                      K_{3,1} K_{2,1} K_{1,1}];                k31 k21 k11];
  // K is of size order^3 x order^3
  // (i.e. order^2 x order^2 Toeplitz blocks of size order x order),
  // K is very close to be Toeplitz itself and even circulant.
  // In order to actually embed K into a circulant matrix C one just
  // needs to insert (ORDER-1) extra lines/columns (to each block).

//  std::cout<< "K=" <<std::endl;
//  for (unsigned int i=0; i<nnodes; ++i){
//    for (unsigned int j=0; j<nnodes; ++j){
//      std::cout<< K[i*nnodes+j]<<", ";
//    }
//    std::cout<<std::endl;
//  }
//  std::cout<<std::endl;

//  // Check matrix node_diff
//  std::cout<< "node_diff=" <<std::endl;
//  for (unsigned int i=0; i<nnodes; ++i){
//    for (unsigned int j=0; j<nnodes; ++j){
//      std::cout<< node_diff[i*nnodes+j] <<", ";
//    }
//    std::cout<<std::endl;
//  }
//  std::cout<<std::endl;

//  // Expected ido for the (2*ORDER-1)^3x(2*ORDER-1)^3 circulant matrix
//  for (unsigned int i=0; i<rc; ++i){
//    for (unsigned int j=0; j<rc; ++j){
//      if(i>j) std::cout<< i-j-1 << ", ";
//      else std::cout<< rc+i-j-1 << ", ";
//    } std::cout<<std::endl;
//  } std::cout<<std::endl;

//  // kernel evaluated at previous ido returns a circulant matrix
//  for (unsigned int i=0; i<rc/2; ++i){
//    for (unsigned int j=0; j<rc/2; ++j){
//      if(i>j) std::cout<< C[i-j-1] << ", ";
//      else std::cout<< C[rc+i-j-1] << ", ";
//    } std::cout<<std::endl;
//  } std::cout<<std::endl;

  // In 1D the Zero Padding consists in
  // inserting ORDER-1 zeros in the multipole exp
  // in order to apply the (ORDER+ORDER-1)x(ORDER+ORDER-1)
  // circulant matrix to it.
  // Let us extend it to the 3D case:
  FReal MultExp[nnodes]; FReal PaddedMultExp[rc];
  for (unsigned int i=0; i<nnodes; ++i) MultExp[i]=W[i];
  FReal LocalExp[nnodes]; FReal PaddedLocalExp[rc];
  FBlas::setzero(nnodes,LocalExp);
  FBlas::setzero(rc,PaddedLocalExp);

  std::cout<< "Expected LocalExp: "<<std::endl;
  for (unsigned int i=0; i<nnodes; ++i)
    std::cout<< F[i] << ", ";
  std::cout<<std::endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // Application of circulant Toeplitz system in PHYSICAL SPACE
  std::cout<< "Apply circulant M2L in Physical space: ";
  time.tic();
  for (unsigned int i=0; i<nnodes; ++i){

    // Pad Multipole Expansion with respect to current row
    FBlas::setzero(rc,PaddedMultExp);
    for (unsigned int j=0; j<nnodes; ++j)
      PaddedMultExp[node_diff[i*nnodes+j]]=MultExp[j];

//    std::cout<< "Padded MultExp for row i=" << i << ": "<<std::endl;
//    for (unsigned int p=0; p<rc; ++p)
//      std::cout<< PaddedMultExp[p] << ", ";
//    std::cout<<std::endl;

    // Application of M2L in PHYSICAL SPACE
    for (unsigned int pj=0; pj<rc; ++pj)
      LocalExp[i]+=C[pj]*PaddedMultExp[pj];

  }// end i
  time.tac();
  std::cout << "took " << time.elapsed() << "sec." << std::endl;

  std::cout<< "LocalExp via product in PHYSICAL SPACE: "<<std::endl;
  for (unsigned int p=0; p<nnodes; ++p)
    std::cout<< LocalExp[p] << ", ";
  std::cout<<std::endl;
  std::cout<<std::endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // Efficient application of the Toeplitz system in FOURIER SPACE
  FComplexe FPMultExp[rc];
  FComplexe FPLocalExp[rc];
  FReal PLocalExp[rc];

  //for (unsigned int n=0; n<rc; ++n) FPLocalExp[n]=FComplexe(0.0,0.0);
  FBlas::c_setzero(rc,reinterpret_cast<FReal*>(FPLocalExp));

  FBlas::setzero(rc,PLocalExp);

  // Init DFT
  //FDft Dft(rc); // direct version
  FFft Dft(rc); // fast version

  // Get first COLUMN of K and Store in T
  FReal T[rc];
  // use permutations
  unsigned int perm[rc];
  for(unsigned int p=0; p<rc; ++p){
    if(p>0) perm[p]=p-1;
    else perm[p]=rc+p-1;
//    std::cout << "perm["<< p << "]="<< perm[p] << std::endl;
  }

  for (unsigned int i=0; i<rc; ++i){
    // keep this lines commented to see what permutation accounts for:
    //  for (unsigned int j=0; j<rc; ++j){
//      if(i>0) T[i]=C[i-0-1];
//      else T[i]=C[rc+i-0-1];
      T[i]=C[perm[i]];
  //  }
  }

//  std::cout<< "First column of C[rc x rc]: "<<std::endl;
//  for (unsigned int p=0; p<rc; ++p)
//    std::cout<< T[p] << ", ";
//  std::cout<<std::endl;

  // Apply DFT to T
  FComplexe FT[rc];
  //  for (unsigned int n=0; n<rc; ++n) FT[n]=FComplexe(0.0,0.0);
  FBlas::c_setzero(rc,reinterpret_cast<FReal*>(FT));

  // if first COLUMN (T) of C is used
  Dft.applyDFT(T,FT);
//  // if first ROW of C is used
//  Dft.applyDFT(C,FT);

  // Pad physical MultExp
  FBlas::setzero(rc,PaddedMultExp); //part of padding
  for (unsigned int j=0; j<nnodes; ++j){
    // if first COLUMN (T) of C is used
    PaddedMultExp[node_diff[j*nnodes]]=MultExp[j];
//    // if first ROW of C is used
//    PaddedMultExp[node_diff[j]]=MultExp[j];
  }

//    std::cout<< "Padded MultExp for row i=" << i << ": "<<std::endl;
//    for (unsigned int p=0; p<rc; ++p)
//      std::cout<< PaddedMultExp[p] << ", ";
//    std::cout<<std::endl;


  // Set transformed MultExp to 0
  //  for (unsigned int n=0; n<rc; ++n) FPMultExp[n]=FComplexe(0.0,0.0);
  FBlas::c_setzero(rc,reinterpret_cast<FReal*>(FPMultExp));

  // Transform PaddedMultExp
  Dft.applyDFT(PaddedMultExp,FPMultExp);

  std::cout<< "Apply M2L in Fourier space: ";
  time.tic();

  // Application of M2L in FOURIER SPACE
  // > Use FMkl::c_had for hadamard product
  // if mkl is used as blas (TODO otherwise use FBlas::c_had())
//  FMkl::c_had(rc,reinterpret_cast<FReal*>(FT),
//            reinterpret_cast<FReal*>(FPMultExp),
//            reinterpret_cast<FReal*>(FPLocalExp));
  // > or perform entrywise product manually
  for (unsigned int pj=0; pj<rc; ++pj){
    FPLocalExp[pj]=FT[pj];
    FPLocalExp[pj]*=FPMultExp[pj];
  }

  time.tac();
  std::cout << "took " << time.elapsed() << "sec." << std::endl;

//    std::cout<< "Transfo Padded LocalExp: "<<std::endl;
//    for (unsigned int p=0; p<rc; ++p)
//      std::cout<< FPLocalExp[p] << ", ";
//    std::cout<<std::endl;

  Dft.applyIDFT(FPLocalExp,PLocalExp);

//  std::cout<< "Padded LocalExp: "<<std::endl;
//  for (unsigned int p=0; p<rc; ++p)
//    std::cout<< PLocalExp[p] << ", ";
//  std::cout<<std::endl;

  // Unpad
  for (unsigned int j=0; j<nnodes; ++j){
    // if first COLUMN (T) of C is used
    LocalExp[j]=PLocalExp[node_diff[nnodes-j-1]];
//    // if first ROW of C is used
//    LocalExp[j]=PLocalExp[node_diff[j*nnodes]];
  }

//  std::cout<< "Mask to be applied to Padded LocalExp: "<<std::endl;
//  for (unsigned int j=0; j<nnodes; ++j)
//    std::cout<< node_diff[nnodes-j-1] << ", ";
//  std::cout<<std::endl;

  std::cout<< "LocalExp via product in FOURIER SPACE: "<<std::endl;
  for (unsigned int p=0; p<nnodes; ++p)
    std::cout<< LocalExp[p] << ", ";
  std::cout<<std::endl;


  /////////////////////////////////////////////////////////////////////////////////////

  std::cout << "L2P (potential) ... " << std::flush;
  time.tic();
  // Interpolate p_i = \sum_m^L S(x_i,\bar x_m) * F_m
  S.applyL2P(cx, width, F, X.getTargets());
  std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

  std::cout << "L2P (forces) ... " << std::flush;
  time.tic();
  // Interpolate f_i = \sum_m^L P(x_i,\bar x_m) * F_m
  S.applyL2PGradient(cx, width, F, X.getTargets());
  std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

  ////////////////////////////////////////////////////////////////////
  // direct computation
  std::cout << "Compute interactions directly ..." << std::endl;
  time.tic();

  FReal* approx_f = new FReal [M * 3];
  FReal*        f = new FReal [M * 3];
  FBlas::setzero(M*3, f);

  FReal* approx_p = new FReal[M];
  FReal*        p = new FReal[M];
  FBlas::setzero(M, p);

  { // start direct computation
    unsigned int counter = 0;

    for(int idxPartX = 0 ; idxPartX < X.getSrc()->getNbParticles() ; ++idxPartX){
      const FPoint x = FPoint(X.getSrc()->getPositions()[0][idxPartX],
                              X.getSrc()->getPositions()[1][idxPartX],
                              X.getSrc()->getPositions()[2][idxPartX]);
      const FReal  wx = X.getSrc()->getPhysicalValues()[idxPartX];

      for(int idxPartY = 0 ; idxPartY < Y.getSrc()->getNbParticles() ; ++idxPartY){
        const FPoint y = FPoint(Y.getSrc()->getPositions()[0][idxPartY],
                                Y.getSrc()->getPositions()[1][idxPartY],
                                Y.getSrc()->getPositions()[2][idxPartY]);
        const FReal  wy = Y.getSrc()->getPhysicalValues()[idxPartY];

        const FReal one_over_r = MatrixKernel.evaluate(x, y);
        // potential
        p[counter] += one_over_r * wy;
        // force
        FPoint force(y - x);
        force *= one_over_r*one_over_r*one_over_r;
        f[counter*3 + 0] += force.getX() * wx * wy;
        f[counter*3 + 1] += force.getY() * wx * wy;
        f[counter*3 + 2] += force.getZ() * wx * wy;
      }

      counter++;
    }
  } // end direct computation


  time.tac();
  std::cout << "Done in " << time.elapsed() << "sec." << std::endl;


  ////////////////////////////////////////////////////////////////////
  unsigned int counter = 0;
  for(int idxPartX = 0 ; idxPartX < X.getSrc()->getNbParticles() ; ++idxPartX){
    approx_p[counter] = X.getSrc()->getPotentials()[idxPartX];
    const FPoint force = FPoint(X.getSrc()->getForcesX()[idxPartX],
                                X.getSrc()->getForcesY()[idxPartX],
                                X.getSrc()->getForcesZ()[idxPartX]);
    approx_f[counter*3 + 0] = force.getX();
    approx_f[counter*3 + 1] = force.getY();
    approx_f[counter*3 + 2] = force.getZ();

    counter++;
  }

  std::cout << "\nPotential error:" << std::endl;
  std::cout << "Relative error   = " << FMath::FAccurater( p, approx_p, M) << std::endl;

  std::cout << "\nForce error:" << std::endl;
  std::cout << "Relative L2 error   = " << FMath::FAccurater( f, approx_f, M*3) << std::endl;
  std::cout << std::endl;

  // free memory
  delete [] approx_p;
  delete [] p;
  delete [] approx_f;
  delete [] f;


  return 0;
}
