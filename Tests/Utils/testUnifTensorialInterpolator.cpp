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
// @FUSE_FFT
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
//  typedef FInterpMatrixKernelR MatrixKernelClass;
  typedef FInterpMatrixKernel_R_IJ MatrixKernelClass;
  typedef FInterpMatrixKernel_R_IJK RIJKMatrixKernelClass; // PB: force computation


  ///////////////////////What we do/////////////////////////////
  std::cout << "\nTask: Compute interactions between source particles in leaf Y and target\n";
  std::cout << " particles in leaf X. Compare the fast summation K ~ Lx K Ly' with the\n";
  std::cout << " direct computation.\n" << std::endl;
  //////////////////////////////////////////////////////////////

  const unsigned int dim = MatrixKernelClass::DIM;
  const unsigned int nrhs = MatrixKernelClass::NRHS;
  const unsigned int nlhs = MatrixKernelClass::NLHS;

  const FReal FRandMax = FReal(RAND_MAX);
  FTic time;


  // Leaf size
  FReal width = FReal(3.723);

  ////////////////////////////////////////////////////////////////////
  LeafClass X;
  FPoint cx(0., 0., 0.);
  const unsigned long M = 5000;
  std::cout << "Fill the leaf X of width " << width
            << " centered at cx=" << cx << " with M=" << M << " target particles" << std::endl;
  {
    for(unsigned long i=0; i<M; ++i){
      FReal x = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getX();
      FReal y = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getY();
      FReal z = (FReal(rand())/FRandMax - FReal(.5)) * width + cx.getZ();
      // TODO: find a way to push vectorial attribute of arbitrary dim (nlhs)
      X.push(FPoint(x, y, z), FReal(rand())/FRandMax);
    }
  }


  ////////////////////////////////////////////////////////////////////
  LeafClass Y;
  //  FPoint cy(FReal(2.)*width, 0., 0.);
  FPoint cy(FReal(2.)*width, FReal(2.)*width, 0.);

  const unsigned long N = 5000;
  std::cout << "Fill the leaf Y of width " << width
            << " centered at cy=" << cy	<< " with N=" << N << " source particles" << std::endl;
  {
    for(unsigned long i=0; i<N; ++i){
      FReal x = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getX();
      FReal y = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getY();
      FReal z = (FReal(rand())/FRandMax - FReal(.5)) * width + cy.getZ();
      // TODO: find a way to push vectorial attribute of arbitrary dim (nrhs)
      Y.push(FPoint(x, y, z), FReal(rand())/FRandMax);
    }
  }



  ////////////////////////////////////////////////////////////////////
  // approximative computation
  const unsigned int ORDER = 4;
  const unsigned int nnodes = TensorTraits<ORDER>::nnodes;
  typedef FUnifInterpolator<ORDER,MatrixKernelClass> InterpolatorClass;
  InterpolatorClass S;
  MatrixKernelClass MatrixKernel;
  RIJKMatrixKernelClass RIJKMatrixKernel;

  std::cout << "\nCompute interactions approximatively, interpolation order = " << ORDER << " ..." << std::endl;

  std::cout << "\nP2M ... " << std::flush;
  time.tic();
  // Anterpolate: W_n = \sum_j^N S(y_j,\bar y_n) * w_j
  FReal W[nrhs*nnodes]; // multipole expansion
  // tensorial case interpolate same Y for each component
  for (unsigned int idRhs=0; idRhs<nrhs; ++idRhs)
    S.applyP2M(cy, width, W + idRhs*nnodes, Y.getSrc()); // the multipole expansions are set to 0 in S.applyP2M
  std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

  std::cout << "M2L ... " << std::flush;
  time.tic();
  // Multipole to local: F_m = \sum_n^L K(\bar x_m, \bar y_n) * W_n
  FPoint rootsX[nnodes], rootsY[nnodes];
  FUnifTensor<ORDER>::setRoots(cx, width, rootsX);
  FUnifTensor<ORDER>::setRoots(cy, width, rootsY);

  FReal F[nlhs*nnodes]; // local expansion
  for (unsigned int i=0; i<nnodes*nlhs; ++i) F[i] = FReal(0.);

  for (unsigned int i=0; i<nnodes; ++i) {
    for (unsigned int j=0; j<nnodes; ++j){
      
      for (unsigned int idxLhs=0; idxLhs<nlhs; ++idxLhs)
        for (unsigned int idxRhs=0; idxRhs<nrhs; ++idxRhs){
          unsigned idxK = idxLhs*3+idxRhs; // or counter
          unsigned int d = MatrixKernel.getPosition(idxK);

          F[i+idxLhs*nnodes] += MatrixKernelClass(d).evaluate(rootsX[i], rootsY[j]) * W[j+idxRhs*nnodes];

        }
    }
  }
  std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

//  std::cout<< "F via direct applyM2L: "<<std::endl;
//  for (unsigned int d=0; d<nlhs; ++d){
//    for (unsigned int i=0; i<nnodes; ++i)
//      std::cout<< F[i+d*nnodes] << ", ";
//    std::cout<<std::endl;
//  }
//  std::cout<<std::endl;

  ////////////////////////////////////////////////////////////////////////////
  // Store M2L in K and apply K
  FReal K[dim*nnodes*nnodes]; // local expansion
  for (unsigned int i=0; i<nnodes; ++i) {
    for (unsigned int j=0; j<nnodes; ++j){

      for (unsigned int d=0; d<dim; ++d){
        K[d*nnodes*nnodes+i*nnodes+j] = MatrixKernelClass(d).evaluate(rootsX[i], rootsY[j]);        
      }

    }
  }
  std::cout<< "Apply M2L in usual sense: ";
  time.tic();
  for (unsigned int i=0; i<nnodes*nlhs; ++i) F[i] = FReal(0.);

  for (unsigned int i=0; i<nnodes; ++i)
    for (unsigned int j=0; j<nnodes; ++j){

      for (unsigned int idxLhs=0; idxLhs<nlhs; ++idxLhs)
        for (unsigned int idxRhs=0; idxRhs<nrhs; ++idxRhs){
          unsigned idxK = idxLhs*3+idxRhs; // or counter
          unsigned int d = MatrixKernel.getPosition(idxK);

          F[i+idxLhs*nnodes] += K[d*nnodes*nnodes+i*nnodes+j] * W[j+idxRhs*nnodes];

    }
  }

  time.tac();
  std::cout << "took " << time.elapsed() << "sec." << std::endl;

//  std::cout<< "F via store and applyM2L: "<<std::endl;
//  for (unsigned int d=0; d<nlhs; ++d){
//    for (unsigned int i=0; i<nnodes; ++i)
//      std::cout<< F[i+d*nnodes] << ", ";
//    std::cout<<std::endl;
//  }
//  std::cout<<std::endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // PB: Verify storage improvement works (indexing etc...)
  // 1) store circulant matrix
  const unsigned int rc = (2*ORDER-1)*(2*ORDER-1)*(2*ORDER-1);
  FReal C[dim*rc];

  typedef FUnifTensor<ORDER> TensorType;
  unsigned int node_diff[nnodes*nnodes];
  TensorType::setNodeIdsDiff(node_diff);
  unsigned int node_ids_pairs[rc][2];
  TensorType::setNodeIdsPairs(node_ids_pairs);

  unsigned int ido=0;

  for(unsigned int l=0; l<2*ORDER-1; ++l)
    for(unsigned int m=0; m<2*ORDER-1; ++m)
      for(unsigned int n=0; n<2*ORDER-1; ++n){

        for (unsigned int d=0; d<dim; ++d){

          C[d*rc + ido]
            = MatrixKernelClass(d).evaluate(rootsX[node_ids_pairs[ido][0]], 
                                    rootsY[node_ids_pairs[ido][1]]);
        }
        
        ido++;
      }

//  // Display C (gathers every values of K that need to be stored,
//  // corresponds to the first line of the padded matrix (reverse order?))
//  std::cout<<"C="<<std::endl;
//    for (unsigned int d=0; d<dim; ++d){
//      for (unsigned int n=0; n<rc; ++n)
//        std::cout<< C[n + d*rc] << ", ";
//      std::cout<<std::endl;
//    }
//  std::cout<<std::endl;

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
  FReal MultExp[nrhs*nnodes]; FReal PaddedMultExp[nrhs*rc];
  for (unsigned int i=0; i<nrhs*nnodes; ++i) MultExp[i]=W[i];
  FReal LocalExp[nlhs*nnodes]; FReal PaddedLocalExp[nlhs*rc];
  FBlas::setzero(nlhs*nnodes,LocalExp);
  FBlas::setzero(nlhs*rc,PaddedLocalExp);

//  std::cout<< "Expected LocalExp: "<<std::endl;
//  for (unsigned int d=0; d<nlhs; ++d){
//    for (unsigned int i=0; i<nnodes; ++i)
//      std::cout<< F[i+d*nnodes] << ", ";
//    std::cout<<std::endl;
//  }

  /////////////////////////////////////////////////////////////////////////////////////
  // Application of circulant Toeplitz system in PHYSICAL SPACE
  std::cout<< "Apply circulant M2L in Physical space: ";
  time.tic();
  for (unsigned int i=0; i<nnodes; ++i){

    // Pad Multipole Expansion with respect to current row
    FBlas::setzero(nrhs*rc,PaddedMultExp);
    for (unsigned int j=0; j<nnodes; ++j)
        for (unsigned int d=0; d<nrhs; ++d)
          PaddedMultExp[node_diff[i*nnodes+j] + d*rc]=MultExp[j + d*nnodes];

//    std::cout<< "Padded MultExp for row i=" << i << ": "<<std::endl;
//    for (unsigned int p=0; p<rc; ++p)
//      std::cout<< PaddedMultExp[p] << ", ";
//    std::cout<<std::endl;

    // Application of M2L in PHYSICAL SPACE
    for (unsigned int pj=0; pj<rc; ++pj)
      for (unsigned int idxLhs=0; idxLhs<nlhs; ++idxLhs)
        for (unsigned int idxRhs=0; idxRhs<nrhs; ++idxRhs){
          unsigned idxK = idxLhs*3+idxRhs; // or counter

          unsigned int d = MatrixKernel.getPosition(idxK);

          LocalExp[i + idxLhs*nnodes]+=C[pj + d*rc]*PaddedMultExp[pj + idxRhs*rc];

        }

  }// end i
  time.tac();
  std::cout << "took " << time.elapsed() << "sec." << std::endl;

//  std::cout<< "LocalExp via product in PHYSICAL SPACE: "<<std::endl;
//  for (unsigned int d=0; d<nlhs; ++d){
//    for (unsigned int i=0; i<nnodes; ++i)
//      std::cout<< LocalExp[i+d*nnodes] << ", ";
//    std::cout<<std::endl;
//  }
//  std::cout<<std::endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // Efficient application of the Toeplitz system in FOURIER SPACE

  // Init DFT
  //FDft Dft(rc); // direct version
  FFft<1/*(TODO: fix MultidimFFT) nK*/> DftK(rc); // fast version
  FFft<1/*nrhs*/> DftRhs(rc);
  FFft<1/*nlhs*/> DftLhs(rc);

  // Get first COLUMN of K and Store in T
  FReal T[dim*rc];
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
    for (unsigned int d=0; d<dim; ++d)
      T[i + d*rc]=C[perm[i] + d*rc];
  //  }
  }

//  std::cout<< "First column of C[rc x rc]: "<<std::endl;
//  for (unsigned int p=0; p<rc; ++p)
//    std::cout<< T[p] << ", ";
//  std::cout<<std::endl;

  // Apply DFT to T
  FComplexe FT[dim*rc];
  //  for (unsigned int n=0; n<rc; ++n) FT[n]=FComplexe(0.0,0.0);
  FBlas::c_setzero(dim*rc,reinterpret_cast<FReal*>(FT));

  // if first COLUMN (T) of C is used
  for (unsigned int d=0; d<dim; ++d)
    DftK.applyDFT(T+d*rc,FT+d*rc);
//  // if first ROW of C is used
//  DftK.applyDFT(C,FT);

  FComplexe FPMultExp[nrhs*rc];
  FComplexe FPLocalExp[nlhs*rc];
  FReal PLocalExp[nlhs*rc];

  //for (unsigned int n=0; n<rc; ++n) FPLocalExp[n]=FComplexe(0.0,0.0);
  FBlas::c_setzero(nlhs*rc,reinterpret_cast<FReal*>(FPLocalExp));

  FBlas::setzero(nlhs*rc,PLocalExp);

  // Pad physical MultExp
  FBlas::setzero(nrhs*rc,PaddedMultExp); //part of padding
  for (unsigned int idRhs=0; idRhs<nrhs; ++idRhs)
    for (unsigned int j=0; j<nnodes; ++j){
      // if first COLUMN (T) of C is used
      PaddedMultExp[node_diff[j*nnodes]+idRhs*rc]=MultExp[j+idRhs*nnodes];
//    // if first ROW of C is used
//    PaddedMultExp[node_diff[j]]=MultExp[j];
    }

//    std::cout<< "Padded MultExp: "<<std::endl;
//    for (unsigned int p=0; p<rc; ++p)
//      std::cout<< PaddedMultExp[p] << ", ";
//    std::cout<<std::endl;


  // Set transformed MultExp to 0
  //  for (unsigned int n=0; n<rc; ++n) FPMultExp[n]=FComplexe(0.0,0.0);
  FBlas::c_setzero(nrhs*rc,reinterpret_cast<FReal*>(FPMultExp));

  // Transform PaddedMultExp
  for (unsigned int idxRhs=0; idxRhs<nrhs; ++idxRhs) // apply nrhs 1 dimensionnal FFT
    DftRhs.applyDFT(PaddedMultExp+idxRhs*rc,FPMultExp+idxRhs*rc);

  std::cout<< "Apply M2L in Fourier space: ";
  time.tic();

  // Application of M2L in FOURIER SPACE
  // > Use FMkl::c_had for hadamard product
  // if mkl is used as blas (TODO otherwise use FBlas::c_had())
//  FMkl::c_had(rc,reinterpret_cast<FReal*>(FT),
//            reinterpret_cast<FReal*>(FPMultExp),
//            reinterpret_cast<FReal*>(FPLocalExp));
  // > or perform entrywise product manually
  FComplexe tmpFX;
  for (unsigned int idxLhs=0; idxLhs<nlhs; ++idxLhs)
    for (unsigned int idxRhs=0; idxRhs<nrhs; ++idxRhs){
      unsigned int idxK = idxLhs*3+idxRhs; // or counter

      unsigned int d = MatrixKernel.getPosition(idxK);

      for (unsigned int pj=0; pj<rc; ++pj){
        tmpFX=FT[pj + d*rc];
        tmpFX*=FPMultExp[pj+idxRhs*rc];
        FPLocalExp[pj+idxLhs*rc]+=tmpFX; // add new contribution +RijYj
      }

  }
  time.tac();
  std::cout << "took " << time.elapsed() << "sec." << std::endl;

//    std::cout<< "Transfo Padded LocalExp: "<<std::endl;
//    for (unsigned int p=0; p<rc; ++p)
//      std::cout<< FPLocalExp[p] << ", ";
//    std::cout<<std::endl;

  for (unsigned int idxLhs=0; idxLhs<nlhs; ++idxLhs) // apply nrhs 1 dimensionnal FFT
    DftLhs.applyIDFT(FPLocalExp+idxLhs*rc,PLocalExp+idxLhs*rc);

//  std::cout<< "Padded LocalExp: "<<std::endl;
//  for (unsigned int p=0; p<rc; ++p)
//    std::cout<< PLocalExp[p] << ", ";
//  std::cout<<std::endl;

  // Unpad
  for (unsigned int idLhs=0; idLhs<nlhs; ++idLhs)
    for (unsigned int j=0; j<nnodes; ++j){
      // if first COLUMN (T) of C is used
      LocalExp[j+idLhs*nnodes]=PLocalExp[node_diff[nnodes-j-1]+idLhs*rc];
//    // if first ROW of C is used
//    LocalExp[j]=PLocalExp[node_diff[j*nnodes]];
    }

//  std::cout<< "Mask to be applied to Padded LocalExp: "<<std::endl;
//  for (unsigned int j=0; j<nnodes; ++j)
//    std::cout<< node_diff[nnodes-j-1] << ", ";
//  std::cout<<std::endl;

//  std::cout<< "LocalExp via product in FOURIER SPACE: "<<std::endl;
//  for (unsigned int d=0; d<nlhs; ++d){
//    for (unsigned int p=0; p<nnodes; ++p)
//      std::cout<< LocalExp[p + d*nnodes] << ", ";
//    std::cout<<std::endl;
//  }
//  std::cout<<std::endl;


  /////////////////////////////////////////////////////////////////////////////////////

  std::cout << "L2P (potential) ... " << std::flush;
  time.tic();
  // Interpolate p_i = \sum_m^L S(x_i,\bar x_m) * F_m
  for (unsigned int idLhs=0; idLhs<nlhs; ++idLhs)
    S.applyL2P(cx, width, F+idLhs*nnodes, X.getTargets());
  std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

  std::cout << "L2P (forces) ... " << std::flush;
  time.tic();
  // Interpolate f_i = \sum_m^L P(x_i,\bar x_m) * F_m
  for (unsigned int idLhs=0; idLhs<nlhs; ++idLhs)
    S.applyL2PGradient(cx, width, F+idLhs*nnodes, X.getTargets());
  std::cout << "took " << time.tacAndElapsed() << "s" << std::endl;

  ////////////////////////////////////////////////////////////////////
  // direct computation 
  // Only scalar phys val => only compute first compo and no derivative
  // TODO add multidim phys val !!!!!!
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

//        // 1/R
//        const FReal one_over_r = MatrixKernel.evaluate(x, y);
//
//        // potential
//        p[counter] += one_over_r * wy;
//        // force
//        FPoint force(y - x);
//        force *= one_over_r*one_over_r*one_over_r;
//        f[counter*3 + 0] += force.getX() * wx * wy;
//        f[counter*3 + 1] += force.getY() * wx * wy;
//        f[counter*3 + 2] += force.getZ() * wx * wy;

        // R,ij and (R,ij),k
        for (unsigned int i=0; i<nlhs; ++i) // sum all compo
          for (unsigned int j=0; j<nrhs; ++j){
            unsigned int d = MatrixKernel.getPosition(i*nrhs+j);
            const FReal rij = MatrixKernelClass(d).evaluate(x, y);
            // potential
            p[counter] += rij * wy;
            // force
            FReal force[3];
            for (unsigned int k=0; k<3; ++k){
              //std::cout << "i,j,k,=" << i << ","<< j << ","<< k << std::endl;
              unsigned int dk = RIJKMatrixKernel.getPosition(i*3*3+j*3+k);
              force[k] = RIJKMatrixKernelClass(dk).evaluate(x, y);
              f[counter*3 + k] += force[k] * wx * wy;
            }
          }

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
