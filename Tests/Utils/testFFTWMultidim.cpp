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
#include <stdlib.h>

// include for libfftw3
//#include <fftw3.h>

// include for mkl_fftw3
#include <fftw/fftw3.h>

#include "../../Src/Utils/FGlobal.hpp"
#include "../../Src/Utils/FComplexe.hpp"

#include "../../Src/Utils/FTic.hpp"


//#include "../../Src/Utils/FDft.hpp"


int main()
{
  const FReal FRandMax = FReal(RAND_MAX);

  FTic time;

  //////////////////////////////////////////////////////////////////////////////
  // INITIALIZATION

  // size (pick a power of 2 for better performance of the FFT algorithm)
  unsigned int rank = 2;
  unsigned int nsteps_ = 500; 
  unsigned int dim = 10;
  const int steps[2]={static_cast<int>(dim),
                      static_cast<int>(nsteps_)};
  unsigned int size = dim*nsteps_;

  //////////////////////////////////////////////////////////////////////////////
  // Multidimensionnal FFT PLANS
  // fftw arrays
  FReal* fftR_;
  FComplexe* fftC_;

  fftR_ = (FReal*) fftw_malloc(sizeof(FReal) * size );
  fftC_ = (FComplexe*) fftw_malloc(sizeof(FComplexe) * size );

  // fftw plans
  // use routine defined in file:
  // /PATH/TO/mkl/interfaces/fftw3xf/wrappers/fftw_plan_dft_c2r_1d.c
  fftw_plan plan_c2r_; // backward FFT plan
  fftw_plan plan_r2c_; // forward FFT plan

  std::cout<< "Init FFTW plans: ";
  time.tic();

  plan_c2r_ =
    fftw_plan_dft_c2r(rank, steps, 
                      reinterpret_cast<fftw_complex*>(fftC_),
                      fftR_, 
                      FFTW_MEASURE);
  plan_r2c_ =
    fftw_plan_dft_r2c(rank, steps, 
                      fftR_, 
                      reinterpret_cast<fftw_complex*>(fftC_), 
                      FFTW_MEASURE);

  std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

  //////////////////////////////////////////////////////////////////////////////
  // EXECUTION
  // generate random physical data
  for(unsigned int s=0; s<size; ++s)
    fftR_[s] = FReal(rand())/FRandMax; 

// // display data in  physical space
// std::cout<< "Physical data: "<<std::endl;
// for(unsigned int d=0; d<dim; ++d){
//   for(unsigned int s=0; s<nsteps_; ++s)
//     std::cout<< fftR_[s+d*nsteps_] << ", ";
//   std::cout<<std::endl;
// }
// std::cout<<std::endl;

  // perform fft
  std::cout<< "Perform Forward FFT: ";
  time.tic();
  fftw_execute( plan_r2c_ );
  std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

//  // display transform in Fourier space
//  // beware the real data FFT stores only N/2+1 complex output values
//  std::cout<< "Transformed data : "<<std::endl;
//  for(unsigned int d=0; d<dim; ++d){
//    for(unsigned int s=0; s<nsteps_/2+1; ++s)
//      std::cout<< fftC_[s+d*(nsteps_/2+1)] << ", ";
//    std::cout<<std::endl;
//  }
//  std::cout<<std::endl;

//  for(unsigned int s=0; s<nsteps_/2+1; ++s){
//    fftC_[nsteps_-s]=FComplexe(fftC_[s].getReal(),-fftC_[s].getImag());
//  }
//
//  std::cout<< "Full Transformed data : "<<std::endl;
//  for(unsigned int s=0; s<nsteps_; ++s)
//    std::cout<< fftC_[s] << ", ";
//  std::cout<<std::endl;

  // perform ifft of tranformed data (in order to get physical values back)
  std::cout<< "Perform Backward FFT: ";
  time.tic();
  fftw_execute( plan_c2r_ );
  std::cout << "took " << time.tacAndElapsed() << "sec." << std::endl;

//  // display data in physical space
//  std::cout<< "Physical data (from 1/N*IFFT(FFT(Physical data))): "<<std::endl;
//  for(unsigned int d=0; d<dim; ++d){
//    for(unsigned int s=0; s<nsteps_; ++s)
//      std::cout<< fftR_[s+d*nsteps_]/(nsteps_*dim) << ", ";
//    std::cout<<std::endl;
//  }
//  std::cout<<std::endl;

  //free memory
  fftw_destroy_plan(plan_r2c_);
  fftw_destroy_plan(plan_c2r_);
  fftw_free(fftR_);
  fftw_free(fftC_);

}// end test
