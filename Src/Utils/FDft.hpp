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
#ifndef FDFT_HPP
#define FDFT_HPP

#include <iostream>
#include <stdlib.h>

// include fftw3 (specify path in cmake)
// if MKL: path/to/mkl/include/fftw/fftw3.h 
// elseif libfftw_dev: usr/include/fftw3.h
#include <fftw3.h>

#include "FGlobal.hpp"
#include "FComplexe.hpp"

#include "FMath.hpp"

#include "FTic.hpp"


/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * @class FDft and @class FFft
 * Please read the license
 *
 * These classes handle the forward and backward Discete Fourier Transform 
 * (DFT).
 * @class FDft implements a direct method while @class FFft uses the Fast 
 * Fourier Transform (FFT). The FFT algorithm can either be provided by the 
 * FFTW(3) library itself or a version that is wrapped in Intel MKL.
 *
 * The aim of writing these specific classes is to allow further customization 
 * of the DFT such as implementing a tensorial variant, a weighted variant 
 * or any other feature.
 *
 * TODO: some copies might slow down the routines, this needs to be optimized.
 * 
 */


/**
 * @class FDft
 *
 * @tparam nsteps number of sampled values \f$N\f$
 */

class FDft
{

  FReal* fftR_;
  FComplexe* fftC_;

  FReal *cosRS_, *sinRS_; 

private:
  unsigned int nsteps_;

public:

  FDft(const unsigned int nsteps)
  : nsteps_(nsteps)
  {
    fftR_ = new FReal[nsteps_];
    fftC_ = new FComplexe[nsteps_];

    // Beware this is extremely HEAVY to store!!! => Use FDft only for debug!
    cosRS_ = new FReal[nsteps_*nsteps_];
    sinRS_ = new FReal[nsteps_*nsteps_];

    for(unsigned int r=0; r<nsteps_; ++r)
      for(unsigned int s=0; s<nsteps_; ++s){
        FReal thetaRS = 2*M_PI*r*s/nsteps_;
        cosRS_[r*nsteps_+s]=FMath::Cos(thetaRS);
        sinRS_[r*nsteps_+s]=FMath::Sin(thetaRS);
      }
  }

  virtual ~FDft()
  {
    delete [] fftR_;
    delete [] fftC_;
    delete [] cosRS_;
    delete [] sinRS_;
  }

  void applyDFT(const FReal* sampledData,
                FComplexe* transformedData) const
  {
//    FTic time;

    // read sampled data
//    std::cout<< "copy(";
//    time.tic();
    FBlas::c_setzero(nsteps_,reinterpret_cast<FReal*>(fftC_));
    FBlas::copy(nsteps_, sampledData,fftR_);
//    std::cout << time.tacAndElapsed() << ")";

//    std::cout<< " - exe(";
//    time.tic();
    for(unsigned int r=0; r<nsteps_; ++r)
      for(unsigned int s=0; s<nsteps_; ++s){
        fftC_[r] += FComplexe(fftR_[s]*cosRS_[r*nsteps_+s], 
                              -fftR_[s]*sinRS_[r*nsteps_+s]);
      }
//    std::cout << time.tacAndElapsed() << ")";

    // write transformed data
//    std::cout<< " - copy(";
//    time.tic();
    FBlas::c_copy(nsteps_,reinterpret_cast<FReal*>(fftC_),
                  reinterpret_cast<FReal*>(transformedData));
//    std::cout << time.tacAndElapsed() << ") ";

  }


  void applyIDFT(const FComplexe* transformedData,
                 FReal* sampledData) const
  {
    // read transformed data
    FBlas::setzero(nsteps_,fftR_);
    FBlas::c_copy(nsteps_,reinterpret_cast<const FReal*>(transformedData),
                  reinterpret_cast<FReal*>(fftC_));

    for(unsigned int r=0; r<nsteps_; ++r){
      for(unsigned int s=0; s<nsteps_; ++s){
        fftR_[r] += fftC_[s].getReal()*cosRS_[r*nsteps_+s]
          -fftC_[s].getImag()*sinRS_[r*nsteps_+s];
      }
      fftR_[r]*=1./nsteps_;
    }  

    // write sampled data
    FBlas::copy(nsteps_,fftR_,sampledData);

  }


};



/**
 * @class FFft
 *
 * @tparam nsteps number of sampled values \f$N\f$
 */

class FFft
{
  // arrays
  FReal* fftR_;
  FComplexe* fftC_;

  // plans
  fftw_plan plan_c2r_; // backward FFT plan
  fftw_plan plan_r2c_; // forward FFT plan 

private:
  unsigned int nsteps_;
  unsigned int nsteps_opt_;

public:

  FFft(const unsigned int nsteps)
  : nsteps_(nsteps),
    nsteps_opt_(nsteps/2+1) // SPECIFIC TO FTT FOR REAL VALUES
  {
    // allocate arrays
    fftR_ = (FReal*) fftw_malloc(sizeof(FReal) * nsteps_);
    fftC_ = (FComplexe*) fftw_malloc(sizeof(FComplexe) * nsteps_opt_);

    // fftw plans
    plan_c2r_ =
      fftw_plan_dft_c2r_1d(nsteps_, 
                           reinterpret_cast<fftw_complex*>(fftC_),
                           fftR_, 
                           FFTW_MEASURE);// TODO: test FFTW_ESTIMATE
    plan_r2c_ =
      fftw_plan_dft_r2c_1d(nsteps_, 
                           fftR_, 
                           reinterpret_cast<fftw_complex*>(fftC_), 
                           FFTW_MEASURE);


  }

  virtual ~FFft()
  {
    fftw_destroy_plan(plan_c2r_);
    fftw_destroy_plan(plan_r2c_);
    fftw_free(fftR_);
    fftw_free(fftC_);
  }

  void applyDFT(const FReal* sampledData,
                FComplexe* transformedData) const
  {
//    FTic time;

    // read sampled data
//    std::cout<< "copy(";
//    time.tic();
    FBlas::c_setzero(nsteps_opt_,reinterpret_cast<FReal*>(fftC_));
    FBlas::copy(nsteps_, sampledData,fftR_);
//    std::cout << time.tacAndElapsed() << ")";

    // perform fft
//    std::cout<< " - exe(";
//    time.tic();
    fftw_execute( plan_r2c_ );
//    std::cout << time.tacAndElapsed() << ")";

    // write transformed data
//    std::cout<< " - copy(";
//    time.tic();
    FBlas::c_copy(nsteps_opt_,reinterpret_cast<FReal*>(fftC_),
                  reinterpret_cast<FReal*>(transformedData));
//    for(unsigned int s=0; s<nsteps_opt_; ++s)
//      transformedData[s]=fftC_[s];

//    std::cout << time.tacAndElapsed() << ") ";

  }


  void applyIDFT(const FComplexe* transformedData,
                 FReal* sampledData) const
  {
    // read transformed data
    FBlas::setzero(nsteps_,fftR_);
    FBlas::c_copy(nsteps_opt_,reinterpret_cast<const FReal*>(transformedData),
                  reinterpret_cast<FReal*>(fftC_));

    // perform ifft
    fftw_execute( plan_c2r_ );

    for(unsigned int s=0; s<nsteps_; ++s)
      fftR_[s]/=nsteps_; // the fft from fftw is not scaled !!!!!!!!

    // write sampled data
    FBlas::copy(nsteps_,fftR_,sampledData);

  }


};

#endif /* FDFT_HPP */

