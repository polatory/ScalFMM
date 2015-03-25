// ===================================================================================
// Copyright ScalFmm 2011 INRIA
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
#include "ScalFmmConfig.h"

#ifdef ScalFMM_USE_FFT
#include <iostream>
#include <stdlib.h>

// include fftw3 (specify path in cmake)
// if MKL: path/to/mkl/include/fftw/fftw3.h 
// elseif libfftw_dev: usr/include/fftw3.h
#include <fftw3.h>

#include "FGlobal.hpp"
#include "FComplex.hpp"

#include "FMath.hpp"

#include "FTic.hpp"


/**
 * @author Pierre Blanchard (pierre.blanchard@inria.fr)
 * @class FDft, @class FFft and @class FCFft
 * Please read the license
 *
 * These classes handle the forward and backward Discete Fourier Transform
 * (DFT).
 * @class FDft implements a direct method while @class FFft uses the Fast
 * Fourier Transform (FFT). The FFT algorithm can either be provided by the
 * FFTW(3) library itself or a version that is wrapped in Intel MKL.
 *
 * The direct DFT is templatized with the input value type (FReal or FComplex<FReal>),
 * while 2 distinct classes resp. @class FFft and @class FCFft are defined
 * resp. for the real and complex valued FFT.
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
template< class FReal, typename ValueType = FReal>
class FDft
{

    ValueType* data_;  //< data in physical space
    FComplex<FReal>* dataF_; //< data in Fourier space

    FReal *cosRS_, *sinRS_;

private:
    unsigned int nsteps_; //< number of steps

    void initDFT()
    {
        // allocate arrays
        data_  = new ValueType[nsteps_];
        dataF_ = new FComplex<FReal>[nsteps_];

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

public:

    FDft(const unsigned int nsteps)
        : nsteps_(nsteps)
    {
        // init DFT
        initDFT();
    }

    FDft(const FDft& other)
        : nsteps_(other.nsteps_)
    {
        // init DFT
        initDFT();
    }

    virtual ~FDft()
    {
        delete [] data_;
        delete [] dataF_;
        delete [] cosRS_;
        delete [] sinRS_;
    }

    /// Forward DFT
    // Real valued DFT
    void applyDFT(const FReal* sampledData,
                  FComplex<FReal>* transformedData) const
    {
        // read sampled data
        FBlas::c_setzero(nsteps_,reinterpret_cast<FReal*>(dataF_));
        FBlas::copy(nsteps_, sampledData,data_);

        // perform direct forward transformation
        for(unsigned int r=0; r<nsteps_; ++r)
            for(unsigned int s=0; s<nsteps_; ++s){
                dataF_[r] += FComplex<FReal>(data_[s]*cosRS_[r*nsteps_+s],
                        -data_[s]*sinRS_[r*nsteps_+s]);
            }

        // write transformed data
        FBlas::c_copy(nsteps_,reinterpret_cast<FReal*>(dataF_),
                      reinterpret_cast<FReal*>(transformedData));
    }
    // Complexe valued DFT
    void applyDFT(const FComplex<FReal>* sampledData,
                  FComplex<FReal>* transformedData) const
    {
        // read sampled data
        FBlas::c_setzero(nsteps_,reinterpret_cast<FReal*>(dataF_));
        FBlas::c_copy(nsteps_,reinterpret_cast<const FReal*>(sampledData),
                      reinterpret_cast<FReal*>(data_));

        // perform direct forward transformation
        for(unsigned int r=0; r<nsteps_; ++r)
            for(unsigned int s=0; s<nsteps_; ++s){
                dataF_[r] += FComplex<FReal>(data_[s].getReal()*cosRS_[r*nsteps_+s]
                        + data_[s].getImag()*sinRS_[r*nsteps_+s],
                        data_[s].getImag()*cosRS_[r*nsteps_+s]
                        - data_[s].getReal()*sinRS_[r*nsteps_+s]);
            }

        // write transformed data
        FBlas::c_copy(nsteps_,reinterpret_cast<FReal*>(dataF_),
                      reinterpret_cast<FReal*>(transformedData));
    }

    /// Backward DFT
    // Real valued IDFT
    void applyIDFT(const FComplex<FReal>* transformedData,
                   FReal* sampledData) const
    {
        // read transformed data
        FBlas::setzero(nsteps_,data_);
        FBlas::c_copy(nsteps_,reinterpret_cast<const FReal*>(transformedData),
                      reinterpret_cast<FReal*>(dataF_));

        // perform direct backward transformation
        for(unsigned int r=0; r<nsteps_; ++r){
            for(unsigned int s=0; s<nsteps_; ++s){
                data_[r] += dataF_[s].getReal()*cosRS_[r*nsteps_+s]
                        + dataF_[s].getImag()*sinRS_[r*nsteps_+s];
            }
            data_[r]*=1./nsteps_;
        }

        // write sampled data
        FBlas::copy(nsteps_,data_,sampledData);
    }

    // Complexe valued IDFT
    void applyIDFT(const FComplex<FReal>* transformedData,
                   FComplex<FReal>* sampledData) const
    {
        // read transformed data
        FBlas::c_setzero(nsteps_,reinterpret_cast<FReal*>(data_));
        FBlas::c_copy(nsteps_,reinterpret_cast<const FReal*>(transformedData),
                      reinterpret_cast<FReal*>(dataF_));

        // perform direct backward transformation
        for(unsigned int r=0; r<nsteps_; ++r){
            for(unsigned int s=0; s<nsteps_; ++s){
                data_[r] += FComplex<FReal>(dataF_[s].getReal()*cosRS_[r*nsteps_+s]
                        - dataF_[s].getImag()*sinRS_[r*nsteps_+s],
                        dataF_[s].getImag()*cosRS_[r*nsteps_+s]
                        + dataF_[s].getReal()*sinRS_[r*nsteps_+s]);
            }
            data_[r]*=1./nsteps_;
        }

        // write sampled data
        FBlas::c_copy(nsteps_,reinterpret_cast<FReal*>(data_),
                      reinterpret_cast<FReal*>(sampledData));
    }

};


/**
 * @class FFft
 *
 * @tparam nsteps number of sampled values \f$N\f$
 */
template<class FReal, int DIM = 1>
class FFft
{
    enum{dim = DIM};

    // arrays
    FReal* data_;      //< data in physical space
    FComplex<FReal>* dataF_; //< data in Fourier space

    // plans
    fftw_plan plan_c2r_; //< backward FFT plan
    fftw_plan plan_r2c_; //< forward FFT plan

private:
    int steps_[dim]; //< number of steps per dimension
    int nsteps_;     //< total number of steps
    int nsteps_opt_; //< reduced number of steps for real valued FFT

    void initDFT()
    {
        // allocate arrays
        data_ = (FReal*) fftw_malloc(sizeof(FReal) * nsteps_);
        dataF_ = (FComplex<FReal>*) fftw_malloc(sizeof(FComplex<FReal>) * nsteps_opt_);

        //    // unidim fftw plans
        //    plan_c2r_ =
        //      fftw_plan_dft_c2r_1d(nsteps_,
        //                           reinterpret_cast<fftw_complex*>(dataF_),
        //                           data_,
        //                           FFTW_MEASURE);// TODO: test FFTW_ESTIMATE
        //    plan_r2c_ =
        //      fftw_plan_dft_r2c_1d(nsteps_,
        //                           data_,
        //                           reinterpret_cast<fftw_complex*>(dataF_),
        //                           FFTW_MEASURE);

        // multidim fftw plans
        plan_c2r_ =
                fftw_plan_dft_c2r(dim, steps_,
                                  reinterpret_cast<fftw_complex*>(dataF_),
                                  data_,
                                  FFTW_MEASURE);// TODO: test FFTW_ESTIMATE
        plan_r2c_ =
                fftw_plan_dft_r2c(dim, steps_,
                                  data_,
                                  reinterpret_cast<fftw_complex*>(dataF_),
                                  FFTW_MEASURE);
    }

public:

    FFft()
        : nsteps_(0), nsteps_opt_(0)
    {
        for(int d = 0; d<dim;++d){steps_[d]=0;}
    }

    void buildDFT(const int steps[dim])
    {
        // init number of steps
        nsteps_=1; nsteps_opt_=1;
        for(int d = 0; d<dim;++d){
            steps_[d]=steps[d];
            nsteps_*=steps[d];
            nsteps_opt_*=steps[d]/2+1; // real valued DFT specific
        }

        // init DFT
        initDFT();

    }

    FFft(const FFft& other)
        : nsteps_(other.nsteps_),nsteps_opt_(other.nsteps_opt_)
    {
        // mem copies
        memcpy(steps_,other.steps_,dim*sizeof(int));
        // init DFT
        initDFT();
    }

    virtual ~FFft()
    {
        fftw_destroy_plan(plan_c2r_);
        fftw_destroy_plan(plan_r2c_);
        fftw_free(data_);
        fftw_free(dataF_);
    }

    void applyDFT(const FReal* sampledData,
                  FComplex<FReal>* transformedData) const
    {
        FBlas::c_setzero(nsteps_opt_,reinterpret_cast<FReal*>(dataF_));
        FBlas::copy(nsteps_, sampledData,data_);

        // perform fft
        fftw_execute( plan_r2c_ );

        // write transformed data
        FBlas::c_copy(nsteps_opt_,reinterpret_cast<FReal*>(dataF_),
                      reinterpret_cast<FReal*>(transformedData));
    }


    void applyIDFT(const FComplex<FReal>* transformedData,
                   FReal* sampledData) const
    {
        // read transformed data
        FBlas::setzero(nsteps_,data_);
        FBlas::c_copy(nsteps_opt_,reinterpret_cast<const FReal*>(transformedData),
                      reinterpret_cast<FReal*>(dataF_));

        // perform ifft
        fftw_execute( plan_c2r_ );

        for( int s=0; s<dim*nsteps_; ++s)
            data_[s]/=nsteps_; // the fft from fftw is not scaled !!!!!!!!

        // write sampled data
        FBlas::copy(nsteps_,data_,sampledData);

    }


};


/**
 * @class FCFft
 *
 * @tparam nsteps number of sampled values \f$N\f$
 */
template<class FReal, int DIM = 1>
class FCFft
{
    enum{dim = DIM};

    // arrays
    FComplex<FReal>* data_;  //< data in physical space
    FComplex<FReal>* dataF_; //< data in Fourier space

    // plans
    fftw_plan plan_b_; // backward FFT plan
    fftw_plan plan_f_; // forward FFT plan

private:
    int steps_[dim]; //< number of steps per dimension
    int nsteps_;     //< total number of steps
public:

    FCFft(const int steps[dim])
    {
        // init number of steps
        nsteps_=1;
        for(int d = 0; d<dim;++d){
            steps_[d]=steps[d];
            nsteps_*=steps[d];
        }

        // allocate arrays
        data_  = (FComplex<FReal>*) fftw_malloc(sizeof(FComplex<FReal>) * nsteps_);
        dataF_ = (FComplex<FReal>*) fftw_malloc(sizeof(FComplex<FReal>) * nsteps_);

        // multidim fftw plans
        plan_b_ =
                fftw_plan_dft(dim, steps_,
                              reinterpret_cast<fftw_complex*>(dataF_),
                              reinterpret_cast<fftw_complex*>(data_),
                              FFTW_BACKWARD,
                              FFTW_MEASURE);// TODO: test FFTW_ESTIMATE
        plan_f_ =
                fftw_plan_dft(dim, steps_,
                              reinterpret_cast<fftw_complex*>(data_),
                              reinterpret_cast<fftw_complex*>(dataF_),
                              FFTW_FORWARD,
                              FFTW_MEASURE);

    }

    virtual ~FCFft()
    {
        fftw_destroy_plan(plan_b_);
        fftw_destroy_plan(plan_f_);
        fftw_free(data_);
        fftw_free(dataF_);
    }

    void applyDFT(const FComplex<FReal>* sampledData,
                  FComplex<FReal>* transformedData) const
    {
        FBlas::c_setzero(nsteps_,reinterpret_cast<FReal*>(dataF_));
        FBlas::c_copy(nsteps_, reinterpret_cast<const FReal*>(sampledData),
                      reinterpret_cast<FReal*>(data_));

        // perform fft
        fftw_execute( plan_f_ );

        // write transformed data
        FBlas::c_copy(nsteps_,reinterpret_cast<FReal*>(dataF_),
                      reinterpret_cast<FReal*>(transformedData));
    }


    void applyIDFT(const FComplex<FReal>* transformedData,
                   FComplex<FReal>* sampledData) const
    {
        // read transformed data
        FBlas::c_setzero(nsteps_,reinterpret_cast<FReal*>(data_));
        FBlas::c_copy(nsteps_,reinterpret_cast<const FReal*>(transformedData),
                      reinterpret_cast<FReal*>(dataF_));

        // perform ifft
        fftw_execute( plan_b_ );

        for(unsigned int s=0; s<nsteps_; ++s)
            data_[s]*=1./(nsteps_); // the fft from fftw is not scaled !!!!!!!!

        // write sampled data
        FBlas::c_copy(nsteps_,reinterpret_cast<FReal*>(data_),
                      reinterpret_cast<FReal*>(sampledData));

    }


};
#endif /*ScalFMM_USE_FFT*/
#endif /* FDFT_HPP */

