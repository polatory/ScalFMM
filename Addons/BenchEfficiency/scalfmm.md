ScalFMM with StarPU+CUDA
========================

## Installing the libraries

For some installation steps, we provide a "valid-if" command which provide a basic test to ensure it should work.
In case of success `STEP-OK` will be print-out.
In addition, if a library is already installed on the system, it is possible to set the output variables directly and test with the "valid-if" command if it should work.

The installation and configuration to have the execution traces and executions times are marked as __Optional__ but higly recommended since they let have the efficiencies.
However, if one wants to execute without any overhead, it might need to remove the usage of FXT.

### Pre-requiste:
In order to follow this tutorial, it is needed to have the following applications installed:

* autoconf (>= 2.69)
* gawk (Awk >= 4.0.1)
* make (>= 3.81) 
* cmake (>= 3.2.2)
* gcc/g++ (>= 4.9) and the gcc/g++ names should point to the correct binaries
* BLAS/LAPACK (The configure of ScalFMM is different if the MKL is used or not, but with the MKL it is recommended to set environment variable `MKLROOT`)
* CUDA (>= 7) and `CUDA_PATH` must be set. In our case, `CUDA_PATH=/usr/local/cuda-7.5/`
* __Optional__ Vite (from `sudo apt-get install vite` or see [http://vite.gforge.inria.fr/download.php](http://vite.gforge.inria.fr/download.php))

> Some installations of CUDA does not have libcuda file.
> In this case, one needs to create a link : `sudo ln /usr/local/cuda-7.5/lib64/libcudart.so /usr/local/cuda-7.5/lib64/libcuda.so`

> [Plafrim-Developers] 
>
> Alloc a node : salloc -N 1 --time=03:00:00 --exclusive -p court_sirocco -CHaswell --gres=gpu:4 -x sirocco06
> 
> Find it: squeue and ssh on it
>
> Load modules : module load compiler/gcc/4.9.2 cuda75/toolkit/7.5.18 intel/mkl/64/11.2/2016.0.0 build/cmake/3.2.1

### Working directory

The variable `SCALFMM_TEST_DIR` is used to specify the working directory:
```bash
export SCALFMM_TEST_DIR=~/scalfmm_test
mkdir $SCALFMM_TEST_DIR
cd $SCALFMM_TEST_DIR
```

*Output variables:* `$SCALFMM_TEST_DIR`

Valid-if
```bash
if [[ -n $SCALFMM_TEST_DIR ]] && [[ -d $SCALFMM_TEST_DIR ]] ; then
   echo “STEP-OK”
fi
```

### Downloading Packages

In case the node used for compiling/testing do not have access to internet, then download the following packages first.

```bash
cd $SCALFMM_TEST_DIR
wget https://www.open-mpi.org/software/hwloc/v1.11/downloads/hwloc-1.11.2.tar.gz
wget http://download.savannah.gnu.org/releases/fkt/fxt-0.2.11.tar.gz # Optional
wget http://www.fftw.org/fftw-3.3.4.tar.gz
svn export svn://scm.gforge.inria.fr/svnroot/starpu/trunk starpu
git clone --depth=1 https://scm.gforge.inria.fr/anonscm/git/scalfmm-public/scalfmm-public.git
```

### HWLOC
```bash
cd $SCALFMM_TEST_DIR
if [[ ! -f hwloc-1.11.2.tar.gz ]] ; then
    wget https://www.open-mpi.org/software/hwloc/v1.11/downloads/hwloc-1.11.2.tar.gz
fi
tar xvf hwloc-1.11.2.tar.gz
cd hwloc-1.11.2/
SCALFMM_HWLOC_DIR=$SCALFMM_TEST_DIR/hwlocinstall
./configure --prefix=$SCALFMM_HWLOC_DIR
make install
```

*Output variables:* `$SCALFMM_HWLOC_DIR`

Valid-if:
```bash
if [[ -n $SCALFMM_HWLOC_DIR ]] && [[ -d $SCALFMM_HWLOC_DIR/lib/ ]] && [[ -f  $SCALFMM_HWLOC_DIR/lib/libhwloc.so ]]; then
   echo “OK”
fi
```

### FXT (__Optional__)
```bash
cd $SCALFMM_TEST_DIR
if [[ ! -f fxt-0.2.11.tar.gz ]] ; then
    wget http://download.savannah.gnu.org/releases/fkt/fxt-0.2.11.tar.gz
fi
tar xvf fxt-0.2.11.tar.gz
cd fxt-0.2.11/
SCALFMM_FXT_DIR=$SCALFMM_TEST_DIR/fxtinstall
./configure --prefix=$SCALFMM_FXT_DIR
make install
```

*Output variables:* `$SCALFMM_FXT_DIR`

Valid-if:
```bash
if [[ -n $SCALFMM_FXT_DIR ]] && [[ -d $SCALFMM_FXT_DIR/lib/ ]] && [[ -f  $SCALFMM_FXT_DIR/lib/libfxt.so ]]; then
   echo “OK”
fi
```

### FFTW
The MKL can be used otherwise we need the FFTW lib:
```bash
cd $SCALFMM_TEST_DIR
if [[ ! -f fftw-3.3.4.tar.gz ]] ; then
    wget http://www.fftw.org/fftw-3.3.4.tar.gz
fi    
tar xvf fftw-3.3.4.tar.gz
cd fftw-3.3.4/
SCALFMM_FFTW_DIR=$SCALFMM_TEST_DIR/fftinstall
./configure --prefix=$SCALFMM_FFTW_DIR
make install
./configure --prefix=$SCALFMM_FFTW_DIR --enable-float
make install
```

*Output variables:* `$SCALFMM_FFTW_DIR`

Valid-if:
```bash
if [[ -n $SCALFMM_FFTW_DIR ]] && [[ -d $SCALFMM_FFTW_DIR/lib/ ]] && [[ -f  $SCALFMM_FFTW_DIR/lib/libfftw3.a ]] && [[ -f  $SCALFMM_FFTW_DIR/lib/libfftw3f.a ]]; then
   echo “OK”
fi
```

### StarPU
```bash
cd $SCALFMM_TEST_DIR
if [[ ! -d starpu ]] ; then
	svn export svn://scm.gforge.inria.fr/svnroot/starpu/trunk starpu
fi    
cd starpu/
SCALFMM_STARPU_DIR=$SCALFMM_TEST_DIR/starpuinstall
./autogen.sh
./configure --prefix=$SCALFMM_STARPU_DIR --with-fxt=$SCALFMM_FXT_DIR --with-hwloc=$SCALFMM_HWLOC_DIR --with-cuda-dir=$CUDA_PATH --disable-opencl
make install
```
> __Optional__ In case you do not want to use trace (FXT) please remove the `--with-fxt=$SCALFMM_FXT_DIR` parameter from the command

*Output variables:* `$SCALFMM_STARPU_DIR`

Valid-if:
```bash
if [[ -n $SCALFMM_STARPU_DIR ]] && [[ -d $SCALFMM_STARPU_DIR/lib/ ]] && [[ -f  $SCALFMM_STARPU_DIR/lib/libstarpu.so ]] ; then
   echo “OK”
fi
```

### ScalFMM

#### Configure
+ Getting the source from the last commit:
```bash
cd $SCALFMM_TEST_DIR
if [[ ! -d scalfmm-public ]] ; then
    git clone --depth=1 https://scm.gforge.inria.fr/anonscm/git/scalfmm-public/scalfmm-public.git
fi    
cd scalfmm-public/Build/
export SCALFMM_BUILD_DIR=`pwd`
```

*Output variables:* `SCALFMM_BUILD_DIR`

+ Configure (No MKL):
```bash
cmake .. -DSCALFMM_BUILD_DEBUG=OFF -DSCALFMM_USE_MPI=OFF -DSCALFMM_BUILD_TESTS=ON -DSCALFMM_BUILD_UTESTS=OFF -DSCALFMM_USE_BLAS=ON -DSCALFMM_USE_MKL_AS_BLAS=OFF -DSCALFMM_USE_LOG=ON -DSCALFMM_USE_STARPU=ON -DSCALFMM_USE_CUDA=ON -DSCALFMM_USE_OPENCL=OFF -DHWLOC_DIR=$SCALFMM_HWLOC_DIR -DSTARPU_DIR=$SCALFMM_STARPU_DIR -DSCALFMM_USE_FFT=ON -DFFT_DIR=$SCALFMM_FFT_DIR
```
+ Configure (MKL BLAS/LAPACK and FFTW):
```bash
cmake .. -DSCALFMM_BUILD_DEBUG=OFF -DSCALFMM_USE_MPI=OFF -DSCALFMM_BUILD_TESTS=ON -DSCALFMM_BUILD_UTESTS=OFF -DSCALFMM_USE_BLAS=ON -DSCALFMM_USE_MKL_AS_BLAS=ON -DSCALFMM_USE_LOG=ON -DSCALFMM_USE_STARPU=ON -DSCALFMM_USE_CUDA=ON -DSCALFMM_USE_OPENCL=OFF -DHWLOC_DIR=$SCALFMM_HWLOC_DIR -DSTARPU_DIR=$SCALFMM_STARPU_DIR -DSCALFMM_USE_FFT=ON -DFFT_DIR=$SCALFMM_FFT_DIR
```
+ Configure (MKL BLAS/LAPACK/FFT and No FFTW):

> [Plafrim-Developers] Should use that one

```bash
cmake .. -DSCALFMM_BUILD_DEBUG=OFF -DSCALFMM_USE_MPI=OFF -DSCALFMM_BUILD_TESTS=ON -DSCALFMM_BUILD_UTESTS=OFF -DSCALFMM_USE_BLAS=ON -DSCALFMM_USE_MKL_AS_BLAS=ON -DSCALFMM_USE_LOG=ON -DSCALFMM_USE_STARPU=ON -DSCALFMM_USE_CUDA=ON -DSCALFMM_USE_OPENCL=OFF -DHWLOC_DIR=$SCALFMM_HWLOC_DIR -DSTARPU_DIR=$SCALFMM_STARPU_DIR -DSCALFMM_USE_FFT=ON -DSCALFMM_USE_MKL_AS_FFTW=ON
```

Valid-if:
```
cmake .. ; if [[ "$?" == "0" ]] ; then echo "OK" ; fi
```

#### Build

```bash
cd $SCALFMM_BUILD_DIR
make testBlockedUnifCudaBench
```

Valid-if:
```
ls ./Tests/Release/testBlockedUnifCudaBench ; if [[ "$?" == "0" ]] ; then echo "OK" ; fi
```

#### Basic Executions

Information for scalfmm binaries

* Passing `--help` as parameter provide the possible/valid parameters
* Simulation properties are choosen by :
  * `-h` : height of the tree
  * `-bs` : granularity/size of the group
  * `-nb` : number of particles generated
* Execution properties are choosen by the StarPU environment variables :
  * `STARPU_NCPUS` : the number of CPU workers
  * `STARPU_NCUDA` : the number of GPU workers (for heterogeneous binary)
* By default the application will not compare the FMM interactions against the direct method (which is N^2) and so it is recommended to avoid the validation for large test cases. But to get the accuracy one must pass the parameter `-validation`
* `-p2p-m2l-cuda-only` : to compute the P2P and the M2L only on GPU (the rest on the CPU)

Examples:

```
STARPU_NCPUS=3 STARPU_NCUDA=1 ./Tests/Release/testBlockedUnifCudaBench -nb 10000 -h 3
```

+ Visualize the execution trace (__Optional__)

Convert the fxt file
```bash
$SCALFMM_STARPU_DIR/bin/starpu_fxt_tool -i "/tmp/prof_file_"$USER"_0"
```
Then visualize the output with vite
```bash
vite ./paje.trace
```

+ Get execution times

```bash
python $SCALFMM_STARPU_DIR/bin/starpu_trace_state_stats.py -t trace.rec
```

Should give something like:
```
"Name","Count","Type","Duration"
"Initializing",3,"Runtime",5.027746
"Overhead",37,"Runtime",0.110073
"Idle",13,"Other",0.03678
"Scheduling",24,"Runtime",16.529527
"Sleeping",17,"Other",2197.255516
"FetchingInput",10,"Runtime",0.012637
"execute_on_all_wrapper",6,"Task",8.431909
"PushingOutput",10,"Runtime",16.505568
"P2P",1,"Task",105.131112
"Callback",4,"Runtime",0.001048
"Deinitializing",3,"Runtime",0.014547
"P2M",1,"Task",2.543303
"L2P",1,"Task",5.649106
"M2L-level-2",1,"Task",2.167273
```

## Generating Execution Results

For test case `-nb 10000000` (10 million) and `-h 6` (height of the tree equal to 6),
we first want to know the best granularity `-bs`.

This parameter will certainly not be the same for sequential/parallel/heterogenous configurations.

```bash
SCALFMM_NB=10000000
SCALFMM_H=7
SCALFMM_MIN_BS=100
SCALFMM_MAX_BS=3000
SCALFMM_MAX_NB_CPU=24
SCALFMM_MAX_NB_GPU=4
```

```bash
SCALFMM_BS_CPU_SEQ=`scalfmm_bench_get_best_bs -nb $SCALFMM_NB -h $SCALFMM_H -start $SCALFMM_MIN_BS -end $SCALFMM_MAX_BS | scalfmm_extract_key "@BEST BS" `

SCALFMM_BS_CPU_PAR=`scalfmm_bench_get_best_bs -nb $SCALFMM_NB -h $SCALFMM_H -start $SCALFMM_MIN_BS -end $SCALFMM_MAX_BS | scalfmm_extract_key "@BEST BS" `

SCALFMM_BS_CPU_GPU=`scalfmm_bench_get_best_bs -nb $SCALFMM_NB -h $SCALFMM_H -start $SCALFMM_MIN_BS -end $SCALFMM_MAX_BS | scalfmm_extract_key "@BEST BS" `
```

Then, we can execute three best configurations, and keep .rec for each of them:
```bash
STARPU_NCPUS=1
STARPU_NCUDA=0
./Tests/Release/testBlockedUnifCudaBench -nb $SCALFMM_NB -h $SCALFMM_H -bs $SCALFMM_CPU_SEQ
SCALFMM_SEQ_REC="trace-nb_$SCALFMM_NB-h_$SCALFMM_H-bs_$SCALFMM_CPU_SEQ-CPU_$STARPU_NCPUS-GPU_$STARPU_NCUDA.rec"
mv trace.rec $SCALFMM_SEQ_REC

STARPU_NCPUS=$SCALFMM_MAX_NB_CPU
STARPU_NCUDA=0
./Tests/Release/testBlockedUnifCudaBench -nb $SCALFMM_NB -h $SCALFMM_H -bs $SCALFMM_BS_CPU_PAR
SCALFMM_PAR_REC="trace-nb_$SCALFMM_NB-h_$SCALFMM_H-bs_$SCALFMM_CPU_SEQ-CPU_$STARPU_NCPUS-GPU_$STARPU_NCUDA.rec"
mv trace.rec $SCALFMM_PAR_REC

STARPU_NCPUS=$SCALFMM_MAX_NB_CPU
STARPU_NCUDA=$SCALFMM_MAX_NB_GPU
./Tests/Release/testBlockedUnifCudaBench -nb $SCALFMM_NB -h $SCALFMM_H -bs $SCALFMM_BS_CPU_GPU
SCALFMM_PAR_CPU_GPU_REC="trace-nb_$SCALFMM_NB-h_$SCALFMM_H-bs_$SCALFMM_CPU_SEQ-CPU_$STARPU_NCPUS-GPU_$STARPU_NCUDA.rec"
mv trace.rec $SCALFMM_PAR_CPU_GPU_REC
```

And we also want the GPU tasks only on GPU
```bash
STARPU_NCPUS=$SCALFMM_MAX_NB_CPU
STARPU_NCUDA=$SCALFMM_MAX_NB_GPU
./Tests/Release/testBlockedUnifCudaBench -nb $SCALFMM_NB -h $SCALFMM_H -bs $SCALFMM_BS_CPU_GPU -p2p-m2l-cuda-only
SCALFMM_PAR_GPU_REC="trace-nb_$SCALFMM_NB-h_$SCALFMM_H-bs_$SCALFMM_CPU_SEQ-CPU_$STARPU_NCPUS-GPU_$STARPU_NCUDA-GPUONLY.rec"
mv trace.rec $SCALFMM_PAR_GPU_REC
```

And we want the sequential version with parallel granularity:
```bash
STARPU_NCPUS=1
STARPU_NCUDA=0

./Tests/Release/testBlockedUnifCudaBench -nb $SCALFMM_NB -h $SCALFMM_H -bs $SCALFMM_BS_CPU_PAR
SCALFMM_SEQ_CPU_BS_REC="trace-nb_$SCALFMM_NB-h_$SCALFMM_H-bs_$SCALFMM_CPU_SEQ-CPU_$STARPU_NCPUS-GPU_$STARPU_NCUDA.rec"
mv trace.rec SCALFMM_SEQ_CPU_BS_REC

./Tests/Release/testBlockedUnifCudaBench -nb $SCALFMM_NB -h $SCALFMM_H -bs $SCALFMM_BS_CPU_GPU
SCALFMM_SEQ_GPU_BS_REC="trace-nb_$SCALFMM_NB-h_$SCALFMM_H-bs_$SCALFMM_CPU_SEQ-CPU_$STARPU_NCPUS-GPU_$STARPU_NCUDA.rec"
mv trace.rec $SCALFMM_SEQ_GPU_BS_REC
```

From these files, we are able to get the different efficencies.

## Post-processing and Plot

Getting all the efficency
Solving the linear programming problem

Plotting the results


## Automatization

```bash
SCALFMM_NB=10000000
SCALFMM_H=7
SCALFMM_MIN_BS=100
SCALFMM_MAX_BS=3000
SCALFMM_MAX_NB_CPU=24
SCALFMM_MAX_NB_GPU=4

scalfmm_generate_efficiency -nb $SCALFMM_NB -h $SCALFMM_H -start $SCALFMM_MIN_BS -end $SCALFMM_MAX_BS
```