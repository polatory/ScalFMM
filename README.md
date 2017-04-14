# ScalFMM: Fast Multipole Method

----

:warning: ScalFMM has moved to Inria's GitLab: https://gitlab.inria.fr/solverstack/ScalFMM

----

**ScalFMM** is a C++ library that implements a kernel independent Fast Multipole Method.


Copyright Inria, please read the licence.

### Requirements

  - CMake v2.8.12 or later
  - C++ compiler that supports
    - C++14 [compiler support list](http://en.cppreference.com/w/cpp/compiler_support)
    - [OpenMP](http://www.openmp.org/resources/openmp-compilers/)

The following are optional:

  - [Doxygen](http://www.stack.nl/~dimitri/doxygen/) to build the documentation.
  - An MPI implementation to build the distributed files.
  - Custom BLAS, FFT implementations.
  - [StarPU](http://starpu.gforge.inria.fr/) for the relevant FMM implementations.

### Build

``` bash
# Move to the build folder
cd scalfmm/Build
# Use cmake, with relevant options
cmake .. # -DSCALFMM_USE_MPI=ON
```

The build may be configured after the first CMake invocation using, for instance, `ccmake` or `cmake-gui`.

```bash
# Still in the Build folder
ccmake .
# Or
cmake-gui .
```

The binaries are then compiled calling `make`. They can be found in `scalfmm/Build/Tests/{Release,Debug}/...`

An example build using StarPU:

```bash
cmake .. -DSCALFMM_USE_STARPU=ON -DSCALFMM_USE_CUDA=OFF -DSCALFMM_USE_OPENCL=OFF  \
               -DHWLOC_DIR=/home/berenger/Download/hwloc-1.10.0/install/      \
               -DSTARPU_DIR=/home/berenger/Download/starpu-work/StarPU/installwithfxt
```


#### Build the doc:

```bash
cd scalfmm/Build
cmake .. -DSCALFMM_BUILD_DOC=ON # or if cmake has already been called, ccmake .
make doc
```

This will generate the documentation in HTML format in the `Build/Doc/html` folder. You can create a local server to access it using Python

```bash
# From the Build folder
cd Doc/html
python3 -m http.server # or python2 -m SimpleHTTPServer
```

The documentation can then be accessed from an internet browser at the address `localhost:8000`.


### Help and News

You can subscribe to the scalfmm-public-users@lists.gforge.inria.fr mailing list (http://lists.gforge.inria.fr/cgi-bin/mailman/listinfo/scalfmm-public-users). The list is very low trafic (~ 2 mails per year), we will let you know of improvements and releases.

Contact the developers at : scalfmm-public-support@lists.gforge.inria.fr

### Folder structure
  - Src : library core.
  - Data : particle distribution examples.
  - Examples : common usage examples.
  - Doc : documentation configuration.
  - UTests : unit tests.
  - Tests : examples to know how to use scalfmm/put particles in the tree/iterate on the tree...
  - Utils : some scripts and binaries to handle data files.

