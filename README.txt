
ScalFmm, Inria, Please read the licence.

To compile:
---------------------------------------------------
# Go to 
cd scalfmm/Build
# Use cmake first by
cmake ..
# Or if you want to use MPI
cmake .. -DScalFMM_USE_MPI=ON

# Configure, for example with:
ccmake ..
# turn on/off the options you want
# For example, it is advised to build the tests and have a look in scalfmm/Tests/
# Finally you can build by
make
# And access executables in scalfmm/Build/Tests/{Release,Debug}/.....

Build the doc
---------------------------------------------------
In scalfmm/Doc you can find several pdf and .tex file about
the implementation, kernels and data structure.

# The developer documentation is generated using DOxygen: 
cd scalfmm/Build
# Tape 
make doc
# This will create a Html dir
browser scalfmm/Build/Doc/html/index.html


