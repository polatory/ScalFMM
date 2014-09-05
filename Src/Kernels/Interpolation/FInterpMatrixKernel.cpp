#include "FInterpMatrixKernel.hpp"


/// R_IJ
const unsigned int FInterpMatrixKernel_R_IJ::indexTab[]={0,0,0,1,1,2,
                                                         0,1,2,1,2,2};
                   
const unsigned int FInterpMatrixKernel_R_IJ::applyTab[]={0,1,2,
                                                         1,3,4,
                                                         2,4,5};

