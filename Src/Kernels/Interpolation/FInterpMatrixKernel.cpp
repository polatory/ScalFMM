#include "FInterpMatrixKernel.hpp"


/// ID_OVER_R
const unsigned int FInterpMatrixKernel_IOR::indexTab[]={0,0,0,1,1,2,
                                                        0,1,2,1,2,2};
                   
const unsigned int FInterpMatrixKernel_IOR::applyTab[]={0,1,2,
                                                        1,3,4,
                                                        2,4,5};

/// R_IJ
const unsigned int FInterpMatrixKernel_R_IJ::indexTab[]={0,0,0,1,1,2,
                                                         0,1,2,1,2,2};
                   
const unsigned int FInterpMatrixKernel_R_IJ::applyTab[]={0,1,2,
                                                         1,3,4,
                                                         2,4,5};

/// R_IJK
const unsigned int FInterpMatrixKernel_R_IJK::indexTab[]={0,0,0,1,1,1,2,2,2,0,
                                                          0,1,2,0,1,2,0,1,2,1,
                                                          0,1,2,0,1,2,0,1,2,2};
                   
const unsigned int FInterpMatrixKernel_R_IJK::applyTab[]={0,3,6,3,1,9,6,9,2,
                                                          3,1,9,1,4,7,9,7,5,
                                                          6,9,2,9,7,5,2,5,8};
