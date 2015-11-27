
// @SCALFMM_PRIVATE

#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <random> // for normal distribution generator
#include <algorithm> // for sort


// ScalFMM includes
#include "Files/FFmaGenericLoader.hpp"
#include "Utils/FGlobal.hpp"
#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FMemUtils.hpp"
#include "Utils/FBlas.hpp" // for FBlas::potrf (and QR,SVD...)
#include "Utils/FMath.hpp"
#include "Utils/FParameterNames.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp" // for kernel matrices

// not mandatory but useful to define some flags
#include "Core/FFmmAlgorithm.hpp"
#include "Core/FFmmAlgorithmThread.hpp"

// ScalFMM/HMat includes
#include "../Src/Utils/FMatrixIO.hpp"

/**
* In this file we explicitely build a kernel matrix and then store it in a binary file. 
*
* Author: Pierre Blanchard (pierre.blanchard@inria.fr)
* Date created: November 26th, 2015 
*/


int main(int argc, char* argv[])
{
    std::cout << "Addons/HMat: Build a kernel matrix and then store it in a binary file." << std::endl;

    // Example for unitCube5000
    // ./Addons/HMat/testGenerateKernelMatrix -binin -N 5000 -d unitCube5000

    ////////////////////////////////////////////////////////////////////
    /// Timers 
    FTic time; 
    time.tic();

    ////////////////////////////////////////////////////////////////////
    /// Parameters
    // Verbose (print matrices)
    const int verbose = FParameters::getValue(argc, argv, "-v", 0);
    std::cout<< "Verbose level: " << verbose <<std::endl;
    // size of the grid
    const FSize matrixSize = FParameters::getValue(argc,argv,"-N", 1000); 
    std::cout<< "Full matrix size: " << matrixSize <<std::endl;
    // Precision
    typedef double FReal;
    if(sizeof(FReal)==4)
        std::cout<< "Precision: Single Float" <<std::endl;
    if(sizeof(FReal)==8)
        std::cout<< "Precision: Double Float" <<std::endl;

    // Memory
    std::cout<< "Memory requirement: " << sizeof(FReal) * matrixSize * matrixSize *1e-6 << " MBytes" <<std::endl;

    // Data path
    const std::string ioPath = FParameters::getStr(argc,argv,"-path", std::string("../Addons/HMat/Data/").c_str());

    // Read geometry
    const std::string distributionName(FParameters::getStr(argc,argv,"-d",   "unitCube1000"));
    std::string distributionFileName = ioPath + distributionName;
    if(  FParameters::existParameter(argc, argv, FParameterDefinitions::InputBinFormat.options)){
        distributionFileName += ".bfma";
    }
    else {
        distributionFileName += ".fma";
    }

    // open particle file
    FFmaGenericLoader<FReal> loader(distributionFileName.c_str());
    if(!loader.isOpen()) throw std::runtime_error("Particle distribution file couldn't be opened!");

    ////////////////////////////////////////////////////////////////////
    /// Load 3D grid from distribution file
    FPoint<FReal>* grid = new FPoint<FReal>[matrixSize];
    FReal* pV = new FReal[matrixSize];

    for(FSize idxPart = 0 ; idxPart < matrixSize ; ++idxPart){
        FPoint<FReal> position;
        FReal physicalValue = 0.0;
        loader.fillParticle(&position,&physicalValue);
        grid[idxPart]=position;
        pV[idxPart]=physicalValue;
    }



    ////////////////////////////////////////////////////////////////////
    /// Build kernel matrix K
    
    // Interaction kernel evaluator
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
    const MatrixKernelClass MatrixKernel;
    const std::string MatrixKernelID = MatrixKernelClass::getID();//.c_str();

    // Allocate memory
    FReal* K = new FReal[matrixSize*matrixSize];

    // Build kernel matrix
    FTic timeAssK;
    for(FSize idxRow = 0 ; idxRow < matrixSize  ; ++idxRow)
        for(FSize idxCol = 0 ; idxCol < matrixSize  ; ++idxCol)
            K[idxRow*matrixSize+idxCol] = MatrixKernel.evaluate(grid[idxRow],
                                                                grid[idxCol]);
    double tAssK = timeAssK.tacAndElapsed();
    std::cout << "... took @tAssK = "<< tAssK <<"\n";

    ////////////////////////////////////////////////////////////////////
    /// Write kernel matrix in binary file

    // Output file name
    const std::string matrixName((distributionName + "_" + MatrixKernelID).c_str()); // if nothing specified then use file associated with n=50
    const std::string fileName = ioPath + matrixName + ".bin";

    // Write 
    std::cout<< "Write matrix in binary file: " << fileName << "\n";
    FTic timeWriteMat;
    double tWriteMat;
    timeWriteMat.tic();

    FMatrixIO::write<FReal>(matrixSize,matrixSize,K,fileName);

    tWriteMat = timeWriteMat.tacAndElapsed();
    std::cout << "... took @tWriteMat = "<< tWriteMat <<"\n";

    /// Free memory
    delete[] K;

    return 0;
}
