
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
#include "Kernels/Interpolation/FInterpMatrixKernel_Covariance.hpp" // for kernel matrices

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


{
    ////////////////////////////////////////////////////////////////////
    /// Build kernel matrix K
    
    // Interaction kernel evaluator
    //typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;
    typedef CK_Gauss<FReal> MatrixKernelClass;
    const FReal lengthScale = FReal(0.01)*FParameters::getValue(argc,argv,"-lengthscale", FReal(100.)); 
    std::ostringstream oss; oss << 100*lengthScale;
    const MatrixKernelClass MatrixKernel(lengthScale);
    const std::string MatrixKernelID = MatrixKernelClass::getID() + oss.str() ;//.c_str();

    // Allocate memory
    FReal* K = new FReal[matrixSize*matrixSize];
    FBlas::setzero(int(matrixSize*matrixSize),K);

    // Build (symmetric) kernel matrix
    FTic timeAssK;
    for(FSize idxRow = 0 ; idxRow < matrixSize  ; ++idxRow)
        for(FSize idxCol = idxRow+1 ; idxCol < matrixSize ; ++idxCol) {
            K[idxRow*matrixSize+idxCol] = MatrixKernel.evaluate(grid[idxRow],
                                                                grid[idxCol]);
            K[idxCol*matrixSize+idxRow] = K[idxRow*matrixSize+idxCol];
        }

    double tAssK = timeAssK.tacAndElapsed();
    std::cout << "... took @tAssK = "<< tAssK <<"\n";

    // Display matrix
    const FSize displaySize = 10;
    if(verbose==2) {
        std::cout<<"\nK=["<<std::endl;
        for ( FSize i=0; i<displaySize; ++i) {
            for ( FSize j=0; j<displaySize; ++j)
                std::cout << K[i*matrixSize+j] << " ";
            std::cout<< std::endl;
        }
        std::cout<<"]"<<std::endl;
    }

    ////////////////////////////////////////////////////////////////////
    /// Write kernel matrix in binary file

    // Output file name
    const std::string matrixName((distributionName + "_" + MatrixKernelID).c_str());
    const std::string fileName = ioPath + matrixName + ".bin";

    // Write 
    std::cout<< "Write matrix in binary file: " << fileName << "\n";
    FTic timeWriteMat;
    double tWriteMat;
    timeWriteMat.tic();

    FMatrixIO::write<FReal>(matrixSize,matrixSize,K,fileName);

    tWriteMat = timeWriteMat.tacAndElapsed();
    std::cout << "... took @tWriteMat = "<< tWriteMat <<"\n";
}
{
    ////////////////////////////////////////////////////////////////////
    /// Build distance matrix D

    // Allocate memory
    FReal* D = new FReal[matrixSize*matrixSize];
    FBlas::setzero(int(matrixSize*matrixSize),D);

    // Build (symmetric) kernel matrix
    FTic timeAssK;

    for(FSize idxRow = 0 ; idxRow < matrixSize  ; ++idxRow)
        for(FSize idxCol = idxRow ; idxCol < matrixSize ; ++idxCol) {
            FReal diffX(grid[idxRow].getX()-grid[idxCol].getX());
            FReal diffY(grid[idxRow].getY()-grid[idxCol].getY());
            FReal diffZ(grid[idxRow].getZ()-grid[idxCol].getZ());
            D[idxRow*matrixSize+idxCol] = FMath::Sqrt(diffX*diffX+diffY*diffY+diffZ*diffZ);
            if(idxCol!=idxRow)
                D[idxCol*matrixSize+idxRow] = D[idxRow*matrixSize+idxCol];
        }

    double tAssK = timeAssK.tacAndElapsed();
    std::cout << "... took @tAssK = "<< tAssK <<"\n";

    // Display matrix
    const FSize displaySize = 10;
    if(verbose==2) {
        std::cout<<"\nD=["<<std::endl;
        for ( FSize i=0; i<displaySize; ++i) {
            for ( FSize j=0; j<displaySize; ++j)
                std::cout << D[i*matrixSize+j] << " ";
            std::cout<< std::endl;
        }
        std::cout<<"]"<<std::endl;
    }

    ////////////////////////////////////////////////////////////////////
    /// Write distance matrix in binary file

    // Output file name
    const std::string fileName = ioPath + distributionName + ".bin";

    // Write 
    std::cout<< "Write matrix in binary file: " << fileName << "\n";
    FTic timeWriteMat;
    double tWriteMat;
    timeWriteMat.tic();

    FMatrixIO::write<FReal>(matrixSize,matrixSize,D,fileName);

    tWriteMat = timeWriteMat.tacAndElapsed();
    std::cout << "... took @tWriteMat = "<< tWriteMat <<"\n";


    double tOverall = time.tacAndElapsed();
    std::cout << "... took @tOverall = "<< tOverall <<"\n";
}

    return 0;
}
