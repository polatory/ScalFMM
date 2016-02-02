// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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

// @SCALFMM_PRIVATE

#include "../Src/Containers/FStaticDiagonalBisection.hpp"
#include "../Src/Viewers/FMatDensePerm.hpp"
#include "../Src/Blocks/FDenseBlock.hpp"
#include "../Src/Blocks/FSVDBlock.hpp"
#include "../Src/Blocks/FACABlock.hpp"

#include "../Src/Clustering/FMaxDistCut.hpp"
#include "../Src/Clustering/FCCLTreeCluster.hpp"
#include "../Src/Clustering/FGraphThreshold.hpp"

#include "../Src/Utils/FSvgRect.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#include "Utils/FTic.hpp"

#include <memory>

template <class FReal, class MatrixClass>
void CheckRank(const FClusterTree<FReal>& ctree, MatrixClass& matrix,
               const char outputdir[], const char configName[], const int height){

    typedef FDenseBlock<FReal> LeafClass;
    typedef FACABlock<FReal,7> CellClass;
    typedef FStaticDiagonalBisection<FReal, LeafClass, CellClass> GridClass;

    const int dim = matrix.getDim();
    std::unique_ptr<int[]> permutationOrgToNew(new int[dim]);

    for( int idxLevel = 1 ; idxLevel < height ; ++idxLevel){
        const int nbPartitions = FMath::pow2(idxLevel);
        std::unique_ptr<int[]> partitions(new int[nbPartitions]);
        ctree.getPartitions(idxLevel+1, nbPartitions, partitions.get());

        ctree.fillPermutations(permutationOrgToNew.get());
        matrix.setPermutOrigToNew(permutationOrgToNew.get());

        std::cout << "\tLevel " << idxLevel << " build blocks\n";
        FTic timer;
        GridClass grid(dim, idxLevel+1, partitions.get(), nbPartitions);
        grid.fillBlocks(matrix);
        std::cout << "\tdone in " << timer.tacAndElapsed() << "s\n";

        {
            char svgName[1024];
            sprintf(svgName, "%s/%s-%d.svg", outputdir, configName, idxLevel);
            std::cout << "\tSave svg to " << svgName << "\n";

            FSvgRect output(svgName, dim);

            grid.forAllCellBlocks([&](const FBlockDescriptor& info, const CellClass& cell){
                output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level, cell.getRank());
            });

            grid.forAllLeafBlocks([&](const FBlockDescriptor& info, const LeafClass& /*leaf*/){
                output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
            });
        }
    }
}

int main(int argc, char** argv){

    static const FParameterNames SvgOutParam = {
        {"-fout", "--out", "-out"} ,
         "Svg output directory."
    };

    FHelpDescribeAndExit(argc, argv, "Test the rank for different clustering.",
                         FParameterDefinitions::InputFileOne,
                         FParameterDefinitions::InputFileTwow,
                         FParameterDefinitions::OctreeHeight,
                         SvgOutParam);

    ////////////////////////////////////////////////////////////////////

    const char* outputdir = FParameters::getStr(argc, argv, SvgOutParam.options, "/tmp/");
    //const char* distanceFilename = FParameters::getStr(argc, argv, FParameterDefinitions::InputFileOne.options, "../Addons/HMat/Data/unitCube1000.bin");
    //const char* matrixFilename = FParameters::getStr(argc, argv, FParameterDefinitions::InputFileTwow.options, "../Addons/HMat/Data/unitCube1000_ONE_OVER_R.bin");
    const char* distanceFilename = FParameters::getStr(argc, argv, FParameterDefinitions::InputFileOne.options, "../Addons/HMat/Data/unitSphere1000.bin");
    const char* matrixFilename = FParameters::getStr(argc, argv, FParameterDefinitions::InputFileTwow.options, "../Addons/HMat/Data/unitSphere1000_GAUSS100.bin");
    typedef double FReal;

    //const char* distanceFilename = FParameters::getStr(argc, argv, FParameterDefinitions::InputFileOne.options, "../Addons/HMat/Data/first_reads_1k-fmr-dissw.bin");
    //const char* matrixFilename = FParameters::getStr(argc, argv, FParameterDefinitions::InputFileTwow.options, "../Addons/HMat/Data/first_reads_1k-fmr-covar.bin");
    //typedef float FReal;

    const int height = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 4);

    std::cout << "Check until height = " << height << "\n";

    typedef FMatDensePerm<FReal> MatrixClass;

    ////////////////////////////////////////////////////////////////////

    std::cout << "Load distances file " << distanceFilename << "\n";

    std::unique_ptr<FReal[]> distanceValues;
    int distanceValuesDim = 0;
    {
        int readNbRows = 0;
        int readNbCols = 0;
        FReal* distanceValuesPtr = nullptr;
        FAssertLF(FMatrixIO::read(distanceFilename, &distanceValuesPtr, &readNbRows, &readNbCols));
        FAssertLF(readNbRows == readNbCols);
        distanceValuesDim = readNbRows;
        distanceValues.reset(distanceValuesPtr);
    }

    ////////////////////////////////////////////////////////////////////

    std::cout << "Load matrix file " << matrixFilename << "\n";

    MatrixClass matrix(matrixFilename);
    const int matrixDim = matrix.getDim();

    // Display covariance matrix
    const FSize displaySize = 10;
        std::cout<<"\nC=["<<std::endl;
        for ( int i=0; i<displaySize; ++i) {
            for ( int j=0; j<displaySize; ++j)
                std::cout << matrix.getVal(i,j) << " ";
            std::cout<< std::endl;
        }
        std::cout<<"]"<<std::endl;


    FAssertLF(distanceValuesDim == matrixDim);
    std::cout << "Matrices dim = " << matrixDim << "\n";

    ////////////////////////////////////////////////////////////////////

    {
        std::cout << "Test FMaxDistCut\n";
        FMaxDistCut<FReal> partitioner(matrixDim, distanceValues.get());

        FClusterTree<FReal> tclusters;
        partitioner.fillClusterTree(&tclusters);
        tclusters.checkData();
        CheckRank<FReal, MatrixClass>(tclusters, matrix, outputdir, "FMaxDistCut", height);
    }
    ////////////////////////////////////////////////////////////////////
    {
        std::cout << "Test FGraphThreshold\n";
        FGraphThreshold<FReal> partitioner(matrixDim, distanceValues.get(), FGraphThreshold<FReal>::GetDefaultRadius(matrixDim, distanceValues.get()));

        FClusterTree<FReal> tclusters;
        partitioner.fillClusterTree(&tclusters);
        tclusters.checkData();
        CheckRank<FReal, MatrixClass>(tclusters, matrix, outputdir, "FGraphThreshold", height);
    }
    ////////////////////////////////////////////////////////////////////
    {
        std::cout << "Test FCCLTreeCluster\n";
        FCCLTreeCluster<double> partitioner(matrixDim, distanceValues.get(), CCL::CCL_TM_MAXIMUM);

        FClusterTree<double> tclusters;
        partitioner.fillClusterTree(&tclusters);
        tclusters.checkData();
        CheckRank<FReal, MatrixClass>(tclusters, matrix, outputdir, "FCCLTreeCluster", height);
    }

    return 0;
}


