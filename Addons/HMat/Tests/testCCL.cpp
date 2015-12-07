
// @SCALFMM_PRIVATE

#include "../Src/Clustering/FCCLTreeCluster.hpp"
#include "../Src/Clustering/FClusterTree.hpp"
#include "../Src/Utils/FMatrixIO.hpp"

#include "../Src/Containers/FStaticDiagonalBisection.hpp"
#include "../Src/Utils/FSvgRect.hpp"
#include "../Src/Viewers/FDenseBlockWrapper.hpp"
#include "../Src/Blocks/FDenseBlock.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#include <memory>

int main(int argc, char** argv){
    static const FParameterNames SvgOutParam = {
        {"-fout", "--out", "-out"} ,
         "Svg output directory."
    };
    static const FParameterNames DimParam = {
        {"-N", "-nb", "-dim"} ,
         "Dim of the matrix."
    };

    FHelpDescribeAndExit(argc, argv,"Test the bisection.",SvgOutParam,DimParam,FParameterDefinitions::OctreeHeight);

    const char* filename = FParameters::getStr(argc, argv, FParameterDefinitions::InputFile.options, "../Addons/HMat/Data/unitCube1000.bin");
    const int height = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 4);
    const char* outputdir = FParameters::getStr(argc, argv, SvgOutParam.options, "/tmp/");

    int readNbRows = 0;
    int readNbCols = 0;
    // Read distances
    double* distances = nullptr;
    FAssertLF(FMatrixIO::read(filename, &distances, &readNbRows, &readNbCols));
    FAssertLF(readNbRows == readNbCols);
    const int dim = readNbRows;

    FCCLTreeCluster<double> partitioner(dim, distances, CCL::CCL_TM_MAXIMUM /*CCL::CCL_TM_AVG_LINKAGE*/);

    FClusterTree<double> tclusters;
    partitioner.fillClusterTree(&tclusters);
    tclusters.checkData();

    std::unique_ptr<int[]> permutations(new int[dim]);
    tclusters.fillPermutations(permutations.get());

    const int nbPartitions = FMath::pow2(height-1);
    std::unique_ptr<int[]> partitions(new int[nbPartitions]);
    tclusters.getPartitions(height, nbPartitions, partitions.get());

    {
        typedef double FReal;
        typedef FDenseBlock<FReal> LeafClass;
        typedef FDenseBlock<FReal> CellClass;
        typedef FStaticDiagonalBisection<FReal, LeafClass, CellClass> GridClass;

        GridClass bissection(dim, height, partitions.get(), nbPartitions);

        FSvgRect output(outputdir, "ccl.svg", dim);

        bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
            output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
        });
    }

    tclusters.saveToXml(outputdir, "ccl.xml");

    tclusters.saveToDot(outputdir, "gccl.dot");

    return 0;
}

