
// @SCALFMM_PRIVATE

#include "../Src/Clustering/FCCLKCluster.hpp"
#include "../Src/Utils/FMatrixIO.hpp"

#include "../Src/Containers/FPartitionsMapping.hpp"
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
    static const FParameterNames NbPartitionsParam = {
        {"-P", "-nbPartitions", "-nbClusters"} ,
         "Dim of the matrix."
    };


    FHelpDescribeAndExit(argc, argv,"Test the k-medoid algorithm from CCL (using 2 methods for cluster center computation).",SvgOutParam,DimParam,NbPartitionsParam,FParameterDefinitions::InputFile);

    const char* filename = FParameters::getStr(argc, argv, FParameterDefinitions::InputFile.options, "../Addons/HMat/Data/unitCube1000.bin");

    typedef float FReal;

    const char* outputdir = FParameters::getStr(argc, argv, SvgOutParam.options, "/tmp/");


    int readNbRows = 0;
    int readNbCols = 0;
    // Read distances
    FReal* distances = nullptr;
    FAssertLF(FMatrixIO::read(filename, &distances, &readNbRows, &readNbCols));
    FAssertLF(readNbRows == readNbCols);
    const int dim = readNbRows;

    // Define number of partitions
    const int nbPartitions = FParameters::getValue(argc, argv, NbPartitionsParam.options, 4);
   
    // Test Cluster Center Method ARITHMETIC_MEAN
    {
        FCCLKCluster<FReal> partitioner(nbPartitions, dim, distances, CCL::CCL_CCM_ARITHMETIC_MEAN);

        std::unique_ptr<int[]> partitions(new int[nbPartitions]);
        partitioner.getPartitions(nbPartitions,partitions.get());
        {
            typedef FDenseBlock<FReal> CellClass;
            typedef FPartitionsMapping<FReal, CellClass> GridClass;

            GridClass bissection(dim, partitions.get(), nbPartitions);

            char svgName[1024];
            sprintf(svgName, "%s/%s-%d.svg", outputdir, "CCL_CCM_ARITHMETIC_MEAN", nbPartitions);
            FSvgRect output(svgName, dim);
            std::cout << "\tSave svg to " << svgName << "\n";

            bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
                output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
            });
        }

    }

    // Test Cluster Center Method MEDIAN
    {
        FCCLKCluster<FReal> partitioner(nbPartitions, dim, distances, CCL::CCL_CCM_MEDIAN);

        std::unique_ptr<int[]> partitions(new int[nbPartitions]);
        partitioner.getPartitions(nbPartitions,partitions.get());
        {
            typedef FDenseBlock<FReal> CellClass;
            typedef FPartitionsMapping<FReal, CellClass> GridClass;

            GridClass bissection(dim, partitions.get(), nbPartitions);

            char svgName[1024];
            sprintf(svgName, "%s/%s-%d.svg", outputdir, "CCL_CCM_MEDIAN", nbPartitions);
            FSvgRect output(svgName, dim);
            std::cout << "\tSave svg to " << svgName << "\n";

            bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
                output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
            });
        }
    }





    return 0;
}

