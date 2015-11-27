
// @SCALFMM_PRIVATE

#include "../Src/Containers/FDiv2Bissection.hpp"
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

    const int dim = FParameters::getValue(argc, argv, DimParam.options, 100);
    const int height = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 4);
    const char* outputdir = FParameters::getStr(argc, argv, SvgOutParam.options, "/tmp/");

    std::cout << "Config : dim = " << dim << "\n";
    std::cout << "Config : height = " << height << "\n";
    std::cout << "Config : outputdir = " << outputdir << "\n";

    typedef double FReal;
    typedef FDenseBlock<FReal> LeafClass;
    typedef FDenseBlock<FReal> CellClass;
    typedef FDiv2Bissection<FReal, LeafClass, CellClass> GridClass;

    {
        GridClass bissection(dim, height);

        FSvgRect output(outputdir, "classic", dim);

        bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
            output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
        });
    }

    {
        const int nbPartitions = FMath::pow2(height-1);
        std::unique_ptr<int[]> partitions(new int[nbPartitions]);
        {
            int nbValuesLeft = dim;
            for(int idxPartition = 0 ; idxPartition < nbPartitions-1 ; ++idxPartition){
                partitions[idxPartition] = FMath::Max(1, int(drand48()*(nbValuesLeft-(nbPartitions-idxPartition))));
                nbValuesLeft -= partitions[idxPartition];
            }
            partitions[nbPartitions-1] = nbValuesLeft;
        }

        GridClass bissection(dim, height, partitions.get(), nbPartitions);

        FSvgRect output(outputdir, "partitioned", dim);

        bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
            output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
        });
    }

    return 0;
}
