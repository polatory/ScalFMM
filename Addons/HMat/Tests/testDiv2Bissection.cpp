
// @SCALFMM_PRIVATE

#include "../Src/Containers/FDiv2Bissection.hpp"
#include "../Src/Utils/FSvgRect.hpp"
#include "../Src/Viewers/FDenseBlockWrapper.hpp"
#include "../Src/Blocks/FDenseMatrix.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

int main(int argc, char** argv){
    static const FParameterNames SvgOutParam = {
        {"-fout", "--out"} ,
         "Svg output filename."
    };
    static const FParameterNames DimParam = {
        {"-N", "-nb", "-dim"} ,
         "Dim of the matrix."
    };
    static const FParameterNames HeightParam = {
        {"-h", "-height"} ,
         "Number of dissection (+1)."
    };

    FHelpDescribeAndExit(argc, argv,
            "Test the bisection.",SvgOutParam,
            DimParam,HeightParam);

    const int dim = FParameters::getValue(argc, argv, DimParam.options, 100);
    const int height = FParameters::getValue(argc, argv, HeightParam.options, 4);
    const char* outputfile = FParameters::getStr(argc, argv, SvgOutParam.options, "/tmp/example.svg");

    std::cout << "Config : dim = " << dim << "\n";
    std::cout << "Config : height = " << height << "\n";
    std::cout << "Config : outputfile = " << outputfile << "\n";

    typedef double FReal;
    typedef FDenseMatrix<FReal> LeafClass;
    typedef FDenseMatrix<FReal> CellClass;
    typedef FDiv2Bissection<FReal, LeafClass, CellClass> GridClass;

    GridClass bissection(dim, height);
    {
        FSvgRect output(outputfile, dim);

        bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
            output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
        });
    }

    return 0;
}
