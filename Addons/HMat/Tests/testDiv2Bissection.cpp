
// @SCALFMM_PRIVATE

#include "../Src/Core/FDiv2Bissection.hpp"
#include "../Src/Containers/FMatGrid.hpp"
#include "../Src/Utils/FSvgRect.hpp"

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

    FDiv2Bissection<double, int, int> bissection(dim, height);

    {
        FSvgRect output(outputfile, dim);

        bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
            output.addRectWithLegend(info.x, info.y, info.width, info.height, info.level);
        });
    }

//    const double array[] = {1, 2, 3, 4,
//                             1, 2, 3, 4,
//                             1, 2, 3, 4,
//                             1, 2, 3, 4};
//    FMatGrid<double> denseMat(4 , array);

//    FMatGrid<double>::BlockDescriptor blockDescriptor = denseMat.getBlock(1,1,1,1);


    return 0;
}
