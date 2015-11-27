
// @SCALFMM_PRIVATE

#include "../Src/Containers/FDiv2Bissection.hpp"
#include "../Src/Viewers/FMatGrid.hpp"
#include "../Src/Blocks/FDenseMatrix.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

int main(int argc, char** argv){
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

    std::cout << "Config : dim = " << dim << "\n";
    std::cout << "Config : height = " << height << "\n";

    {
        typedef double FReal;
        typedef FMatGrid<FReal> MatrixClass;
        typedef FDenseMatrix<FReal> LeafClass;
        typedef FDenseMatrix<FReal> CellClass;
        typedef FDiv2Bissection<FReal, LeafClass, CellClass> GridClass;

        GridClass bissection(dim, height);
        bissection.fillBlocks(matrix);
    }

    return 0;
}

