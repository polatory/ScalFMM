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
#include "../Src/Viewers/FMatDense.hpp"
#include "../Src/Blocks/FDenseBlock.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#include "Utils/FTic.hpp"

#include <memory>

int main(int argc, char** argv){
    static const FParameterNames DimParam = {
        {"-N", "-nb", "-dim"} ,
         "Dim of the matrix."
    };

    FHelpDescribeAndExit(argc, argv, "Test the bisection.", DimParam, FParameterDefinitions::OctreeHeight);

    ////////////////////////////////////////////////////////////////////
    /// Timers 
    FTic time; 
    time.tic();

    const int dim = FParameters::getValue(argc, argv, DimParam.options, 100);
    const int height = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 4);

    std::cout << "Config : dim = " << dim << "\n";
    std::cout << "Config : height = " << height << "\n";


    typedef double FReal;
    typedef FMatDense<FReal> MatrixClass;

    MatrixClass matrix(dim);
    for(int idxRow = 0; idxRow < dim ; ++idxRow){
        for(int idxCol = 0; idxCol < dim ; ++idxCol){
            matrix.setVal(idxRow, idxCol, 1./(FMath::Abs(FReal(idxRow-idxCol))+1.));
        }
    }

    std::unique_ptr<FReal[]> vec(new FReal[dim]);
    for(int idxVal = 0 ; idxVal < dim ; ++idxVal){
        vec[idxVal] = 1.0;
    }

    std::unique_ptr<FReal[]> resTest(new FReal[dim]);
    FSetToZeros(resTest.get(), dim);
    {
        for(int idxRow = 0; idxRow < dim ; ++idxRow){
            for(int idxCol = 0; idxCol < dim ; ++idxCol){
                resTest[idxRow] += vec[idxCol] * matrix.getVal(idxRow, idxCol);
            }
        }
    }

    {
        typedef FDenseBlock<FReal> LeafClass;
        typedef FDenseBlock<FReal> CellClass;
        typedef FStaticDiagonalBisection<FReal, LeafClass, CellClass> GridClass;

        GridClass grid(dim, height);
        grid.fillBlocks(matrix);

        std::unique_ptr<FReal[]> resDense(new FReal[dim]);
        FSetToZeros(resDense.get(), dim);

        grid.gemv(resDense.get(), vec.get());

        FMath::FAccurater<FReal> testDense(resTest.get(), resDense.get(), dim);

        std::cout << "Test Dense, Error = " << testDense << "\n";
    }

    {
        typedef FDenseBlock<FReal> LeafClass;
        typedef FDenseBlock<FReal> CellClass;
        typedef FStaticDiagonalBisection<FReal, LeafClass, CellClass> GridClass;

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

        GridClass grid(dim, height, partitions.get(), nbPartitions);
        grid.fillBlocks(matrix);

        std::unique_ptr<FReal[]> resDense(new FReal[dim]);
        FSetToZeros(resDense.get(), dim);

        grid.gemv(resDense.get(), vec.get());

        FMath::FAccurater<FReal> testDense(resTest.get(), resDense.get(), dim);

        std::cout << "Test Dense with partitions, Error = " << testDense << "\n";
    }

    double tOverall = time.tacAndElapsed();
    std::cout << "... took @tOverall = "<< tOverall <<"\n";

    return 0;
}
