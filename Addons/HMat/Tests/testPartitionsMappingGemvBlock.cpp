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

#include "../Src/Containers/FBlockPMapping.hpp"
#include "../Src/Viewers/FMatDensePerm.hpp"
#include "../Src/Blocks/FDenseBlock.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#include <memory>

int main(int argc, char** argv){
    static const FParameterNames DimParam = {
        {"-N", "-nb", "-dim"} ,
         "Dim of the matrix."
    };
    static const FParameterNames PartitionsParam = {
        {"-part", "-parts", "-nbparts"} ,
         "Number of partitions."
    };

    FHelpDescribeAndExit(argc, argv, "Test the bisection.", DimParam, FParameterDefinitions::OctreeHeight);

    const int dim = FParameters::getValue(argc, argv, DimParam.options, 100);
    const int nbPartitions = FParameters::getValue(argc, argv, PartitionsParam.options, 5);

    std::cout << "Config : dim = " << dim << "\n";
    std::cout << "Config : nbPartitions = " << nbPartitions << "\n";


    typedef double FReal;
    typedef FMatDensePerm<FReal> MatrixClass;

    MatrixClass matrix(dim);
    for(int idxRow = 0; idxRow < dim ; ++idxRow){
        for(int idxCol = 0; idxCol < dim ; ++idxCol){
            matrix.setVal(idxRow, idxCol, 1./(FMath::Abs(FReal(idxRow-idxCol))+1.));
        }
    }

    { // Here we fill the block manually
        // We consider a fack permutation
        std::unique_ptr<int[]> permutations(new int[dim]);
        for(int idx = 0 ; idx < dim ; ++idx){
            permutations[idx] = idx;
        }
        // Set permutation to matrix
        matrix.setPermutOrigToNew(permutations.get());
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
        typedef FDenseBlock<FReal> CellClass;
        typedef FBlockPMapping<FReal, CellClass> GridClass;

        std::unique_ptr<int[]> partitions(new int[nbPartitions]);
        {
            int nbValuesLeft = dim;
            for(int idxPartition = 0 ; idxPartition < nbPartitions-1 ; ++idxPartition){
                partitions[idxPartition] = FMath::Max(1, int(drand48()*(nbValuesLeft-(nbPartitions-idxPartition))));
                nbValuesLeft -= partitions[idxPartition];
            }
            partitions[nbPartitions-1] = nbValuesLeft;
        }

        GridClass grid(dim, partitions.get(), nbPartitions);

        // We iterate on the blocks
        // V blocks cover all the rows, but only some columns (based on the clustering)
        for(int idxColBlock = 0 ; idxColBlock < nbPartitions ; ++idxColBlock){
            const MatrixClass::BlockDescriptor colBlock = matrix.getBlock(grid.getVBlockInfo(idxColBlock));
            int rj = -1;
            /// TODO HERE
            /// Compute rj, and the resulting Vj blocks,
            /// use the colBlock (some or all of its values)

            /// TODO END
            // Store the result in grid.getVBlock(idxColBlock)
            CellClass& Vj = grid.getVBlock(idxColBlock);
            Vj.resize(rj, colBlock.getNbCols());
            for(int idxRow = 0 ; idxRow < Vj.getNbRows() ; ++idxRow){
                for(int idxCol = 0 ; idxCol < Vj.getNbCols() ; ++idxCol){
                    /// TODO HERE
                    /// Fill Vj with the result
                    Vj.setValue(idxRow, idxCol, -1);
                    /// TODO END
                }
            }
        }
        // U blocks cover all the columns, but only some rows (based on the clustering)
        for(int idxRowBlock = 0 ; idxRowBlock < nbPartitions ; ++idxRowBlock){
            const MatrixClass::BlockDescriptor rowBlock = matrix.getBlock(grid.getUBlockInfo(idxRowBlock));
            int ri = -1;
            /// TODO HERE
            /// Compute ri, and the resulting Ui blocks
            /// use the rowBlock (some or all of its values)

            /// TODO END
            // Store the result in grid.getUBlock(idxRowBlock)
            CellClass& Ui = grid.getUBlock(idxRowBlock);
            Ui.resize(rowBlock.getNbRows(), ri);
            for(int idxRow = 0 ; idxRow < Ui.getNbRows() ; ++idxRow){
                for(int idxCol = 0 ; idxCol < Ui.getNbCols() ; ++idxCol){
                    /// TODO HERE
                    /// Fill Vj with the result
                    Ui.setValue(idxRow, idxCol, -1);
                    /// TODO END
                }
            }
        }
        // Build the core part
        for(int idxColBlock = 0 ; idxColBlock < nbPartitions ; ++idxColBlock){
            for(int idxRowBlock = 0 ; idxRowBlock < nbPartitions ; ++idxRowBlock){
                const CellClass& Ui = grid.getUBlock(idxRowBlock);
                const CellClass& Vj = grid.getVBlock(idxColBlock);
                // Store the result in grid.getCBlock(idxRowBlock, idxColBlock)
                CellClass& Cij = grid.getCBlock(idxRowBlock, idxColBlock);
                Cij.resize(Vj.getNbRows(), Ui.getNbCols());
                for(int idxRow = 0 ; idxRow < Cij.getNbRows() ; ++idxRow){
                    for(int idxCol = 0 ; idxCol < Cij.getNbCols() ; ++idxCol){
                        /// TODO HERE
                        /// Fill Vj with the result
                        Cij.setValue(idxRow, idxCol, -1);
                        /// TODO END
                    }
                }
            }
        }

        std::unique_ptr<FReal[]> resDense(new FReal[dim]);
        FSetToZeros(resDense.get(), dim);

        //grid.gemv(resDense.get(), vec.get());

        FMath::FAccurater<FReal> testDense(resTest.get(), resDense.get(), dim);

        std::cout << "Test Dense partitions mapping, Error = " << testDense << "\n";
    }

    return 0;
}



