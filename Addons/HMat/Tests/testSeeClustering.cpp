
// @SCALFMM_PRIVATE

#include "../Src/Clustering/FCCLTreeCluster.hpp"
#include "../Src/Clustering/FGraphThreshold.hpp"
#include "../Src/Clustering/FCCLTreeCluster.hpp"
#include "../Src/Clustering/FMaxDistCut.hpp"
#include "../Src/Utils/FMatrixIO.hpp"

#include "../Src/Containers/FStaticDiagonalBisection.hpp"
#include "../Src/Utils/FSvgRect.hpp"
#include "../Src/Viewers/FDenseBlockWrapper.hpp"
#include "../Src/Blocks/FDenseBlock.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#include "Files/FGenerateDistribution.hpp"

#include <memory>

// FUSE_SCOTCH


int main(int argc, char** argv){
    static const FParameterNames SvgOutParam = {
        {"-fout", "--out", "-out"} ,
         "Svg output directory."
    };
    static const FParameterNames DimParam = {
        {"-N", "-nb", "-dim"} ,
         "Dim of the matrix."
    };

    FHelpDescribeAndExit(argc, argv,"Test the bisection.",SvgOutParam,DimParam,FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::NbParticles);

    const int nbParticles = (FParameters::getValue(argc, argv, FParameterDefinitions::NbParticles.options, 3000) & ~1);
    const int height = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 4);
    const char* outputdir = FParameters::getStr(argc, argv, SvgOutParam.options, "/tmp/");

    typedef double FReal;

    const FReal radius = 0.5;
    const FReal distanceBetweenSpheres = 0.4;
    const int nbPointsSphere = nbParticles/2;

    std::unique_ptr<FReal[]> spherePoints(new FReal[nbParticles*4]);
    unifRandonPointsOnSphere(nbPointsSphere, radius, spherePoints.get()) ;

    for(int idxPoint = 0 ; idxPoint < nbPointsSphere ; ++idxPoint){
        spherePoints[(idxPoint+nbPointsSphere)*4 + 0] = spherePoints[idxPoint*4 + 0] + radius*2 + distanceBetweenSpheres;
        spherePoints[(idxPoint+nbPointsSphere)*4 + 1] = spherePoints[idxPoint*4 + 1] + radius*2 + distanceBetweenSpheres;
        spherePoints[(idxPoint+nbPointsSphere)*4 + 2] = spherePoints[idxPoint*4 + 2] + radius*2 + distanceBetweenSpheres;
    }

    const int dim = nbParticles;
    std::unique_ptr<FReal[]> distances(new FReal[nbParticles*nbParticles]);
    for(int idxCol = 0 ; idxCol < nbParticles ; ++idxCol){
        for(int idxRow = 0 ; idxRow < nbParticles ; ++idxRow){
            const FReal diffx = spherePoints[idxCol*4 + 0] - spherePoints[idxRow*4 + 0];
            const FReal diffy = spherePoints[idxCol*4 + 1] - spherePoints[idxRow*4 + 1];
            const FReal diffz = spherePoints[idxCol*4 + 2] - spherePoints[idxRow*4 + 2];
            distances[idxCol*nbParticles+idxRow] = FMath::Sqrt((diffx*diffx) + (diffy*diffy) + (diffz*diffz));
        }
    }

    //FGraphThreshold<FReal> partitioner(dim, distances.get(), radius/4);
    //FCCLTreeCluster<FReal> partitioner(dim, distances.get(), CCL::CCL_TM_MAXIMUM);
    FMaxDistCut<FReal> partitioner(dim, distances.get());

    FClusterTree<double> tclusters;
    partitioner.fillClusterTree(&tclusters);
    tclusters.checkData();

    tclusters.saveToXml(outputdir, "scotch.xml");
    tclusters.saveToDot(outputdir, "gscotch.dot");

    std::unique_ptr<int[]> permutations(new int[dim]);
    std::unique_ptr<int[]> invpermutations(new int[dim]);
    tclusters.fillPermutations(permutations.get(), invpermutations.get());

    for(int idxLevel = 2 ; idxLevel <= height ; ++idxLevel){
        const int nbPartitions = FMath::pow2(idxLevel-1);
        std::unique_ptr<int[]> partitions(new int[nbPartitions]);
        tclusters.getPartitions(idxLevel, nbPartitions, partitions.get());

        {
            char coordfilename[1024];
            sprintf(coordfilename, "%s/coord-%d.csv", outputdir, idxLevel);
            FILE* fcoord = fopen(coordfilename, "w");

            int offsetParticles = 0;
            for(int idxPartition = 0 ; idxPartition < nbPartitions ; ++idxPartition){
                for(int idxPart = 0 ; idxPart < partitions[idxPartition] ; ++idxPart){
                    const int idxUnk = invpermutations[idxPart+offsetParticles];
                    fprintf(fcoord, "%e,%e,%e,%d\n",
                            spherePoints[idxUnk*4 + 0],
                            spherePoints[idxUnk*4 + 1],
                            spherePoints[idxUnk*4 + 2],
                            idxPartition);
                }
                offsetParticles += partitions[idxPartition];
            }
            FAssertLF(offsetParticles == dim);

            fclose(fcoord);
        }
        {
            typedef FDenseBlock<FReal> LeafClass;
            typedef FDenseBlock<FReal> CellClass;
            typedef FStaticDiagonalBisection<FReal, LeafClass, CellClass> GridClass;

            GridClass bissection(dim, idxLevel, partitions.get(), nbPartitions);

            char svgfilename[1024];
            sprintf(svgfilename, "%s/scotch-%d.svg", outputdir, idxLevel);
            FSvgRect output(svgfilename, dim);

            bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
                output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
            });
        }
    }


    return 0;
}



