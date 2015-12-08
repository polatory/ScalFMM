
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

template <class FReal, class PartitionerClass>
void TestPartitions(const PartitionerClass& partitioner, const int height,
                    const char outputdir[], const int dim,
                    const FReal coords[], const char config[]){
    FClusterTree<double> tclusters;
    partitioner.fillClusterTree(&tclusters);
    tclusters.checkData();

    {
        char filenameBuffer[1024];
        sprintf(filenameBuffer, "%s/%s.xml", outputdir, config);
        std::cout << "\t save xml to " << filenameBuffer << "\n";
        tclusters.saveToXml(filenameBuffer);
        sprintf(filenameBuffer, "%s/%s.dot", outputdir, config);
        std::cout << "\t save cdot to " << filenameBuffer << "\n";
        tclusters.saveToDot(filenameBuffer);
    }

    std::unique_ptr<int[]> permutations(new int[dim]);
    std::unique_ptr<int[]> invpermutations(new int[dim]);
    tclusters.fillPermutations(permutations.get(), invpermutations.get());

    for(int idxLevel = 2 ; idxLevel <= height ; ++idxLevel){
        const int nbPartitions = FMath::pow2(idxLevel-1);
        std::unique_ptr<int[]> partitions(new int[nbPartitions]);
        tclusters.getPartitions(idxLevel, nbPartitions, partitions.get());

        {
            char coordfilename[1024];
            sprintf(coordfilename, "%s/%s-coord-%d.csv", outputdir, config, idxLevel);
            std::cout << "\t save coord for level " << idxLevel << " to " << coordfilename << "\n";
            FILE* fcoord = fopen(coordfilename, "w");

            int offsetParticles = 0;
            for(int idxPartition = 0 ; idxPartition < nbPartitions ; ++idxPartition){
                for(int idxPart = 0 ; idxPart < partitions[idxPartition] ; ++idxPart){
                    const int idxUnk = invpermutations[idxPart+offsetParticles];
                    fprintf(fcoord, "%e,%e,%e,%d\n",
                            coords[idxUnk*4 + 0],
                            coords[idxUnk*4 + 1],
                            coords[idxUnk*4 + 2],
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
            sprintf(svgfilename, "%s/%s-%d.svg", outputdir, config, idxLevel);
            std::cout << "\t save svg for level " << idxLevel << " to " << svgfilename << "\n";
            FSvgRect output(svgfilename, dim);

            bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
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
    static const FParameterNames DimParam = {
        {"-N", "-nb", "-dim"} ,
         "Dim of the matrix."
    };

    FHelpDescribeAndExit(argc, argv,"Test the bisection.",SvgOutParam,DimParam,FParameterDefinitions::OctreeHeight,
                         FParameterDefinitions::NbParticles);

    const int nbParticles = (FParameters::getValue(argc, argv, FParameterDefinitions::NbParticles.options, 3000) & ~1);
    const int height = FParameters::getValue(argc, argv, FParameterDefinitions::OctreeHeight.options, 6);
    const char* outputdir = FParameters::getStr(argc, argv, SvgOutParam.options, "/tmp/");

    typedef double FReal;

    /////////////////////////////////////////////////////////////////

    std::cout << "Create spheres\n";

    const FReal radius = 0.5;
    const FReal distanceBetweenSpheres = 0.4;
    const int nbPointsSphere = nbParticles/2;

    std::cout << "\tradius = " << radius << "\n";
    std::cout << "\tdistanceBetweenSpheres = " << distanceBetweenSpheres << "\n";
    std::cout << "\tnbPointsSphere = " << nbPointsSphere << "\n";

    std::unique_ptr<FReal[]> spherePoints(new FReal[nbParticles*4]);
    unifRandonPointsOnSphere(nbPointsSphere, radius, spherePoints.get()) ;

    for(int idxPoint = 0 ; idxPoint < nbPointsSphere ; ++idxPoint){
        spherePoints[(idxPoint+nbPointsSphere)*4 + 0] = spherePoints[idxPoint*4 + 0] + radius*2 + distanceBetweenSpheres;
        spherePoints[(idxPoint+nbPointsSphere)*4 + 1] = spherePoints[idxPoint*4 + 1] + radius*2 + distanceBetweenSpheres;
        spherePoints[(idxPoint+nbPointsSphere)*4 + 2] = spherePoints[idxPoint*4 + 2] + radius*2 + distanceBetweenSpheres;
    }

    /////////////////////////////////////////////////////////////////

    std::cout << "Compute distance\n";

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

    /////////////////////////////////////////////////////////////////
    {
        std::cout << "Test FGraphThreshold\n";
        FGraphThreshold<FReal> partitioner(dim, distances.get(), radius/4);
        TestPartitions<FReal, FGraphThreshold<FReal>>(partitioner,height, outputdir,
                                                      dim, spherePoints.get(), "FGraphThreshold");
    }
    /////////////////////////////////////////////////////////////////
    {
        std::cout << "Test FCCLTreeCluster\n";
        FCCLTreeCluster<FReal> partitioner(dim, distances.get(), CCL::CCL_TM_MAXIMUM);
        TestPartitions<FReal, FCCLTreeCluster<FReal>>(partitioner,height, outputdir,
                                                      dim, spherePoints.get(), "FCCLTreeCluster");
    }
    /////////////////////////////////////////////////////////////////
    {
        std::cout << "Test FMaxDistCut\n";
        FMaxDistCut<FReal> partitioner(dim, distances.get());
        TestPartitions<FReal, FMaxDistCut<FReal>>(partitioner,height, outputdir,
                                                      dim, spherePoints.get(), "FMaxDistCut");
    }

    return 0;
}



