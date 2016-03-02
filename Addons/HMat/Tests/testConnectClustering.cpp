
// @SCALFMM_PRIVATE

#include "../Src/Clustering/FConnexClustering.hpp"
#include "../Src/Utils/FMatrixIO.hpp"

#include "../Src/Containers/FPartitionsMapping.hpp"
#include "../Src/Utils/FSvgRect.hpp"
#include "../Src/Viewers/FDenseBlockWrapper.hpp"
#include "../Src/Blocks/FDenseBlock.hpp"

#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#include "Files/FGenerateDistribution.hpp"

#include <memory>

// FUSE_SCOTCH

template <class FReal, class PartitionerClass>
void TestPartitions(const PartitionerClass& partitioner,
                    const char outputdir[], const int idFile,
                    const int dim,
                    const FReal coords[], const char config[]){

    std::unique_ptr<int[]> permutations(new int[dim]);
    std::unique_ptr<int[]> invpermutations(new int[dim]);
    partitioner.fillPermutations(permutations.get(), invpermutations.get());

    const int nbPartitions = partitioner.getNbPartitions();
    std::unique_ptr<int[]> partitions(new int[nbPartitions]);
    partitioner.getPartitions(nbPartitions, partitions.get());

    if(coords!=nullptr)
    {
        char coordfilename[1024];
        sprintf(coordfilename, "%s/%s-coord-%d.csv", outputdir, config, idFile);
        std::cout << "\t save coord to " << coordfilename << "\n";
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
        typedef FDenseBlock<FReal> CellClass;
        typedef FPartitionsMapping<FReal, CellClass> GridClass;

        GridClass bissection(dim, partitions.get(), nbPartitions);

        char svgfilename[1024];
        sprintf(svgfilename, "%s/%s-%d.svg", outputdir, config, idFile);
        std::cout << "\t save svg to " << svgfilename << "\n";
        FSvgRect output(svgfilename, dim);

        bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
            output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
        });
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

    FHelpDescribeAndExit(argc, argv,"Test the bisection.",SvgOutParam,DimParam,
                         FParameterDefinitions::InputFileOne,
                         FParameterDefinitions::InputFileTwow,
                         FParameterDefinitions::NbParticles);

    const int nbParticles = (FParameters::getValue(argc, argv, FParameterDefinitions::NbParticles.options, 3000) & ~1);
    const char* outputdir = FParameters::getStr(argc, argv, SvgOutParam.options, "/tmp/");

    typedef float FReal;

    ////////////////////////////////////////////////////////////////////

    const char* distanceFilename = FParameters::getStr(argc, argv, FParameterDefinitions::InputFileOne.options, "../Addons/HMat/Data/unitCube1000.bin");

    ////////////////////////////////////////////////////////////////////

    std::cout << "Load distances file " << distanceFilename << "\n";

    std::unique_ptr<FReal[]> particles(nullptr);
    std::unique_ptr<FReal[]> distances;
    int dim = 0;
    {
        int readNbRows = 0;
        int readNbCols = 0;
        FReal* distanceValuesPtr = nullptr;
        FAssertLF(FMatrixIO::read(distanceFilename, &distanceValuesPtr, &readNbRows, &readNbCols));
        FAssertLF(readNbRows == readNbCols);
        dim = readNbRows;
        distances.reset(distanceValuesPtr);
    }


    const FSize displaySize=10;
    std::cout<<"\nD=["<<std::endl;
    for ( FSize i=0; i<displaySize; ++i) {
        for ( FSize j=0; j<displaySize; ++j)
            std::cout << distances.get()[i*nbParticles+j] << " ";
        std::cout<< std::endl;
    }
    std::cout<<"]"<<std::endl;

//    /////////////////////////////////////////////////////////////////
//
//    std::cout << "Create spheres\n";
//
//    const FReal radius = 0.5;
//    const FReal distanceBetweenSpheres = 0.4;
//    const int nbPointsSphere = nbParticles/2;
//
//    std::cout << "\tradius = " << radius << "\n";
//    std::cout << "\tdistanceBetweenSpheres = " << distanceBetweenSpheres << "\n";
//    std::cout << "\tnbPointsSphere = " << nbPointsSphere << "\n";
//
//    std::unique_ptr<FReal[]> particles(new FReal[nbParticles*4]);
//    unifRandonPointsOnSphere(nbPointsSphere, radius, particles.get()) ;
//
//    for(int idxPoint = 0 ; idxPoint < nbPointsSphere ; ++idxPoint){
//        particles[(idxPoint+nbPointsSphere)*4 + 0] = particles[idxPoint*4 + 0] + radius*2 + distanceBetweenSpheres;
//        particles[(idxPoint+nbPointsSphere)*4 + 1] = particles[idxPoint*4 + 1] + radius*2 + distanceBetweenSpheres;
//        particles[(idxPoint+nbPointsSphere)*4 + 2] = particles[idxPoint*4 + 2] + radius*2 + distanceBetweenSpheres;
//    }
//
//    /////////////////////////////////////////////////////////////////
//
//    std::cout << "Compute distance\n";
//
//    const int dim = nbParticles;
//    std::unique_ptr<FReal[]> distances(new FReal[nbParticles*nbParticles]);
//    for(int idxCol = 0 ; idxCol < nbParticles ; ++idxCol){
//        for(int idxRow = 0 ; idxRow < nbParticles ; ++idxRow){
//            const FReal diffx = particles[idxCol*4 + 0] - particles[idxRow*4 + 0];
//            const FReal diffy = particles[idxCol*4 + 1] - particles[idxRow*4 + 1];
//            const FReal diffz = particles[idxCol*4 + 2] - particles[idxRow*4 + 2];
//            distances[idxCol*nbParticles+idxRow] = FMath::Sqrt((diffx*diffx) + (diffy*diffy) + (diffz*diffz));
//        }
//    }
//
//    /////////////////////////////////////////////////////////////////

    FReal distMin = std::numeric_limits<FReal>::max();
    FReal distMax = std::numeric_limits<FReal>::min();

    for(int idxRow = 0 ; idxRow < dim ; ++idxRow){
        for(int idxCol = idxRow ; idxCol < dim ; ++idxCol){
            distMin = FMath::Min(distMin, distances[idxCol*nbParticles+idxRow]);
            distMax = FMath::Max(distMax, distances[idxCol*nbParticles+idxRow]);
        }
    }
    std::cout << "Dist min " << distMin << "\n";
    std::cout << "Dist max " << distMax << "\n";

    /////////////////////////////////////////////////////////////////
    const FReal stepThresh = (distMax-distMin)/100;
    for(FReal thresh = distMin ; thresh <= distMax/10 ; thresh += stepThresh){
        std::cout << "Test FConnexClustering for thresh " << thresh << "\n";
        FConnexClustering<FReal> partitioner(dim, distances.get(), thresh);
        std::cout << "\tGot " << partitioner.getNbPartitions() << " partitions\n";
        TestPartitions<FReal, FConnexClustering<FReal>>(partitioner,outputdir,
                                                        int((thresh-distMin)/stepThresh),
                                                        dim, particles.get(), "FConnexClustering");
    }

    return 0;
}




