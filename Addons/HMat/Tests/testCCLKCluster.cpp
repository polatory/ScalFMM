
// @SCALFMM_PRIVATE

// HMAT specific includes
#include "../Src/Clustering/FCCLKCluster.hpp"
#include "../Src/Utils/FMatrixIO.hpp"

#include "../Src/Containers/FPartitionsMapping.hpp"
#include "../Src/Utils/FSvgRect.hpp"
#include "../Src/Viewers/FDenseBlockWrapper.hpp"
#include "../Src/Blocks/FDenseBlock.hpp"

// ScalFMM's include
#include "Files/FFmaGenericLoader.hpp" // load grid
#include "Utils/FGlobal.hpp"
#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"
#include "Utils/FParameterNames.hpp"

#include <memory>

int main(int argc, char** argv){
    static const FParameterNames SvgOutParam = {
        {"-fout", "--out", "-out"} ,
         "Svg output directory."
    };
    static const FParameterNames ElemParam = {
        {"-nE", "-nbElements", "-elem"} ,
         "Nbof elements to partition."
    };
    static const FParameterNames DimParam = {
        {"-nD", "-nbDim", "-dim"} ,
         "Dimension of ambiant space."
    };
    static const FParameterNames NbPartitionsParam = {
        {"-Part", "-nbPartitions", "-nbClusters"} ,
         "Dim of the matrix."
    };
    static const FParameterNames NbPassesParam = {
        {"-Pass", "-nbPasses", "-nbPass"} ,
         "Dim of the matrix."
    };

    //FHelpDescribeAndExit(argc, argv,"Test the k-means, k-medians and k-medoids algorithms from CCL (using 2 methods for distance computation) on 3D grid (unitCube).",SvgOutParam,ElemParam,DimParam,NbPartitionsParam,NbPassesParam,FParameterDefinitions::InputFile);

    ////////////////////////////////////////////////////////////////////
    /// Command line examples
    // Partition the unitCube with 1000 particles in 8 partitions (use 5 iterations of EM algorithm)
    // ./Addons/HMat/testCCLKCluster -binin -nE 1000 -d unitCube -Part 8 -Pass 5
    // > Results show that k-medoids seems to work better for unitCube1000 (use fewer pass, but still more costly)
    //
    //
    ////////////////////////////////////////////////////////////////////
    /// Parameters
    // size of the grid
    const int nbElements = FParameters::getValue(argc,argv,"-nE", 1000); 

    // size of the grid
    const int nbDim = FParameters::getValue(argc,argv,"-nD", 3); 

    // Define number of partitions
    const int nbPartitions = FParameters::getValue(argc, argv, NbPartitionsParam.options, 4);

    // Define number of passes
    const int nbPasses = FParameters::getValue(argc, argv, NbPassesParam.options, 5);

    // Precision
    typedef double FReal;
    if(sizeof(FReal)==4)
        std::cout<< "Precision: Single Float" <<std::endl;
    if(sizeof(FReal)==8)
        std::cout<< "Precision: Double Float" <<std::endl;

    ////////////////////////////////////////////////////////////////////
    /// Files and Directories

    // Output directory (SVG vizualisation of partitions)
    const char* outputdir = FParameters::getStr(argc, argv, SvgOutParam.options, "/tmp/");

    // Data path
    const std::string ioPath = FParameters::getStr(argc,argv,"-path", std::string("../Addons/HMat/Data/").c_str());

    // Geometry
    const std::string geometryName(FParameters::getStr(argc,argv,"-d",   "unitCube"));

    // Distribution
    std::ostringstream oss_nbElements; oss_nbElements << nbElements;
    const std::string distributionName(geometryName + oss_nbElements.str());

    // Compute partitions from distance matrix
    {

        const char* filename = (ioPath + distributionName + ".bin").c_str();

        std::cout << "Distance matrix read from file:" << filename << std::endl;

        int readNbRows = 0;
        int readNbCols = 0;
        // Read distances
        FReal* distances = nullptr;
        FAssertLF(FMatrixIO::read(filename, &distances, &readNbRows, &readNbCols));
        FAssertLF(readNbRows == readNbCols);
        FAssertLF(readNbRows == nbElements);


        // Test Cluster Center Method K-MEDOIDS
        {
            FTic timeKMedoids;
            double tKMedoids;
            timeKMedoids.tic();

            FCCLKCluster<FReal> partitioner(nbPartitions, nbElements, distances, nbPasses);

            tKMedoids = timeKMedoids.tacAndElapsed();
            std::cout << "... took @tKMedoids = "<< tKMedoids <<"\n";

            std::unique_ptr<int[]> partitions(new int[nbPartitions]);
            partitioner.getPartitions(nbPartitions,partitions.get());
            {
                typedef FDenseBlock<FReal> CellClass;
                typedef FPartitionsMapping<FReal, CellClass> GridClass;

                GridClass bissection(nbElements, partitions.get(), nbPartitions);

                char svgName[1024];
                sprintf(svgName, "%s/%s-%s-P%d.svg", outputdir, distributionName.c_str(), "CCL_KMEDOIDS_DIST_MEAN", nbPartitions);
                FSvgRect output(svgName, nbElements);
                std::cout << "\tSave svg to " << svgName << "\n";

                bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
                    output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
                });
            }

        }

    }


    // Compute partitions from data matrix (i.e. grid)
    {

        // Read geometry
        std::string distributionFileName = ioPath + distributionName;
        if(  FParameters::existParameter(argc, argv, FParameterDefinitions::InputBinFormat.options)){
            distributionFileName += ".bfma";
        }
        else {
            distributionFileName += ".fma";
        }

        std::cout << "Data matrix read from file:" << distributionFileName << std::endl;

        // open particle file
        FFmaGenericLoader<FReal> loader(distributionFileName.c_str());
        if(!loader.isOpen()) throw std::runtime_error("Particle distribution file couldn't be opened!");

        ////////////////////////////////////////////////////////////////////
        /// Load grid from distribution file
        //FPoint<FReal>* grid = new FPoint<FReal>[nbElements];
        FReal* grid = new FReal[nbElements*nbDim];

        for(FSize idxPart = 0 ; idxPart < nbElements ; ++idxPart){
            FPoint<FReal> position;
            FReal physicalValue = 0.0;
            loader.fillParticle(&position,&physicalValue);
            //grid[idxPart]=position;
            grid[idxPart+0*nbElements]=position.getX();
            grid[idxPart+1*nbElements]=position.getY();
            grid[idxPart+2*nbElements]=position.getZ();
        }


        // Test Cluster Center Method K-MEANS
        {
            FTic timeKMeans;
            double tKMeans;
            timeKMeans.tic();

            FCCLKCluster<FReal> partitioner(nbPartitions, nbElements, nbDim, grid, CCL::CCL_CCM_ARITHMETIC_MEAN, CCL::CCL_DIST_MEAN, nbPasses);

            tKMeans = timeKMeans.tacAndElapsed();
            std::cout << "... took @tKMeans = "<< tKMeans <<"\n";


            std::unique_ptr<int[]> partitions(new int[nbPartitions]);
            partitioner.getPartitions(nbPartitions,partitions.get());
            {
                typedef FDenseBlock<FReal> CellClass;
                typedef FPartitionsMapping<FReal, CellClass> GridClass;

                GridClass bissection(nbElements, partitions.get(), nbPartitions);

                char svgName[1024];
                sprintf(svgName, "%s/%s-%s-P%d.svg", outputdir, distributionName.c_str(), "CCL_KMEANS_DIST_MEAN", nbPartitions);
                FSvgRect output(svgName, nbElements);
                std::cout << "\tSave svg to " << svgName << "\n";

                bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
                    output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
                });
            }

        }

        // Test Cluster Center Method K-MEDIANS
        {
            FTic timeKMedians;
            double tKMedians;
            timeKMedians.tic();

            FCCLKCluster<FReal> partitioner(nbPartitions, nbElements, nbDim, grid, CCL::CCL_CCM_MEDIAN, CCL::CCL_DIST_MEAN, nbPasses);

            tKMedians = timeKMedians.tacAndElapsed();
            std::cout << "... took @tKMedians = "<< tKMedians <<"\n";

            std::unique_ptr<int[]> partitions(new int[nbPartitions]);
            partitioner.getPartitions(nbPartitions,partitions.get());
            {
                typedef FDenseBlock<FReal> CellClass;
                typedef FPartitionsMapping<FReal, CellClass> GridClass;

                GridClass bissection(nbElements, partitions.get(), nbPartitions);

                char svgName[1024];
                sprintf(svgName, "%s/%s-%s-P%d.svg", outputdir, distributionName.c_str(), "CCL_KMEDIANS_DIST_MEAN", nbPartitions);
                FSvgRect output(svgName, nbElements);
                std::cout << "\tSave svg to " << svgName << "\n";

                bissection.forAllBlocksDescriptor([&](const FBlockDescriptor& info){
                    output.addRectWithLegend(info.col, info.row, info.nbCols, info.nbRows, info.level);
                });
            }
        }

    }



    return 0;
}

