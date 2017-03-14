
// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_FFT
// ================
// Keep in private GIT


#include "../../Src/Utils/FGlobal.hpp"

#include "../../Src/GroupTree/Core/FGroupTree.hpp"

#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainer.hpp"

#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "../../Src/Kernels/Uniform/FUnifKernel.hpp"

#include "../../Src/GroupTree/Uniform/FUnifCellPOD.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Files/FRandomLoader.hpp"
#include "../../Src/Files/FFmaGenericLoader.hpp"

#include "../../Src/GroupTree/Core/FGroupSeqAlgorithm.hpp"
#include "../../Src/GroupTree/Core/FGroupTaskAlgorithm.hpp"
#ifdef SCALFMM_USE_OMP4
#include "../../Src/GroupTree/Core/FGroupTaskDepAlgorithm.hpp"
#endif
#ifdef SCALFMM_USE_STARPU
#include "../../Src/GroupTree/Core/FGroupTaskStarpuAlgorithm.hpp"
#include "../../Src/GroupTree/StarPUUtils/FStarPUKernelCapacities.hpp"
#endif
#include "../../Src/GroupTree/Core/FP2PGroupParticleContainer.hpp"

#include "../../Src/Utils/FParameterNames.hpp"

#include <memory>


#define RANDOM_PARTICLES

template <class FReal>
class FSphericalRandomLoader : public FAbstractLoader<FReal> {
protected:
    const int nbParticles;            //< the number of particles
    const FReal boxWidth;             //< the box width
    const FPoint<FReal> centerOfBox;    //< The center of box
    const bool nu;
    const bool snu;
    const bool su;
    const bool elu;
    const bool ssnu;
    const bool elsu;


    FReal rotationMatrix[3][3];

    void initRotationMatrix(){
        const FReal alpha = FMath::FPi<FReal>()/8;
        const FReal omega = FMath::FPi<FReal>()/4;

        FReal yrotation[3][3];
        yrotation[0][0] = FMath::Cos(alpha); yrotation[0][1] = 0.0; yrotation[0][2] = FMath::Sin(alpha);
        yrotation[1][0] = 0.0;               yrotation[1][1] = 1.0; yrotation[1][2] = 0.0;
        yrotation[2][0] = -FMath::Sin(alpha); yrotation[2][1] = 0.0;   yrotation[2][2] = FMath::Cos(alpha);

        FReal zrotation[3][3];
        zrotation[0][0] = FMath::Cos(omega); zrotation[0][1] = -FMath::Sin(omega); zrotation[0][2] = 0.0;
        zrotation[1][0] = FMath::Sin(omega); zrotation[1][1] = FMath::Cos(omega); zrotation[1][2] = 0.0;
        zrotation[2][0] = 0.0; zrotation[2][1] = 0.0;   zrotation[2][2] = 1.0;

        for(int i = 0 ; i < 3 ; ++i){
            for(int j = 0 ; j < 3 ; ++j){
                FReal sum = 0.0;
                for(int k = 0 ; k < 3 ; ++k){
                    sum += zrotation[i][k] * yrotation[k][j];
                }
                rotationMatrix[i][j] = sum;
            }
        }
    }

public:
    /**
    * The constructor need the simulation data
    */
    FSphericalRandomLoader(const int inNbParticles, const bool inNu = false,
                           const bool inSnu = false,
                           const bool inSu = false,
                           const bool inElu = false,
                           const bool inSsnu = false,
                           const bool inElsu = false)
        : nbParticles(inNbParticles), boxWidth(1.0), centerOfBox(0,0,0), nu(inNu),
          snu(inSnu), su(inSu), elu(inElu), ssnu(inSsnu), elsu(inElsu) {
        srand48(static_cast<unsigned int>(0));
        if( !nu && !snu && !su && !elu && !ssnu && !elsu ){
            std::cout << "UNIFORM" << std::endl;
        }
        else if( snu ){
            std::cout << "slightly NON UNIFORM" << std::endl;
        }
        else if( su ){
            std::cout << "SPHERICAL UNIFORM" << std::endl;
        }
        else if( elu ){
            std::cout << "ELLIPSE UNIFORM" << std::endl;
        }
        else if( elsu ){
            std::cout << "ELLIPSE NON UNIFORM" << std::endl;
            initRotationMatrix();
        }
        else if( ssnu ){
            std::cout << "spherical Slightly non UNIFORM" << std::endl;
        }
        else{
            std::cout << "NON UNIFORM" << std::endl;
        }
    }

    /**
    * Default destructor
    */
    virtual ~FSphericalRandomLoader(){
    }

    /**
      * @return true
      */
    bool isOpen() const{
        return true;
    }

    /**
      * To get the number of particles from this loader
      * @param the number of particles the loader can fill
      */
    FSize getNumberOfParticles() const{
        return FSize(this->nbParticles);
    }

    /**
      * The center of the box
      * @return box center
      */
    FPoint<FReal> getCenterOfBox() const{
        return this->centerOfBox;
    }

    /**
      * The box width
      * @return box width
      */
    FReal getBoxWidth() const{
        return this->boxWidth;
    }

    /**
      * Fill a particle
      * @warning to work with the loader, particles has to expose a setPosition method
      * @param the particle to fill
      */
    void fillParticle(FPoint<FReal>* partPtr){
        FPoint<FReal>& inParticle = *partPtr;
        if( !nu && !snu && !su && !elu && !ssnu && !elsu ){
            inParticle.setPosition(
                        (getRandom() * boxWidth) + centerOfBox.getX() - boxWidth/2,
                        (getRandom() * boxWidth) + centerOfBox.getY() - boxWidth/2,
                        (getRandom() * boxWidth) + centerOfBox.getZ() - boxWidth/2);
        }
        else if( snu ){
            const FReal XCenter = centerOfBox.getX();
            const FReal YCenter = centerOfBox.getY();
            const FReal ZCenter = centerOfBox.getZ();

            const FReal rayon = FReal(0.4);
            const FReal thresh = FReal(0.15);
            const FReal threshDiv2 = thresh/2;

            // Generate particles
            const FReal theta = getRandom() * FMath::FPi<FReal>();
            const FReal omega = getRandom() * FMath::FPi<FReal>() * FReal(2);

            const FReal px = rayon * FMath::Cos(omega) * FMath::Sin(theta) + XCenter + thresh * getRandom() - threshDiv2;
            const FReal py = rayon * FMath::Sin(omega) * FMath::Sin(theta) + YCenter + thresh * getRandom() - threshDiv2;
            const FReal pz = rayon * FMath::Cos(theta) + ZCenter + thresh * getRandom() - threshDiv2;

            inParticle.setPosition(px,py,pz);
        }
        else if( su ){
            //http://www.cs.cmu.edu/~mws/rpos.html
            const FReal r = 0.4;

            const FReal pz = getRandom()*2.0*r - r;
            const FReal omega = getRandom() * FMath::FPi<FReal>() * FReal(2);

            const FReal theta = FMath::ASin(pz/r);

            const FReal px = r * cos(theta) * cos(omega);
            const FReal py = r * cos(theta) * sin(omega);

            inParticle.setPosition(px,py,pz);
        }
        else if( elu ){
            const FReal a = 0.4;
            const FReal b = 0.15;

            const FReal maxPerimeter = 2.0 * FMath::FPi<FReal>() * a;

            FReal px   = 0;
            // rayon du cercle pour ce x
            FReal subr = 0;

            do {
                px   = (getRandom() * a * 2) - a;
                subr = FMath::Sqrt( (1.0 - ((px*px)/(a*a))) * (b*b) );
            } while( (getRandom()*maxPerimeter) > subr );

            // on genere un angle
            const FReal omega = getRandom() * FMath::FPi<FReal>() * FReal(2);
            // on recupere py et pz sur le cercle
            const FReal py = FMath::Cos(omega) * subr;
            const FReal pz = FMath::Sin(omega) * subr;

            inParticle.setPosition(px,py,pz);
        }
        else if( elsu ){
            const FReal a = 0.5;
            const FReal b = 0.1;

            const FReal MaxDensity = 10.0;
            const FReal maxPerimeter = 2.0 * FMath::FPi<FReal>() * a ;

            FReal px   = 0;
            // rayon du cercle pour ce x
            FReal subr = 0;
            FReal coef = 1.0;

            do {
                //px   = ( ((getRandom()*8.0+getRandom())/9.0) * a * 2) - a;
                px = (getRandom() * a * 2.0) - a;

                coef = FMath::Abs(px) * MaxDensity/a + 1.0;

                subr = FMath::Sqrt( (1.0 - ((px*px)/(a*a))) * (b*b) );

            } while( (getRandom()*maxPerimeter) > subr * coef );

            // on genere un angle
            const FReal omega = getRandom() * FMath::FPi<FReal>() * FReal(2);
            // on recupere py et pz sur le cercle
            const FReal py = FMath::Cos(omega) * subr;
            const FReal pz = FMath::Sin(omega) * subr;

           // inParticle.setPosition(px,py,pz);
            inParticle.setPosition(px * rotationMatrix[0][0] + py * rotationMatrix[0][1]+ pz * rotationMatrix[0][2],
                                   px * rotationMatrix[1][0] + py * rotationMatrix[1][1]+ pz * rotationMatrix[1][2],
                                   px * rotationMatrix[2][0] + py * rotationMatrix[2][1]+ pz * rotationMatrix[2][2]);

        }
        else if( ssnu ){
            const FReal XCenter = centerOfBox.getX();
            const FReal YCenter = centerOfBox.getY();
            const FReal ZCenter = centerOfBox.getZ();

            const FReal rayon = FReal(0.4);

            // Generate particles
            /*static const int NbAcc = 2;
            FReal acc = 0;
            for(int idx = 0 ; idx < NbAcc ; ++idx){
                acc += getRandom()/FReal(NbAcc);
            }*/
            FReal acc = ((getRandom()*8)+getRandom())/9;

            const FReal theta = acc * FMath::FPi<FReal>();
            const FReal omega = getRandom() * FMath::FPi<FReal>() * FReal(2);

            const FReal px = rayon * FMath::Cos(omega) * FMath::Sin(theta) + XCenter ;
            const FReal py = rayon * FMath::Sin(omega) * FMath::Sin(theta) + YCenter ;
            const FReal pz = rayon * FMath::Cos(theta) + ZCenter ;

            inParticle.setPosition(px,py,pz);
        }
        else{
            const FReal XCenter = centerOfBox.getX();
            const FReal YCenter = centerOfBox.getY();
            const FReal ZCenter = centerOfBox.getZ();

            const FReal rayon = FReal(0.4);

            const FReal theta = getRandom() * FMath::FPi<FReal>();
            const FReal omega = getRandom() * FMath::FPi<FReal>() * FReal(2);

            const FReal px = rayon * FMath::Cos(omega) * FMath::Sin(theta) + XCenter ;
            const FReal py = rayon * FMath::Sin(omega) * FMath::Sin(theta) + YCenter ;
            const FReal pz = rayon * FMath::Cos(theta) + ZCenter ;

            inParticle.setPosition(px,py,pz);
        }
    }

    /** Get a random number between 0 & 1 */
    FReal getRandom() const{
        return FReal(drand48());
    }
};

extern "C" {
int main(int argc, char* argv[]){
    const FParameterNames LocalOptionBlocSize { {"-bs"}, "The size of the block of the blocked tree"};
    const FParameterNames LocalOptionNoValidate { {"-no-validation"}, "To avoid comparing with direct computation"};
    const FParameterNames LocalOptionNonUnif { {"-nonunif"}, "To generate non uniform"};
    const FParameterNames LocalOptionProlateNonUnif { {"-prolate-nonunif"}, "To generate prolate non unif distribution"};
    const FParameterNames LocalOptionProlate { {"-prolate"}, "To generate prolate distribution"};
    const FParameterNames LocalOptionSpherical { {"-spherical"}, "To generate spherical distribution"};
    const FParameterNames LocalOptionSphericalNonUnif { {"-spherical-nonunif"}, "To generate spherical non unif distribution"};
    FHelpDescribeAndExit(argc, argv, "Test the blocked tree by counting the particles.",
                         FParameterDefinitions::OctreeHeight,
#ifdef RANDOM_PARTICLES
                         FParameterDefinitions::NbParticles,
#else
                         FParameterDefinitions::InputFile,
#endif
                         FParameterDefinitions::NbThreads,
                         LocalOptionProlateNonUnif, LocalOptionProlate,
                         LocalOptionSpherical, LocalOptionSphericalNonUnif,LocalOptionNonUnif,
                         LocalOptionBlocSize, LocalOptionNoValidate);

    // Initialize the types
    typedef double FReal;
    static const int ORDER = 5;
    typedef FInterpMatrixKernelR<FReal> MatrixKernelClass;

    typedef FUnifCellPODCore         GroupCellSymbClass;
    typedef FUnifCellPODPole<FReal,ORDER>  GroupCellUpClass;
    typedef FUnifCellPODLocal<FReal,ORDER> GroupCellDownClass;
    typedef FUnifCellPOD<FReal,ORDER>      GroupCellClass;

    typedef FP2PGroupParticleContainer<FReal>          GroupContainerClass;
    typedef FGroupTree< FReal, GroupCellClass, GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupContainerClass, 1, 4, FReal>  GroupOctreeClass;
#ifdef SCALFMM_USE_STARPU
    typedef FStarPUAllCpuCapacities<FUnifKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER>> GroupKernelClass;
    typedef FStarPUCpuWrapper<typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass> GroupCpuWrapper;
    typedef FGroupTaskStarPUAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupCpuWrapper, GroupContainerClass > GroupAlgorithm;
#elif defined(SCALFMM_USE_OMP4)
    typedef FUnifKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER> GroupKernelClass;
    // Set the number of threads
    omp_set_num_threads(FParameters::getValue(argc,argv,FParameterDefinitions::NbThreads.options, omp_get_max_threads()));
    typedef FGroupTaskDepAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass,
            GroupCellSymbClass, GroupCellUpClass, GroupCellDownClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#else
    typedef FUnifKernel<FReal,GroupCellClass,GroupContainerClass,MatrixKernelClass,ORDER> GroupKernelClass;
    //typedef FGroupSeqAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
    typedef FGroupTaskAlgorithm<GroupOctreeClass, typename GroupOctreeClass::CellGroupClass, GroupCellClass, GroupKernelClass, typename GroupOctreeClass::ParticleGroupClass, GroupContainerClass > GroupAlgorithm;
#endif
    // Get params
    const int NbLevels      = FParameters::getValue(argc,argv,FParameterDefinitions::OctreeHeight.options, 5);
    const int groupSize     = FParameters::getValue(argc,argv,LocalOptionBlocSize.options, 250);

    // Load the particles
#ifdef RANDOM_PARTICLES
    const bool prolate = FParameters::existParameter(argc,argv,LocalOptionProlate.options);
    const bool prolatenonunif = FParameters::existParameter(argc,argv,LocalOptionProlateNonUnif.options);
    const bool spherical = FParameters::existParameter(argc,argv,LocalOptionSpherical.options);
    const bool sphericalnonunif = FParameters::existParameter(argc,argv,LocalOptionSphericalNonUnif.options);
    const bool nonunif = FParameters::existParameter(argc,argv,LocalOptionNonUnif.options);
    FSphericalRandomLoader<FReal> loader(FParameters::getValue(argc,argv,FParameterDefinitions::NbParticles.options, 2000),
                                  nonunif, false, spherical, prolate, sphericalnonunif, prolatenonunif);

#else
    const char* const filename = FParameters::getStr(argc,argv,FParameterDefinitions::InputFile.options, "../Data/test20k.fma");
    FFmaGenericLoader<FReal> loader(filename);
#endif
    FAssertLF(loader.isOpen());
    FTic timer;

    FP2PParticleContainer<FReal> allParticles;
    for(FSize idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        FPoint<FReal> particlePosition;
        FReal physicalValue;
#ifdef RANDOM_PARTICLES
        physicalValue = 0.10;
        loader.fillParticle(&particlePosition);
#else
        loader.fillParticle(&particlePosition, &physicalValue);
#endif
        allParticles.push(particlePosition, physicalValue);
    }
    std::cout << "Particles loaded in " << timer.tacAndElapsed() << "s\n";

    // Put the data into the tree
    timer.tic();
    GroupOctreeClass groupedTree(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), groupSize, &allParticles);
    groupedTree.printInfoBlocks();
    std::cout << "Tree created in " << timer.tacAndElapsed() << "s\n";

    // Run the algorithm
    const MatrixKernelClass MatrixKernel;
    GroupKernelClass groupkernel(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox(), &MatrixKernel);
    GroupAlgorithm groupalgo(&groupedTree,&groupkernel);

    // Extended for Native vs SimGrid makespans comparison
    timer.tic();
    double start_time = starpu_timing_now();
    groupalgo.execute();
    double end_time = starpu_timing_now();
    std::cout << "Kernel executed in in " << timer.tacAndElapsed() << "s\n";
    std::cout << (end_time - start_time)/1000 << "\n";

    // Validate the result
    if(FParameters::existParameter(argc, argv, LocalOptionNoValidate.options) == false){
        FSize offsetParticles = 0;
        FReal*const allPhysicalValues = allParticles.getPhysicalValues();
        FReal*const allPosX = const_cast<FReal*>( allParticles.getPositions()[0]);
        FReal*const allPosY = const_cast<FReal*>( allParticles.getPositions()[1]);
        FReal*const allPosZ = const_cast<FReal*>( allParticles.getPositions()[2]);

        groupedTree.forEachCellLeaf<FP2PGroupParticleContainer<FReal> >([&](GroupCellClass cellTarget, FP2PGroupParticleContainer<FReal> * leafTarget){
            const FReal*const physicalValues = leafTarget->getPhysicalValues();
            const FReal*const posX = leafTarget->getPositions()[0];
            const FReal*const posY = leafTarget->getPositions()[1];
            const FReal*const posZ = leafTarget->getPositions()[2];
            const FSize nbPartsInLeafTarget = leafTarget->getNbParticles();

            for(FSize idxPart = 0 ; idxPart < nbPartsInLeafTarget ; ++idxPart){
                allPhysicalValues[offsetParticles + idxPart] = physicalValues[idxPart];
                allPosX[offsetParticles + idxPart] = posX[idxPart];
                allPosY[offsetParticles + idxPart] = posY[idxPart];
                allPosZ[offsetParticles + idxPart] = posZ[idxPart];
            }

            offsetParticles += nbPartsInLeafTarget;
        });

        FAssertLF(offsetParticles == loader.getNumberOfParticles());

        FReal*const allDirectPotentials = allParticles.getPotentials();
        FReal*const allDirectforcesX = allParticles.getForcesX();
        FReal*const allDirectforcesY = allParticles.getForcesY();
        FReal*const allDirectforcesZ = allParticles.getForcesZ();

        for(int idxTgt = 0 ; idxTgt < offsetParticles ; ++idxTgt){
            for(int idxMutual = idxTgt + 1 ; idxMutual < offsetParticles ; ++idxMutual){
                FP2PR::MutualParticles(
                    allPosX[idxTgt],allPosY[idxTgt],allPosZ[idxTgt], allPhysicalValues[idxTgt],
                    &allDirectforcesX[idxTgt], &allDirectforcesY[idxTgt], &allDirectforcesZ[idxTgt], &allDirectPotentials[idxTgt],
                    allPosX[idxMutual],allPosY[idxMutual],allPosZ[idxMutual], allPhysicalValues[idxMutual],
                    &allDirectforcesX[idxMutual], &allDirectforcesY[idxMutual], &allDirectforcesZ[idxMutual], &allDirectPotentials[idxMutual]
                );
            }
        }

        FMath::FAccurater<FReal> potentialDiff;
        FMath::FAccurater<FReal> fx, fy, fz;
        offsetParticles = 0;
        groupedTree.forEachCellLeaf<FP2PGroupParticleContainer<FReal> >([&](GroupCellClass cellTarget, FP2PGroupParticleContainer<FReal> * leafTarget){
            const FReal*const potentials = leafTarget->getPotentials();
            const FReal*const forcesX = leafTarget->getForcesX();
            const FReal*const forcesY = leafTarget->getForcesY();
            const FReal*const forcesZ = leafTarget->getForcesZ();
            const FSize nbPartsInLeafTarget = leafTarget->getNbParticles();

            for(int idxTgt = 0 ; idxTgt < nbPartsInLeafTarget ; ++idxTgt){
                potentialDiff.add(allDirectPotentials[idxTgt + offsetParticles], potentials[idxTgt]);
                fx.add(allDirectforcesX[idxTgt + offsetParticles], forcesX[idxTgt]);
                fy.add(allDirectforcesY[idxTgt + offsetParticles], forcesY[idxTgt]);
                fz.add(allDirectforcesZ[idxTgt + offsetParticles], forcesZ[idxTgt]);
            }

            offsetParticles += nbPartsInLeafTarget;
        });

        std::cout << "Error : Potential " << potentialDiff << "\n";
        std::cout << "Error : fx " << fx << "\n";
        std::cout << "Error : fy " << fy << "\n";
        std::cout << "Error : fz " << fz << "\n";
    }

    return 0;
}
}
