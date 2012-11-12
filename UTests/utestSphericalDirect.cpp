// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
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

#include "../Src/Utils/FGlobal.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../Src/Kernels/Spherical/FSphericalParticle.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"
#include "../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../Src/Kernels/Spherical/FSphericalRotationKernel.hpp"
#include "../Src/Kernels/Spherical/FSphericalBlasKernel.hpp"
#include "../Src/Kernels/Spherical/FSphericalBlockBlasKernel.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FTreeIO.hpp"

#include "../Src/Core/FFmmAlgorithmThread.hpp"
#include "../Src/Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"

/*
  In this test we compare the spherical fmm results and the direct results.
  */

/** We need to know the position of the particle in the array */
class IndexedParticle : public FSphericalParticle {
    int index;
public:
    IndexedParticle(): index(-1){}

    int getIndex() const{
        return index;
    }
    void setIndex( const int inIndex ){
        index = inIndex;
    }
};

/** the test class
  *
  */
class TestSphericalDirect : public FUTester<TestSphericalDirect> {
    /** The test method to factorize all the test based on different kernels */
    template <class ParticleClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass,
              class OctreeClass, class FmmClass>
    void RunTest(const bool isBlasKernel){
        // Warning in make test the exec dir it Build/UTests
        // Load particles
        const char* const filename = (sizeof(FReal) == sizeof(float))?
                                        "../../Data/utestDirect.bin.fma.single":
                                        "../../Data/utestDirect.bin.fma.double";
        FFmaBinLoader<ParticleClass> loader(filename);
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            uassert(false);
            return;
        }
        Print("Number of particles:");
        Print(loader.getNumberOfParticles());

        const int NbLevels      = 4;
        const int SizeSubLevels = 2;
        const int DevP = 9;
        FSphericalCell::Init(DevP, isBlasKernel);

        // Create octree
        OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        ParticleClass* const particles = new ParticleClass[loader.getNumberOfParticles()];
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(particles[idxPart]);
            particles[idxPart].setIndex( idxPart );
            tree.insert(particles[idxPart]);
        }


        // Run FMM
        Print("Fmm...");
        //KernelClass kernels(NbLevels,loader.getBoxWidth());
        KernelClass kernels(DevP,NbLevels,loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClass algo(&tree,&kernels);
        algo.execute();

        // Run direct computation
        Print("Direct...");
        for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
            for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                kernels.directInteractionMutual(&particles[idxTarget], &particles[idxOther]);
            }
        }

        // Compare
        Print("Compute Diff...");
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;
        { // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                typename ContainerClass::BasicIterator leafIter(*octreeIterator.getCurrentListTargets());

                while( leafIter.hasNotFinished() ){
                    const ParticleClass& other = particles[leafIter.data().getIndex()];

                    potentialDiff.add(other.getPotential(),leafIter.data().getPotential());

                    fx.add(other.getForces().getX(),leafIter.data().getForces().getX());

                    fy.add(other.getForces().getY(),leafIter.data().getForces().getY());

                    fz.add(other.getForces().getZ(),leafIter.data().getForces().getZ());

                    leafIter.gotoNext();
                }
            } while(octreeIterator.moveRight());
        }

        delete[] particles;

        // Print for information
        Print("Potential diff is = ");
        Print(potentialDiff.getL2Norm());
        Print(potentialDiff.getInfNorm());
        Print("Fx diff is = ");
        Print(fx.getL2Norm());
        Print(fx.getInfNorm());
        Print("Fy diff is = ");
        Print(fy.getL2Norm());
        Print(fy.getInfNorm());
        Print("Fz diff is = ");
        Print(fz.getL2Norm());
        Print(fz.getInfNorm());

        // Assert
        const FReal MaximumDiff = FReal(0.0001);
        uassert(potentialDiff.getL2Norm() < MaximumDiff);
        uassert(potentialDiff.getInfNorm() < MaximumDiff);
        uassert(fx.getL2Norm()  < MaximumDiff);
        uassert(fx.getInfNorm() < MaximumDiff);
        uassert(fy.getL2Norm()  < MaximumDiff);
        uassert(fy.getInfNorm() < MaximumDiff);
        uassert(fz.getL2Norm()  < MaximumDiff);
        uassert(fz.getInfNorm() < MaximumDiff);
    }

    /** If memstas is running print the memory used */
    void PostTest() {
        if( FMemStats::controler.isUsed() ){
            std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated() << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated() << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated() << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }

    ///////////////////////////////////////////////////////////
    // The tests!
    ///////////////////////////////////////////////////////////

    /** Classic */
    void TestSpherical(){
        typedef IndexedParticle         ParticleClass;
        typedef FSphericalCell            CellClass;
        typedef FVector<ParticleClass>  ContainerClass;

        typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest<ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass,
                OctreeClass, FmmClass>(false);
    }

    /** Rotation */
    void TestRotation(){
        typedef IndexedParticle         ParticleClass;
        typedef FSphericalCell            CellClass;
        typedef FVector<ParticleClass>  ContainerClass;

        typedef FSphericalRotationKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest<ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass,
                OctreeClass, FmmClass>(false);
    }

#ifdef SCALFMM_USE_BLAS
    /** Blas */
    void TestSphericalBlas(){
        typedef IndexedParticle         ParticleClass;
        typedef FSphericalCell            CellClass;
        typedef FVector<ParticleClass>  ContainerClass;

        typedef FSphericalBlasKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest<ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass,
                OctreeClass, FmmClass>(true);
    }

    /** Block blas */
    void TestSphericalBlockBlas(){
        typedef IndexedParticle         ParticleClass;
        typedef FSphericalCell            CellClass;
        typedef FVector<ParticleClass>  ContainerClass;

        typedef FSphericalBlockBlasKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        RunTest<ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass,
                OctreeClass, FmmClass>(true);
    }
#endif

    ///////////////////////////////////////////////////////////
    // Set the tests!
    ///////////////////////////////////////////////////////////

    /** set test */
    void SetTests(){
        AddTest(&TestSphericalDirect::TestSpherical,"Test Spherical Kernel");
        AddTest(&TestSphericalDirect::TestRotation,"Test Rotation Spherical Kernel");
#ifdef SCALFMM_USE_BLAS
        AddTest(&TestSphericalDirect::TestSphericalBlas,"Test Spherical Blas Kernel");
        AddTest(&TestSphericalDirect::TestSphericalBlockBlas,"Test Spherical Block Blas Kernel");
#endif
    }
};


// You must do this
TestClass(TestSphericalDirect)



