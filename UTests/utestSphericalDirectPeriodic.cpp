// ===================================================================================
// Logiciel initial: ScalFmm Version 0.5
// Co-auteurs : Olivier Coulaud, Bérenger Bramas.
// Propriétaires : INRIA.
// Copyright © 2011-2012, diffusé sous les termes et conditions d’une licence propriétaire.
// Initial software: ScalFmm Version 0.5
// Co-authors: Olivier Coulaud, Bérenger Bramas.
// Owners: INRIA.
// Copyright © 2011-2012, spread under the terms and conditions of a proprietary license.
// ===================================================================================
#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../Src/Kernels/Spherical/FSphericalParticle.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"
#include "../Src/Core/FFmmAlgorithmPeriodic.hpp"

#include "FUTester.hpp"

#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestKernels.hpp"
#include "../Src/Components/FTestParticle.hpp"

/*
  In this test we compare the fmm results and the direct results.
  */
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

/** The class to run the test */
class TestSphericalDirectPeriodic : public FUTester<TestSphericalDirectPeriodic> {
    /** Here we test only the P2P */
    void TestPeriodicP2P(){
        typedef IndexedParticle         ParticleClass;
        typedef FSphericalCell            CellClass;
        typedef FVector<ParticleClass>  ContainerClass;

        typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >   KernelClass;

        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithmPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        // Parameters
        const int NbLevels      = 2;
        const int SizeSubLevels = 1;
        const int DevP = 12;
        const int PeriodicDeep = 0;

        const FReal BoxWidth = 1;
        const long NbSmallBoxesPerSide = (1 << (NbLevels-1));
        const FReal SmallBoxWidth = BoxWidth / FReal(NbSmallBoxesPerSide);
        const FReal SmallBoxWidthDiv2 = SmallBoxWidth / 2;
        const F3DPosition CenterOfBox = F3DPosition(0.5,0.5,0.5);

        FSphericalCell::Init(DevP);

        // Create octree
        OctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, CenterOfBox);
        { // Insert one particle in each leaf
            int idxPart = 0;
            for(int idxX = 0 ; idxX < 2 ; ++idxX){
                for(int idxY = 0 ; idxY < 2 ; ++idxY){
                    for(int idxZ = 0 ; idxZ < 2 ; ++idxZ){
                        ParticleClass particleToFill;
                        particleToFill.setPosition(FReal(idxX)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                   FReal(idxY)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                   FReal(idxZ)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001));
                        particleToFill.setPhysicalValue(FReal(0.01));
                        particleToFill.setIndex(idxPart++);
                        tree.insert(particleToFill);
                    }
                }
            }
        }

        // Run P2P
        Print("Fmm...");
        KernelClass kernels( DevP, NbLevels, BoxWidth, CenterOfBox, PeriodicDeep + 1);
        FmmClass algo(&tree,&kernels,PeriodicDeep);
        algo.directPass();

        Print("Prepare direct...");
        // Run direct computation
        const int directNbPart =  4 * 4 * 4;
        ParticleClass* const particles = new ParticleClass[directNbPart];
        {
            for(int idxBoxX = 0 ; idxBoxX < 4 ; ++idxBoxX){
                for(int idxBoxY = 0 ; idxBoxY < 4 ; ++idxBoxY){
                    for(int idxBoxZ = 0 ; idxBoxZ < 4 ; ++idxBoxZ){
                        const int indexPart = ((4 * idxBoxX) + idxBoxY) * 4 + idxBoxZ;

                        particles[indexPart].setPosition(FReal(idxBoxX - 1)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                               FReal(idxBoxY - 1)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                               FReal(idxBoxZ - 1)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001));
                        particles[indexPart].setPhysicalValue(FReal(0.01));

                    }
                }
            }
        }

        // Run direct
        Print("Direct...");
        for(int idxBoxX = 1 ; idxBoxX < 3 ; ++idxBoxX){
            for(int idxBoxY = 1 ; idxBoxY < 3 ; ++idxBoxY){
                for(int idxBoxZ = 1 ; idxBoxZ < 3 ; ++idxBoxZ){
                    const int indexPart = ((4 * idxBoxX) + idxBoxY) * 4 + idxBoxZ;

                    for(int idxBoxXSrc = idxBoxX - 1 ; idxBoxXSrc <= idxBoxX + 1 ; ++idxBoxXSrc){
                        for(int idxBoxYSrc = idxBoxY - 1 ; idxBoxYSrc <= idxBoxY + 1 ; ++idxBoxYSrc){
                            for(int idxBoxZSrc = idxBoxZ - 1 ; idxBoxZSrc <= idxBoxZ + 1 ; ++idxBoxZSrc){
                                const int indexPartSrc = ((4 * idxBoxXSrc) + idxBoxYSrc) * 4 + idxBoxZSrc;

                                if( indexPart != indexPartSrc ){
                                    kernels.directInteraction(&particles[indexPart], particles[indexPartSrc]);
                                }
                            }
                        }
                    }

                }
            }
        }

        // Compare
        Print("Compute Diff...");
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;
        { // Check that each particle has been summed with all other
            OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                ContainerClass::BasicIterator leafIter(*octreeIterator.getCurrentListTargets());

                while( leafIter.hasNotFinished() ){
                    const ParticleClass& other = particles[((4 * (octreeIterator.getCurrentGlobalCoordinate().getX()+1)) +
                                                            (octreeIterator.getCurrentGlobalCoordinate().getY()+1)) * 4 +
                                                            (octreeIterator.getCurrentGlobalCoordinate().getZ()+1)];

                    printf("Tree x %e y %e z %e physical %e potential %e fx %e fy %e fz %e\n",
                           leafIter.data().getPosition().getX(),leafIter.data().getPosition().getY(),leafIter.data().getPosition().getZ(),
                           leafIter.data().getPhysicalValue(),leafIter.data().getPotential(),
                           leafIter.data().getForces().getX(),leafIter.data().getForces().getY(),leafIter.data().getForces().getZ());// todo delete
                    printf("Direct x %e y %e z %e physical %e potential %e fx %e fy %e fz %e\n",
                           other.getPosition().getX(),other.getPosition().getY(),other.getPosition().getZ(),
                           other.getPhysicalValue(),other.getPotential(),
                           other.getForces().getX(),other.getForces().getY(),other.getForces().getZ());// todo delete


                    potentialDiff.add(other.getPotential(),leafIter.data().getPotential());

                    fx.add(other.getForces().getX(),leafIter.data().getForces().getX());

                    fy.add(other.getForces().getY(),leafIter.data().getForces().getY());

                    fz.add(other.getForces().getZ(),leafIter.data().getForces().getZ());

                    leafIter.gotoNext();
                }
            } while(octreeIterator.moveRight());
        }


        delete[] particles;

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

    /** Test real Periodic FMM */
    void TestDirectHigh(){
        typedef IndexedParticle         ParticleClass;
        typedef FSphericalCell            CellClass;
        typedef FVector<ParticleClass>  ContainerClass;

        typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >   KernelClass;

        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithmPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        // Parameters
        const int NbLevels      = 2;
        const int SizeSubLevels = 1;
        const int DevP = 8;
        const int PeriodicDeep = 2;

        const FReal BoxWidth = 1;
        const long NbSmallBoxesPerSide = (1 << (NbLevels-1));
        const FReal SmallBoxWidth = BoxWidth / FReal(NbSmallBoxesPerSide);
        const FReal SmallBoxWidthDiv2 = SmallBoxWidth / 2;
        const F3DPosition CenterOfBox = F3DPosition(0.5,0.5,0.5);

        FSphericalCell::Init(DevP);

        // Create octree
        OctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, CenterOfBox);
        {
            int idxPart = 0;
            for(int idxX = 0 ; idxX < NbSmallBoxesPerSide ; ++idxX){
                for(int idxY = 0 ; idxY < NbSmallBoxesPerSide ; ++idxY){
                    for(int idxZ = 0 ; idxZ < NbSmallBoxesPerSide ; ++idxZ){
                        ParticleClass particleToFill;
                        particleToFill.setPosition(FReal(idxX)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                   FReal(idxY)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                   FReal(idxZ)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001));
                        particleToFill.setPhysicalValue(FReal(0.01));
                        particleToFill.setIndex(idxPart++);
                        tree.insert(particleToFill);
                    }
                }
            }
        }

        // Run FMM
        Print("Fmm...");
        KernelClass kernels( DevP, NbLevels, BoxWidth, CenterOfBox, PeriodicDeep);
        FmmClass algo( &tree, &kernels, 0);
        algo.execute();

        Print("Prepare direct...");
        // Run direct computation
        const int NbBoxPerPeriodicSide = 12;
        const int directNbPart = 1 * NbBoxPerPeriodicSide * NbBoxPerPeriodicSide * NbBoxPerPeriodicSide;
        ParticleClass* const particles = new ParticleClass[directNbPart];
        {
            for(int idxBoxX = 0 ; idxBoxX < NbBoxPerPeriodicSide ; ++idxBoxX){
                for(int idxBoxY = 0 ; idxBoxY < NbBoxPerPeriodicSide ; ++idxBoxY){
                    for(int idxBoxZ = 0 ; idxBoxZ < NbBoxPerPeriodicSide ; ++idxBoxZ){
                        const int indexPart = ((NbBoxPerPeriodicSide * idxBoxX) + idxBoxY) * NbBoxPerPeriodicSide + idxBoxZ;

                        particles[indexPart].setPosition(FReal(idxBoxX - 1)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                               FReal(idxBoxY - 1)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                               FReal(idxBoxZ - 1)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001));
                        particles[indexPart].setPhysicalValue(FReal(0.01));

                    }
                }
            }
        }

        Print("Direct...");
        for(int idxBoxX = 4 ; idxBoxX <= 5 ; ++idxBoxX){
            for(int idxBoxY = 4 ; idxBoxY <= 5 ; ++idxBoxY){
                for(int idxBoxZ = 4 ; idxBoxZ <= 5 ; ++idxBoxZ){
                    const int indexPart = ((NbBoxPerPeriodicSide * idxBoxX) + idxBoxY) * NbBoxPerPeriodicSide + idxBoxZ;

                    for(int idxBoxXSrc = 0 ; idxBoxXSrc <= NbBoxPerPeriodicSide ; ++idxBoxXSrc){
                        for(int idxBoxYSrc = 0 ; idxBoxYSrc <= NbBoxPerPeriodicSide ; ++idxBoxYSrc){
                            for(int idxBoxZSrc = 0 ; idxBoxZSrc <= NbBoxPerPeriodicSide ; ++idxBoxZSrc){
                                const int indexPartSrc = ((NbBoxPerPeriodicSide * idxBoxXSrc) + idxBoxYSrc) * NbBoxPerPeriodicSide + idxBoxZSrc;

                                if( indexPart != indexPartSrc ){
                                    kernels.directInteraction(&particles[indexPart], particles[indexPartSrc]);
                                }
                            }
                        }
                    }

                }
            }
        }

        // Compare
        Print("Compute Diff...");
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;
        { // Check that each particle has been summed with all other
            ParticleClass*const partBox = 0;//&particles[boxStartIdx];

            OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                ContainerClass::BasicIterator leafIter(*octreeIterator.getCurrentListTargets());

                while( leafIter.hasNotFinished() ){
                    const ParticleClass& other = partBox[leafIter.data().getIndex()];

                    printf("Tree x %e y %e z %e physical %e potential %e fx %e fy %e fz %e\n",
                           leafIter.data().getPosition().getX(),leafIter.data().getPosition().getY(),leafIter.data().getPosition().getZ(),
                           leafIter.data().getPhysicalValue(),leafIter.data().getPotential(),
                           leafIter.data().getForces().getX(),leafIter.data().getForces().getY(),leafIter.data().getForces().getZ());// todo delete
                    printf("Direct x %e y %e z %e physical %e potential %e fx %e fy %e fz %e\n",
                           other.getPosition().getX(),other.getPosition().getY(),other.getPosition().getZ(),
                           other.getPhysicalValue(),other.getPotential(),
                           other.getForces().getX(),other.getForces().getY(),other.getForces().getZ());// todo delete


                    potentialDiff.add(other.getPotential(),leafIter.data().getPotential());

                    fx.add(other.getForces().getX(),leafIter.data().getForces().getX());

                    fy.add(other.getForces().getY(),leafIter.data().getForces().getY());

                    fz.add(other.getForces().getZ(),leafIter.data().getForces().getZ());

                    leafIter.gotoNext();
                }
            } while(octreeIterator.moveRight());
        }


        delete[] particles;

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


    /** The test */
    void TestPeriodicFmm(){
        typedef FTestParticle         TestParticleClass;
        typedef FTestCell             TestCellClass;
        typedef FVector<TestParticleClass>  TestContainerClass;

        typedef FTestKernels<TestParticleClass, TestCellClass, TestContainerClass >          TestKernelClass;

        typedef FSimpleLeaf<TestParticleClass, TestContainerClass >                     TestLeafClass;
        typedef FOctree<TestParticleClass, TestCellClass, TestContainerClass , TestLeafClass >  TestOctreeClass;

        typedef FFmmAlgorithmPeriodic<TestOctreeClass, TestParticleClass, TestCellClass, TestContainerClass, TestKernelClass, TestLeafClass > TestFmmClass;

        const int NbLevels      = 2;
        const int SizeSubLevels = 1;
        const int PeriodicDeep = 2;

        const FReal BoxWidth = 1;
        const long NbSmallBoxesPerSide = (1 << (NbLevels-1));
        const FReal SmallBoxWidth = BoxWidth / FReal(NbSmallBoxesPerSide);
        const FReal SmallBoxWidthDiv2 = SmallBoxWidth / 2;
        const F3DPosition CenterOfBox = F3DPosition(0.5,0.5,0.5);

        const int NbPart = NbSmallBoxesPerSide * NbSmallBoxesPerSide * NbSmallBoxesPerSide;

        // Create octree
        TestOctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, CenterOfBox);
        {
            for(int idxX = 0 ; idxX < NbSmallBoxesPerSide ; ++idxX){
                for(int idxY = 0 ; idxY < NbSmallBoxesPerSide ; ++idxY){
                    for(int idxZ = 0 ; idxZ < NbSmallBoxesPerSide ; ++idxZ){
                        TestParticleClass particleToFill;
                        particleToFill.setPosition(FReal(idxX)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                   FReal(idxY)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                   FReal(idxZ)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001));
                        tree.insert(particleToFill);
                    }
                }
            }
        }

        // Run FMM
        Print("Fmm...");
        TestKernelClass kernels;
        TestFmmClass algo(&tree,&kernels,PeriodicDeep);
        algo.execute();

        // Compare
        Print("Check down data...");
        { // Check that each particle has been summed with all other
            const int NbBoxPerPeriodicSide = 6 * PeriodicDeep;
            const int directNbPart = NbPart * NbBoxPerPeriodicSide * NbBoxPerPeriodicSide * NbBoxPerPeriodicSide;

            TestOctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                TestContainerClass::BasicIterator leafIter(*octreeIterator.getCurrentListTargets());

                while( leafIter.hasNotFinished() ){
                    uassert(leafIter.data().getDataDown() == directNbPart - 1);// todo delete


                    leafIter.gotoNext();
                }
            } while(octreeIterator.moveRight());
        }
    }

    /** After check the memory if needed */
    void After() {
        if( FMemStats::controler.isUsed() ){
            std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated() << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated() << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated() << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }

    // set test
    void SetTests(){        
        //AddTest(&TestSphericalDirectPeriodic::TestPeriodicFmm,"Test Simu and with direct compare to Test fmm periodic");
        //AddTest(&TestSphericalDirectPeriodic::TestPeriodicP2P,"Test direct compare to real fmm periodic (P2P only)");
        //AddTest(&TestSphericalDirectPeriodic::TestDirectHigh,"Test direct compare to real fmm periodic");
    }
};


// You must do this
TestClass(TestSphericalDirectPeriodic)



