// ===================================================================================
// Ce LOGICIEL "ScalFmm" est couvert par le copyright Inria 20xx-2012.
// Inria détient tous les droits de propriété sur le LOGICIEL, et souhaite que
// la communauté scientifique l'utilise afin de le tester et de l'évaluer.
// Inria donne gracieusement le droit d'utiliser ce LOGICIEL. Toute utilisation
// dans un but lucratif ou à des fins commerciales est interdite sauf autorisation
// expresse et préalable d'Inria.
// Toute utilisation hors des limites précisées ci-dessus et réalisée sans l'accord
// expresse préalable d'Inria constituerait donc le délit de contrefaçon.
// Le LOGICIEL étant un produit en cours de développement, Inria ne saurait assurer
// aucune responsabilité et notamment en aucune manière et en aucun cas, être tenu
// de répondre d'éventuels dommages directs ou indirects subits par l'utilisateur.
// Tout utilisateur du LOGICIEL s'engage à communiquer à Inria ses remarques
// relatives à l'usage du LOGICIEL
// ===================================================================================
#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Fmb/FFmbComponents.hpp"

#include "../Src/Kernels/FSphericalCell.hpp"

#include "../Src/Kernels/FSphericalKernel.hpp"

#include "../Src/Core/FFmmAlgorithmPeriodic.hpp"

#include "FUTester.hpp"

#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestPeriodicKernels.hpp"
#include "../Src/Components/FTestParticle.hpp"

/*
  In this test we compare the fmm results and the direct results.
  */
class IndexedParticle : public FmbParticle {
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
class TestFmbDirectPeriodic : public FUTester<TestFmbDirectPeriodic> {
    /** The test */
    void TestDirect(){
        typedef IndexedParticle         ParticleClass;
        typedef FSphericalCell            CellClass;
        typedef FVector<ParticleClass>  ContainerClass;

        typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >   KernelClass;

        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithmPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        const int NbLevels      = 2;
        const int SizeSubLevels = 1;
        const int DevP = 12;
        const int PeriodicDeep = 2;

        const FReal BoxWidth = 1;
        const long NbSmallBoxesPerSide = (1 << (NbLevels-1));
        const FReal SmallBoxWidth = BoxWidth / FReal(NbSmallBoxesPerSide);
        const FReal SmallBoxWidthDiv2 = SmallBoxWidth / 2;

        const int NbPart = NbSmallBoxesPerSide * NbSmallBoxesPerSide * NbSmallBoxesPerSide;

        FSphericalCell::Init(DevP);

        // Create octree
        OctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, F3DPosition(0.5,0.5,0.5));
        {
            int idxPart = 0;
            for(int idxX = 0 ; idxX < NbSmallBoxesPerSide ; ++idxX){
                for(int idxY = 0 ; idxY < NbSmallBoxesPerSide ; ++idxY){
                    for(int idxZ = 0 ; idxZ < NbSmallBoxesPerSide ; ++idxZ){
                        ParticleClass particleToFill;
                        particleToFill.setPosition(FReal(idxX)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                   FReal(idxY)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001),
                                                   FReal(idxZ)*SmallBoxWidth + SmallBoxWidthDiv2 + FReal(0.001));
                        particleToFill.setPhysicalValue(FReal(0.01) /*+ FReal(idxX) * FReal(0.01) + FReal(idxY) * FReal(0.07) + FReal(idxZ) * FReal(0.013)*/);
                        particleToFill.setIndex(idxPart++);
                        tree.insert(particleToFill);
                        idxX = idxY = idxZ =  NbSmallBoxesPerSide;//todo delete
                    }
                }
            }
        }

        // Run FMM
        Print("Fmm...");
        KernelClass kernels( DevP, NbLevels, BoxWidth, PeriodicDeep);
        FmmClass algo(&tree,&kernels,PeriodicDeep);
        algo.execute();

        Print("Prepare direct...");
        // Run direct computation
        const int NbBoxPerPeriodicSide = 6 * PeriodicDeep;
        const int directNbPart = NbPart * NbBoxPerPeriodicSide * NbBoxPerPeriodicSide * NbBoxPerPeriodicSide;
        ParticleClass* const particles = new ParticleClass[directNbPart];
        {
            for(int idxBoxX = 0 ; idxBoxX < NbBoxPerPeriodicSide ; ++idxBoxX){
                for(int idxBoxY = 0 ; idxBoxY < NbBoxPerPeriodicSide ; ++idxBoxY){
                    for(int idxBoxZ = 0 ; idxBoxZ < NbBoxPerPeriodicSide ; ++idxBoxZ){
                        const FReal xoffset = FReal(idxBoxX) * BoxWidth;
                        const FReal yoffset = FReal(idxBoxY) * BoxWidth;
                        const FReal zoffset = FReal(idxBoxZ) * BoxWidth;

                        ParticleClass*const partBox = &particles[ NbPart * ((idxBoxX * NbBoxPerPeriodicSide * NbBoxPerPeriodicSide) + (idxBoxY * NbBoxPerPeriodicSide) + idxBoxZ)];
                        int idxPart = 0;

                        for(int idxX = 0 ; idxX < NbSmallBoxesPerSide ; ++idxX){
                            for(int idxY = 0 ; idxY < NbSmallBoxesPerSide ; ++idxY){
                                for(int idxZ = 0 ; idxZ < NbSmallBoxesPerSide ; ++idxZ){
                                    partBox[idxPart].setPosition(FReal(idxX)*SmallBoxWidth + SmallBoxWidthDiv2 + xoffset + FReal(0.001),
                                                               FReal(idxY)*SmallBoxWidth + SmallBoxWidthDiv2 + yoffset + FReal(0.001),
                                                               FReal(idxZ)*SmallBoxWidth + SmallBoxWidthDiv2 + zoffset + FReal(0.001)   );
                                    partBox[idxPart].setPhysicalValue(FReal(0.01) /*+ FReal(idxX) * FReal(0.01) + FReal(idxY) * FReal(0.07) + FReal(idxZ) * FReal(0.013)*/);
                                    idxX = idxY = idxZ =  NbSmallBoxesPerSide;//todo delete
                                    ++idxPart;
                                }
                            }
                        }
                    }
                }
            }
        }

        Print("Direct...");
        for(int idxTarget = 0 ; idxTarget < directNbPart ; idxTarget += 8){
            for(int idxOther = idxTarget + 8 ; idxOther < directNbPart ; idxOther += 8){
                kernels.directInteractionMutual(&particles[idxTarget], &particles[idxOther]);
            }
        }//todo exchange
        /*for(int idxTarget = 0 ; idxTarget < directNbPart ; ++idxTarget){
            for(int idxOther = idxTarget + 1 ; idxOther < directNbPart ; ++idxOther){
                kernels.directInteractionMutual(&particles[idxTarget], &particles[idxOther]);
            }
        }*/

        // Compare
        Print("Compute Diff...");
        FReal potentialDiff = 0;
        FReal fx = 0, fy = 0, fz = 0;
        { // Check that each particle has been summed with all other
            const int middle = (NbBoxPerPeriodicSide/2) - PeriodicDeep;
            const int boxStartIdx = NbPart * ((middle * NbBoxPerPeriodicSide * NbBoxPerPeriodicSide) + (middle * NbBoxPerPeriodicSide) + middle);
            ParticleClass*const partBox = &particles[boxStartIdx];

            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                typename ContainerClass::BasicIterator leafIter(*octreeIterator.getCurrentListTargets());

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


                    const FReal currentPotentialDiff = FMath::RelativeDiff(other.getPotential(),leafIter.data().getPotential());
                    if( potentialDiff < currentPotentialDiff ){
                        potentialDiff = currentPotentialDiff;
                    }

                    const FReal currentFx = FMath::RelativeDiff(other.getForces().getX() , leafIter.data().getForces().getX());
                    if( fx < currentFx ){
                        fx = currentFx;
                    }

                    const FReal currentFy = FMath::RelativeDiff(other.getForces().getY() , leafIter.data().getForces().getY());
                    if( fy < currentFy ){
                        fy = currentFy;
                    }

                    const FReal currentFz = FMath::RelativeDiff(other.getForces().getZ() , leafIter.data().getForces().getZ());
                    if( fz < currentFz ){
                        fz = currentFz;
                    }

                    leafIter.gotoNext();
                }
            } while(octreeIterator.moveRight());
        }


        delete[] particles;

        Print("Potential diff is = ");
        Print(potentialDiff);
        Print("Fx diff is = ");
        Print(fx);
        Print("Fy diff is = ");
        Print(fy);
        Print("Fz diff is = ");
        Print(fz);

        const FReal MaximumDiff = FReal(0.5);
        assert(potentialDiff < MaximumDiff);
        assert(fx < MaximumDiff);
        assert(fy < MaximumDiff);
        assert(fz < MaximumDiff);
    }


    /** The test */
    void TestDirectSum(){
        typedef FTestParticle         TestParticleClass;
        typedef FTestCell             TestCellClass;
        typedef FVector<TestParticleClass>  TestContainerClass;

        typedef FTestPeriodicKernels<TestParticleClass, TestCellClass, TestContainerClass >          TestKernelClass;

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

        const int NbPart = NbSmallBoxesPerSide * NbSmallBoxesPerSide * NbSmallBoxesPerSide;

        // Create octree
        TestOctreeClass tree(NbLevels, SizeSubLevels, BoxWidth, F3DPosition(0.5,0.5,0.5));
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

            typename TestOctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                typename TestContainerClass::BasicIterator leafIter(*octreeIterator.getCurrentListTargets());

                while( leafIter.hasNotFinished() ){
                    assert(leafIter.data().getDataDown() == directNbPart - 1);// todo delete


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
        AddTest(&TestFmbDirectPeriodic::TestDirectSum,"Test Simu and with direct compare to Test fmm periodic");
        AddTest(&TestFmbDirectPeriodic::TestDirect,"Test Simu and with direct compare to real fmm periodic");
    }
};


// You must do this
TestClass(TestFmbDirectPeriodic)



