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

#include "../Src/Kernels/Rotation/FRotationCell.hpp"
#include "../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../Src/Kernels/Rotation/FRotationParticle.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"
#include "../Src/Core/FFmmAlgorithmPeriodic.hpp"
#include "../Src/Files/FRandomLoader.hpp"

#include "FUTester.hpp"

#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestKernels.hpp"
#include "../Src/Components/FTestParticle.hpp"

/*
  In this test we compare the fmm results and the direct results.
  */
class IndexedParticle : public FRotationParticle {
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
class TestRotationDirectPeriodic : public FUTester<TestRotationDirectPeriodic> {
    /** Here we test only the P2P */
    void TestPeriodicFmm(){
        static const int P = 9;
        typedef IndexedParticle         ParticleClass;
        typedef FRotationCell<P>            CellClass;
        typedef FVector<ParticleClass>  ContainerClass;

        typedef FRotationKernel<ParticleClass, CellClass, ContainerClass, P >   KernelClass;

        typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
        typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithmPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        // Parameters
        const int NbLevels      = 4;
        const int SizeSubLevels = 2;
        const int PeriodicDeep  = -1;
        const int NbParticles   = 100;
        ParticleClass* const particles = new ParticleClass[NbParticles];

        FRandomLoader<ParticleClass> loader(NbParticles);
        OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        for( int idxPart = 0 ; idxPart < NbParticles ; ++idxPart ){
            loader.fillParticle(particles[idxPart]);
            particles[idxPart].setIndex(idxPart);
            particles[idxPart].setPhysicalValue(FReal(0.80));
            tree.insert(particles[idxPart]);
        }

        // Run FMM
        Print("Fmm...");
        FmmClass algo(&tree,PeriodicDeep, DirX | DirMinusY  | DirPlusZ);
        KernelClass kernels( algo.extendedTreeHeight(), algo.extendedBoxWidth(), algo.extendedBoxCenter());
        algo.setKernel(&kernels);
        algo.execute();

        // Run Direct
        Print("Run direct...");
        FTreeCoordinate min, max;
        algo.repetitionsIntervals(&min, &max);

        for(int idxTarget = 0 ; idxTarget < NbParticles ; ++idxTarget){
            for(int idxSource = idxTarget + 1 ; idxSource < NbParticles ; ++idxSource){
                kernels.particlesMutualInteraction(&particles[idxTarget], &particles[idxSource]);
            }
            for(int idxX = min.getX() ; idxX <= max.getX() ; ++idxX){
                for(int idxY = min.getY() ; idxY <= max.getY() ; ++idxY){
                    for(int idxZ = min.getZ() ; idxZ <= max.getZ() ; ++idxZ){
                        if(idxX == 0 && idxY == 0 && idxZ == 0) continue;

                        const FPoint offset(loader.getBoxWidth() * FReal(idxX),
                                            loader.getBoxWidth() * FReal(idxY),
                                            loader.getBoxWidth() * FReal(idxZ));

                        for(int idxSource = 0 ; idxSource < NbParticles ; ++idxSource){
                            ParticleClass source = particles[idxSource];
                            source.incPosition(offset.getX(),offset.getY(),offset.getZ());
                            kernels.particlesInteraction(&particles[idxTarget], source);
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
                    const ParticleClass& other = particles[leafIter.data().getIndex()];

                    /*printf("Tree x %e y %e z %e physical %e potential %e fx %e fy %e fz %e\n",
                           leafIter.data().getPosition().getX(),leafIter.data().getPosition().getY(),leafIter.data().getPosition().getZ(),
                           leafIter.data().getPhysicalValue(),leafIter.data().getPotential(),
                           leafIter.data().getForces().getX(),leafIter.data().getForces().getY(),leafIter.data().getForces().getZ());// todo delete
                    printf("Direct x %e y %e z %e physical %e potential %e fx %e fy %e fz %e\n",
                           other.getPosition().getX(),other.getPosition().getY(),other.getPosition().getZ(),
                           other.getPhysicalValue(),other.getPotential(),
                           other.getForces().getX(),other.getForces().getY(),other.getForces().getZ());//*/

                    potentialDiff.add(other.getPotential(),leafIter.data().getPotential());

                    fx.add(other.getForces().getX(),leafIter.data().getForces().getX());

                    fy.add(other.getForces().getY(),leafIter.data().getForces().getY());

                    fz.add(other.getForces().getZ(),leafIter.data().getForces().getZ());

                    leafIter.gotoNext();
                }
            } while(octreeIterator.moveRight());
        }

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

        delete[] particles;
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
        AddTest(&TestRotationDirectPeriodic::TestPeriodicFmm,"Test Simu and with direct compare to Test fmm periodic");
    }
};


// You must do this
TestClass(TestRotationDirectPeriodic)



