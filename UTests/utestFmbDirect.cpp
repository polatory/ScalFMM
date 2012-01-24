// [--License--]
#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Fmb/FFmbKernels.hpp"
#include "../Src/Fmb/FFmbComponents.hpp"

#include "../Src/Kernels/FComputeCell.hpp"

#include "../Src/Kernels/FElecForcesKernels.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FTreeIO.hpp"

#include "../Src/Core/FFmmAlgorithmThread.hpp"
#include "../Src/Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"

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


class TestFmbDirect : public FUTester<TestFmbDirect> {
    typedef IndexedParticle         ParticleClass;
    //typedef FmbCell            CellClass;
    typedef FComputeCell            CellClass;
    typedef FVector<ParticleClass>  ContainerClass;

    //typedef FFmbKernels<ParticleClass, CellClass, ContainerClass >          KernelClass;
    typedef FElecForcesKernels<ParticleClass, CellClass, ContainerClass >          KernelClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

    typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;


    void TestDirect(){
        FFmaBinLoader<ParticleClass> loader("../Data/utestFmbDirect.bin.fma");
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            return;
        }
        Print("Number of particles:");
        Print(loader.getNumberOfParticles());

        const int NbLevels      = 4;
        const int SizeSubLevels = 2;
        const int DevP = 12;
        FComputeCell::Init(DevP);
        OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());


        ParticleClass* const particles = new ParticleClass[loader.getNumberOfParticles()];
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(particles[idxPart]);
            particles[idxPart].setIndex( idxPart );
            tree.insert(particles[idxPart]);
        }


        Print("Fmm...");
        //KernelClass kernels(NbLevels,loader.getBoxWidth());
        KernelClass kernels(DevP,NbLevels,loader.getBoxWidth());
        FmmClass algo(&tree,&kernels);
        algo.execute();


        Print("Direct...");
        for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
            for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                //kernels.DIRECT_COMPUTATION_MUTUAL_SOFT(particles[idxTarget], particles[idxOther]);
                kernels.directInteractionMutual(&particles[idxTarget], &particles[idxOther]);
            }
        }

        Print("Compute Diff...");
        FReal potentialDiff = 0;
        FReal fx = 0, fy = 0, fz = 0;

        { // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();

            do{
                typename ContainerClass::BasicIterator leafIter(*octreeIterator.getCurrentListTargets());

                while( leafIter.hasNotFinished() ){
                    const ParticleClass& other = particles[leafIter.data().getIndex()];

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

    // set test
    void SetTests(){
        AddTest(&TestFmbDirect::TestDirect,"Test Simu and with direct");
    }
};


// You must do this
TestClass(TestFmbDirect)


// [--END--]
