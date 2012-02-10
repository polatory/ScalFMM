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

#include "../Src/Kernels/FSphericalCell.hpp"
#include "../Src/Kernels/FSphericalKernel.hpp"
#include "../Src/Kernels/FSphericalParticle.hpp"
#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FTreeIO.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"


/*
  This test compare a previous FMM result with a simulation.
  */

/** The result of a previous simulation has been saved and will be oponed */
class ParticleSerial : public FSphericalParticle, public FTreeIO::FAbstractSerial {
public:
    void write(std::ofstream*const stream) const{
        save(stream, getPosition().getX());
        save(stream, getPosition().getY());
        save(stream, getPosition().getZ());
        save(stream, getForces().getX());
        save(stream, getForces().getY());
        save(stream, getForces().getZ());
        save(stream, getPotential());
        save(stream, getPhysicalValue());
    }

    void read(std::ifstream*const stream) {
        const FReal posX = restore<FReal>(stream);
        const FReal posY = restore<FReal>(stream);
        const FReal posZ = restore<FReal>(stream);
        setPosition(posX, posY, posZ);

        const FReal forceX = restore<FReal>(stream);
        const FReal forceY = restore<FReal>(stream);
        const FReal forceZ = restore<FReal>(stream);
        setForces(forceX, forceY, forceZ);

        setPotential(restore<FReal>(stream));
        setPhysicalValue(restore<FReal>(stream));
    }
};

/** The result of a previous simulation has been saved and will be oponed */
class ComputeCellSerial : public FSphericalCell, public FTreeIO::FAbstractSerial {
public:
    void write(std::ofstream*const stream) const{
        saveArray(stream, FSphericalCell::getMultipole(), FSphericalCell::GetPoleSize());
        saveArray(stream, FSphericalCell::getLocal(), FSphericalCell::GetLocalSize());
    }

    void read(std::ifstream*const stream){
        restoreArray(stream, FSphericalCell::getMultipole(), FSphericalCell::GetPoleSize());
        restoreArray(stream, FSphericalCell::getLocal(), FSphericalCell::GetLocalSize());
    }
};

typedef ParticleSerial       ParticleClass;
typedef ComputeCellSerial             CellClass;
typedef FVector<ParticleClass>  ContainerClass;

typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >          KernelClass;

typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;

typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

class TestFmb : public FUTester<TestFmb> {

    void TestTree(){
        // Warning in make test the exec dir it Build/UTests
        const char* const DataFile = "../../Data/utestFmb.data";
        const char* const ParticleFile = "../../Data/utestFmb.bin.fma";

        const int NbLevels      = 5;
        const int SizeSubLevels = 3;
        const int DevP = 12;

        // Load the particles file
        FFmaBinLoader<ParticleClass> loader(ParticleFile);
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            assert(false);
            return;
        }

        // Create octree
        FSphericalCell::Init(DevP);
        OctreeClass testTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        loader.fillTree(testTree);

        // Run simulation
        KernelClass kernels(DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClass algo(&testTree,&kernels);
        algo.execute();

        // If needed save the result
        // FTreeIO::Save<OctreeClass, CellClass, ParticleClass, FTreeIO::Serializer<CellClass, ParticleClass> >(DataFile, testTree);

        // Load previous result
        OctreeClass goodTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        FTreeIO::Load<OctreeClass, CellClass, ParticleClass, FTreeIO::Serializer<CellClass, ParticleClass> >(DataFile, goodTree);

        // Compare the two simulations
        Print("Check the particles...");
        { // Check that each particle has been summed with all other
            typename OctreeClass::Iterator testOctreeIterator(&testTree);
            typename OctreeClass::Iterator goodOctreeIterator(&goodTree);

            testOctreeIterator.gotoBottomLeft();
            goodOctreeIterator.gotoBottomLeft();

            do{
                if(testOctreeIterator.getCurrentGlobalIndex() != goodOctreeIterator.getCurrentGlobalIndex()){
                    assert(false);
                    break;
                }

                if(testOctreeIterator.getCurrentListSrc()->getSize() != goodOctreeIterator.getCurrentListSrc()->getSize()){
                    assert(false);
                    break;
                }

                typename ContainerClass::BasicIterator goodIter(*goodOctreeIterator.getCurrentListTargets());
                typename ContainerClass::BasicIterator testIter(*testOctreeIterator.getCurrentListTargets());

                while( goodIter.hasNotFinished() ){
                    assert( memcmp(&goodIter.data(), &testIter.data(), sizeof(ParticleClass)) == 0);

                    goodIter.gotoNext();
                    testIter.gotoNext();
                }


                if(!testOctreeIterator.moveRight()){
                    if(goodOctreeIterator.moveRight()){
                        assert(false);
                    }
                    break;
                }
                if(!goodOctreeIterator.moveRight()){
                    assert(false);
                    break;
                }

            } while(true);
        }
        Print("Check the leaves...");
        { // Ceck if there is number of NbPart summed at level 1
            typename OctreeClass::Iterator testOctreeIterator(&testTree);
            typename OctreeClass::Iterator goodOctreeIterator(&goodTree);

            testOctreeIterator.gotoBottomLeft();
            goodOctreeIterator.gotoBottomLeft();

            for(int idxLevel = NbLevels - 1 ; idxLevel > 1 ; --idxLevel ){
                do{
                    if(testOctreeIterator.getCurrentGlobalIndex() != goodOctreeIterator.getCurrentGlobalIndex()){
                        assert(false);
                        break;
                    }

                    assert( memcmp(testOctreeIterator.getCurrentCell()->getLocal(),
                                   goodOctreeIterator.getCurrentCell()->getLocal(), CellClass::GetLocalSize() * sizeof(FComplexe)) == 0);

                    assert( memcmp(testOctreeIterator.getCurrentCell()->getMultipole(),
                                   goodOctreeIterator.getCurrentCell()->getMultipole(),CellClass::GetPoleSize() * sizeof(FComplexe)) == 0);

                    if(!testOctreeIterator.moveRight()){
                        if(goodOctreeIterator.moveRight()){
                            assert(false);
                        }
                        break;
                    }
                    if(!goodOctreeIterator.moveRight()){
                        assert(false);
                        break;
                    }

                } while(true);

                testOctreeIterator.moveUp();
                testOctreeIterator.gotoLeft();

                goodOctreeIterator.moveUp();
                goodOctreeIterator.gotoLeft();
            }
        }
        Print("Over...");
    }


    // set test
    void SetTests(){
        AddTest(&TestFmb::TestTree,"Test Simu and compare tree");
    }
};



// You must do this
TestClass(TestFmb)



