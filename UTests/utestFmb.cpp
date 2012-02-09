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

#include "../Src/Kernels/FSphericalCell.hpp"
#include "../Src/Kernels/FSphericalKernel.hpp"

#include "../Src/Fmb/FFmbKernels.hpp"
#include "../Src/Fmb/FFmbComponents.hpp"

#include "../Src/Files/FFmaBinLoader.hpp"
#include "../Src/Files/FTreeIO.hpp"

#include "../Src/Core/FFmmAlgorithm.hpp"

#include "FUTester.hpp"


/*
  This test compare a previous FMM result with a simulation.
  */

/** The result of a previous simulation has been saved and will be oponed */
class FmbParticleSerial : public FmbParticle, public FTreeIO::FAbstractSerial {
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
        saveArray(stream, FSphericalCell::getMultipole(), FSphericalCell::ExpP);
        saveArray(stream, FSphericalCell::getLocal(), FSphericalCell::ExpP);
    }

    void read(std::ifstream*const stream){
        restoreArray(stream, FSphericalCell::getMultipole(), FSphericalCell::ExpP);
        restoreArray(stream, FSphericalCell::getLocal(), FSphericalCell::ExpP);
    }
};

typedef FmbParticleSerial       ParticleClass;
typedef ComputeCellSerial             CellClass;
typedef FVector<ParticleClass>  ContainerClass;

//typedef FFmbKernels<ParticleClass, CellClass, ContainerClass >          KernelClass;
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

        FSphericalCell::Init(DevP);

        // Load the particles file
        FFmaBinLoader<ParticleClass> loader(ParticleFile);
        if(!loader.isOpen()){
            Print("Cannot open particles file.");
            assert(false);
            return;
        }

        // Create octree
        OctreeClass testTree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        {
            ParticleClass particleToFill;
            for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                loader.fillParticle(particleToFill);
                testTree.insert(particleToFill);
            }
        }

        // Run simulation
        KernelClass kernels(DevP,NbLevels,loader.getBoxWidth());
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
                                   goodOctreeIterator.getCurrentCell()->getLocal(), CellClass::GetExp() * sizeof(FComplexe)) == 0);

                    assert( memcmp(testOctreeIterator.getCurrentCell()->getMultipole(),
                                   goodOctreeIterator.getCurrentCell()->getMultipole(),CellClass::GetExp() * sizeof(FComplexe)) == 0);

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



