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

// ==== CMAKE =====
// @FUSE_STARPU
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>


#include "../../Src/Starpu/FCommon.hpp"

#include "../../Src/Files/FRandomLoader.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Components/FBasicCell.hpp"
#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"
#include "../../Src/Starpu/FCommon.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#include "../Cuda/eventTimer.h"

// Copy the function to group particles in a memory block

static GLeavesGroup* GLeavesGroupFromMemory(char* memory){
    return (GLeavesGroup*) memory;
}

static GLeaf* GLeavesFromMerory(char* memory){
    return (GLeaf*)(memory + sizeof(GLeavesGroup));
}

static int GOffsetTo8(const int address){
    return ((address + 7) & ~7) - address;
}

static int GHeaderSize(const int nbLeaves){
    return sizeof(GLeavesGroup) + (nbLeaves * sizeof(GLeaf));
}

static GParticle* GParticlesFromMemory(char* memory){
    GLeavesGroup* group = GLeavesGroupFromMemory(memory);
    const int header = GHeaderSize(group->nbLeaves);
    const int offset = GOffsetTo8(header);
    return (GParticle*)  (memory + header + offset);
}

static int GSizeMemory(const int nbLeaves, const int nbParticles){
    const int header = GHeaderSize(nbLeaves);
    const int offset = GOffsetTo8(header);
    return header + offset + (sizeof(GParticle) * nbParticles);
}

static int GGetNbParticles(char* memory, const int memorySize){
    const char*const particlesAddress = (char*)GParticlesFromMemory(memory);
    return ((memorySize + memory) - particlesAddress)/sizeof(GParticle);
}


// Simply create particles and try the kernels
int main(int argc, char* argv[]){
    // parameters to configure the run
    const unsigned int NbPart     = FParameters::getValue(argc, argv, "-nb", 100000);
    const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
    const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);
    const int heightMinusOne = TreeHeight - 1;

    const int GpuFlop = 21;
    const int CpuFlops = GpuFlop + 5 + 14; // mutual + sqrt on X86

    // typedefs for STARPU
    typedef GParticle ParticleClass;
    typedef FVector<ParticleClass> ContainerClass;
    typedef FSimpleLeaf<ParticleClass,ContainerClass> LeafClass;
    typedef FBasicCell CellClass;
    typedef FOctree<ParticleClass,CellClass,ContainerClass,LeafClass> OctreeClass;
    typedef FSphericalKernel<ParticleClass, FSphericalCell, ContainerClass >     KernelClass;

    // What we do //////////////////////////////////////////////////////
    std::cout << ">> Will generate " << NbPart << " particles randomly.\n";
    std::cout << ">> On a tree of height " << TreeHeight << "\n";
    std::cout << ">> So there is a maximum of " << ((1 << (3 * (TreeHeight - 1))) - 1) << " leaves\n";
    std::cout << ">> Theorical averages particles per leaf " << FReal(NbPart)/FReal((1 << (3 * (TreeHeight - 1))) - 1) << "\n";


    // Insert particles in tree  ///////////////////////////////////////
    FRandomLoader<ParticleClass> loader(NbPart);
    OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());
    {
        ParticleClass part;
        part.setPhysicalValue(FReal(0.1));
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(part);
            tree.insert(part);
        }
    }
    std::cout << "\tDone.\n";

    long long counterNbInteraction = 0;
    // Compute nb interactions /////////////////////////////////////////
    {
        std::cout << "Find the number of interactions between particles...\n";
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        // There is a maximum of 26 neighbors
        ContainerClass* neighbors[27];
        // for each leafs
        do{
            // need the current particles and neighbors particles
            tree.getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),heightMinusOne);
            // interaction in current leaf
            const int nbPartInLeaf =  octreeIterator.getCurrentListSrc()->getSize();;
            counterNbInteraction += nbPartInLeaf * (nbPartInLeaf-1);

            for(int idxInter = 0 ; idxInter < 27 ; ++idxInter){
                if( neighbors[idxInter] ){
                    // interactions betwen leaf and neighbors
                    counterNbInteraction += nbPartInLeaf * neighbors[idxInter]->getSize();
                }
            }
        } while(octreeIterator.moveRight());

        std::cout << "\tDone.\n";
        std::cout << "\tThere is " << counterNbInteraction << " interactions.\n";
    }

    // Put cells in array /////////////////////////////////////////
    typename OctreeClass::Iterator* leavesArray = 0;
    int nbLeaves = 0;
    {
        std::cout << "Prepare data for CPUs...\n";
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();

        do{
            ++nbLeaves;
        }while(octreeIterator.moveRight());

        leavesArray = new typename OctreeClass::Iterator[nbLeaves];
        octreeIterator.gotoLeft();
        int iterLeaf = 0;
        do{
            leavesArray[iterLeaf++] = octreeIterator;
        }while(octreeIterator.moveRight());
    }

    // CPU interactions /////////////////////////////////////////
    {
        KernelClass kernel(5, TreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

        std::cout << "Compute on the cpu...\n";
        // init timer
        FTic timer;

        // There is a maximum of 26 neighbors
        ContainerClass* neighbors[27];
        // for each leafs
        for(int idxLeaf = 0 ; idxLeaf < nbLeaves ; ++idxLeaf){
            // need the current particles and neighbors particles
            const int counter = tree.getLeafsNeighbors(neighbors, leavesArray[idxLeaf].getCurrentGlobalCoordinate(),heightMinusOne);
            kernel.P2P(leavesArray[idxLeaf].getCurrentGlobalCoordinate(),leavesArray[idxLeaf].getCurrentListTargets(),
                         leavesArray[idxLeaf].getCurrentListSrc(), neighbors, counter);
        }

        std::cout << "\tDone in " << timer.tacAndElapsed() << "\n";
        std::cout << "\tOne operation cost " << CpuFlops << "flop so total is " << (counterNbInteraction/2) * CpuFlops << "flop\n";
        std::cout << "\tFlop/s = " << double((counterNbInteraction/2) * CpuFlops)/timer.elapsed() << std::endl;
    }
    delete[] leavesArray;
    leavesArray = 0;

    // Compute on gpu and check result /////////////////////////////
    {
        // init timer
        FTic timer;
        std::cout << "Copy data to GPUs format...\n";
        nbLeaves = 0;
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        // count nb leaves
        do{
            ++nbLeaves;
        } while(octreeIterator.moveRight());

        // copy to group
        const int MemorySize = GSizeMemory(nbLeaves, NbPart);
        char*const gpusParticles = new char[MemorySize];
        {
            // Copy data
            GLeavesGroup* leaves = GLeavesGroupFromMemory( gpusParticles );
            leaves->begin = 0;
            leaves->end = (1 << (3 * (TreeHeight - 1))) - 1;
            leaves->nbLeaves = nbLeaves;

            GLeaf* leaf = GLeavesFromMerory(gpusParticles);
            GParticle* particles = GParticlesFromMemory(gpusParticles);

            int particlesPosition = 0;
            octreeIterator.gotoLeft();

            do {
                memcpy(leaf->coordinate, &octreeIterator.getCurrentGlobalCoordinate(), sizeof(int[3]));
                leaf->mindex = octreeIterator.getCurrentGlobalIndex();
                leaf->nbParticles = octreeIterator.getCurrentListSrc()->getSize();
                leaf->positionInArray = particlesPosition;

                FMemUtils::copyall( particles, octreeIterator.getCurrentListSrc()->data(), leaf->nbParticles);
                for(int idxPart = 0 ; idxPart < leaf->nbParticles ; ++idxPart){
                    particles[idxPart].setPotential(FReal(0));
                    particles[idxPart].getForces().setPosition(FReal(0),FReal(0),FReal(0));
                }

                particlesPosition += leaf->nbParticles;
                particles += leaf->nbParticles;
                leaf += 1;
            } while(octreeIterator.moveRight());
        }
        std::cout << "\tDone in " << timer.tacAndElapsed() << "\n";

        // copy data to gpu
        std::cout << "Copy to Gpu...\n";
        timer.tic();
        char* cuda_particles = 0;
        cudaMalloc( (void **)&cuda_particles, MemorySize);
        cudaMemcpy(cuda_particles, gpusParticles, MemorySize, cudaMemcpyHostToDevice);

        const int EmptyMemorySize = GSizeMemory(0, 0);
        char*const emptyGpusParticles = new char[EmptyMemorySize];

        // Create a fack empty group
        GLeavesGroup* emptyleaves = GLeavesGroupFromMemory( emptyGpusParticles );
        emptyleaves->begin = 0;
        emptyleaves->end = 0;
        emptyleaves->nbLeaves = 0;
        char* cuda_emptyparticles = 0;
        cudaMalloc( (void **)&cuda_emptyparticles, EmptyMemorySize);
        cudaMemcpy(cuda_emptyparticles, emptyleaves, EmptyMemorySize, cudaMemcpyHostToDevice);

        std::cout << "\tDone in " << timer.tacAndElapsed() << "\n";

        // compute
        std::cout << "Compute on Gpu...\n";
        timer.tic();

        eventTimerType t;
        initEventTimer(&t);
        startEventTimer(&t);
        /* call the gpu func */
        perf_direct_cuda_func(cuda_particles,
                              cuda_emptyparticles,
                              TreeHeight,
                              nbLeaves);
        stopEventTimer(&t);
        std::cout << "\tDone in " <<  getEventTimer(&t) << " (event timer)\n"; // unit is second
        finalizeEventTimer(&t);

        std::cout << "\tDone in " << timer.tacAndElapsed() << "\n";
        std::cout << "\tOne operation cost " << GpuFlop << "flop so total is " << counterNbInteraction * GpuFlop << "flop\n";
        std::cout << "\tFlop/s = " << double(counterNbInteraction * GpuFlop)/timer.elapsed() << std::endl;

        // Retrieve results
        cudaMemcpy(gpusParticles, cuda_particles, MemorySize, cudaMemcpyDeviceToHost);
        cudaFree(cuda_particles);

        // Check results
        std::cout << "Check result...\n";
        FMath::FAccurater accuracy;
        {
            // Copy data
            const GLeaf* leaf = GLeavesFromMerory(gpusParticles);
            const GParticle*const particles = GParticlesFromMemory(gpusParticles);

            octreeIterator.gotoLeft();

            do {
                const ContainerClass& validParticles = *octreeIterator.getCurrentListSrc();

                for(int idxPart = 0 ; idxPart < leaf->nbParticles ; ++idxPart){
                    accuracy.add( validParticles[idxPart].getPotential(), particles[leaf->positionInArray + idxPart].getPotential());
                    accuracy.add( validParticles[idxPart].getForces().getX(), particles[leaf->positionInArray + idxPart].getForces().getX());
                    accuracy.add( validParticles[idxPart].getForces().getY(), particles[leaf->positionInArray + idxPart].getForces().getY());
                    accuracy.add( validParticles[idxPart].getForces().getZ(), particles[leaf->positionInArray + idxPart].getForces().getZ());

                    /*std::cout << "Pos " << validParticles[idxPart].getPosition() << "\n";
                    std::cout << "GPos " << particles[leaf->positionInArray + idxPart].getPosition() << "\n";
                    std::cout << "\t Potential: " << validParticles[idxPart].getPotential() << " / "
                              << particles[leaf->positionInArray + idxPart].getPotential() << "\n";
                    std::cout << "\t Physical value: " << validParticles[idxPart].getPhysicalValue() << " / "
                              << particles[leaf->positionInArray + idxPart].getPhysicalValue() << "\n";*/
                }

                leaf += 1;
            } while(octreeIterator.moveRight());
        }
        std::cout << "\tDone, inf norm = " << accuracy.getInfNorm() << " l2 norm = " << accuracy.getL2Norm() << "\n";

        delete[] gpusParticles;
    }

    return 0;
}



