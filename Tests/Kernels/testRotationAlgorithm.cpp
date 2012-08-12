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

#include <limits>
#include <iostream>


#include "../../Src/Components/FSimpleLeaf.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Containers/FOctree.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"
#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"

#include "../../Src/Kernels/Rotation/FRotationKernel.hpp"
#include "../../Src/Kernels/Rotation/FRotationCell.hpp"

#include "../../Src/Utils/FMath.hpp"
#include "../../Src/Utils/FMemUtils.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Files/FFmaLoader.hpp"



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


int main(int argc, char** argv){
    static const int P = 8;

    typedef IndexedParticle                ParticleClass;
    typedef FRotationCell<P>               CellClass;
    typedef FVector<ParticleClass>         ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FRotationKernel<ParticleClass, CellClass, ContainerClass , P>   KernelClass;

    typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    typedef FFmmAlgorithmThread<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassThread;
    typedef FFmmAlgorithmTask<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassTask;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    std::cout << ">> You can pass -sequential or -task (thread by default).\n";
    //////////////////////////////////////////////////////////////
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");

    std::cout << "Opening : " << filename << "\n";

    FFmaLoader<ParticleClass> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    //loader.fillTree(tree);
    ParticleClass* const particles = new ParticleClass[loader.getNumberOfParticles()];
    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(particles[idxPart]);
        particles[idxPart].setIndex( idxPart );
        tree.insert(particles[idxPart]);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Create kernel ..." << std::endl;
    counter.tic();

    KernelClass kernels(NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    counter.tac();
    std::cout << "Done  " << " in " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;

    if( FParameters::findParameter(argc,argv,"-sequential") != FParameters::NotFound){
        FmmClass algo(&tree,&kernels);
        counter.tic();
        algo.execute();
    }
    else if( FParameters::findParameter(argc,argv,"-task") != FParameters::NotFound){
        FmmClassTask algo(&tree,&kernels);
        counter.tic();
        algo.execute();
    }
    else {
        FmmClassThread algo(&tree,&kernels);
        counter.tic();
        algo.execute();
    }

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    { // get sum forces&potential
        FReal potential = 0;
        FPoint forces;
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){
                potential += iter.data().getPotential() * iter.data().getPhysicalValue();
                forces += iter.data().getForces();

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());

        std::cout << "Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
        std::cout << "Potential = " << potential << std::endl;
    }

    // -----------------------------------------------------

    {
        std::cout << "Compute direct interaction for all\n";
        for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
            for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                kernels.particlesMutualInteraction(&particles[idxTarget], &particles[idxOther]);
            }
        }

        std::cout << "Compute Diff...\n";
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
//std::cout << leafIter.data().getIndex() << " Direct Potential = " << other.getPotential() << " Fmm Potential = " << leafIter.data().getPotential() << std::endl; // Remove Me
                    fx.add(other.getForces().getX(),leafIter.data().getForces().getX());

                    fy.add(other.getForces().getY(),leafIter.data().getForces().getY());

                    fz.add(other.getForces().getZ(),leafIter.data().getForces().getZ());

                    leafIter.gotoNext();
                }
            } while(octreeIterator.moveRight());
        }

        // Print for information
        std::cout << "Potential diff is = " << potentialDiff.getL2Norm() << " " << potentialDiff.getInfNorm() << std::endl;
        std::cout << "Fx diff is = " << fx.getL2Norm() << " " << fx.getInfNorm() << std::endl;
        std::cout << "Fy diff is = " << fy.getL2Norm() << " " << fy.getInfNorm() << std::endl;
        std::cout << "Fz diff is = " << fz.getL2Norm() << " " << fz.getInfNorm() << std::endl;
    }

    delete[] particles;

    {
        typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass>   ShKernelClass;
        typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, ShKernelClass, LeafClass > ShFmmClass;

        FFmaLoader<ParticleClass> loaderSh(filename);
        OctreeClass treeSh(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            ParticleClass particle;
            loaderSh.fillParticle(particle);
            treeSh.insert(particle);
        }

        ShKernelClass kernelsSh(P, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        ShFmmClass algoSh(&treeSh,&kernelsSh);
        algoSh.execute();

        {
            const int OctreeHeight = tree.getHeight();
            typename OctreeClass::Iterator octreeIterator(&treeSh);
            octreeIterator.gotoBottomLeft();

            typename OctreeClass::Iterator octreeIteratorValide(&tree);
            octreeIteratorValide.gotoBottomLeft();

            for(int level = OctreeHeight - 1 ; level > 1 ; --level){
                FMath::FAccurater pole, local;

                do {
                    for(int idx = 0 ; idx < ((P+2)*(P+1))/2 ; ++idx){
                        pole.add(octreeIteratorValide.getCurrentCell()->getMultipole()[idx].getReal(),
                                 octreeIterator.getCurrentCell()->getMultipole()[idx].getReal());
                        pole.add(octreeIteratorValide.getCurrentCell()->getMultipole()[idx].getImag(),
                                 octreeIterator.getCurrentCell()->getMultipole()[idx].getImag());
                    }
                    for(int idx = 0 ; idx < ((P+2)*(P+1))/2 ; ++idx){
                        local.add(octreeIteratorValide.getCurrentCell()->getLocal()[idx].getReal(),
                                  octreeIterator.getCurrentCell()->getLocal()[idx].getReal());
                        local.add(octreeIteratorValide.getCurrentCell()->getLocal()[idx].getImag(),
                                  octreeIterator.getCurrentCell()->getLocal()[idx].getImag());
                    }
                } while(octreeIterator.moveRight() && octreeIteratorValide.moveRight());

                octreeIterator.moveUp();
                octreeIterator.gotoLeft();

                octreeIteratorValide.moveUp();
                octreeIteratorValide.gotoLeft();

                printf("At level %d:", level);
                printf(">> pole %f, %f\n", pole.getL2Norm(), pole.getInfNorm() );
                printf(">> local %f, %f\n", local.getL2Norm(), local.getInfNorm() );
            }
        }
        {
            // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(&treeSh);
            octreeIterator.gotoBottomLeft();

            typename OctreeClass::Iterator octreeIteratorValide(&tree);
            octreeIteratorValide.gotoBottomLeft();

            FMath::FAccurater potential, forces;

            do {
                typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());
                typename ContainerClass::BasicIterator iterValide(*octreeIteratorValide.getCurrentListTargets());

                while( iter.hasNotFinished() ){

                   potential.add( iterValide.data().getPotential() , iter.data().getPotential());
                   forces.add(iterValide.data().getForces().getX(),iter.data().getForces().getX());
                   forces.add(iterValide.data().getForces().getY(),iter.data().getForces().getY());
                   forces.add(iterValide.data().getForces().getZ(),iter.data().getForces().getZ());

                    iter.gotoNext();
                    iterValide.gotoNext();
                }

            } while(octreeIterator.moveRight() && octreeIteratorValide.moveRight());
        }
    }

    return 0;
}

