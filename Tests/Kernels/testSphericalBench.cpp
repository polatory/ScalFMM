
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

#include <iostream>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmTask.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Files/FFmaLoader.hpp"




void printDetails(const char name[], const int NbLevels, const int DevP, const FMath::FAccurater diffs[]){
    printf("Print %s:\n", name);
    printf("H\\P\t|");
    for(int idxP = 2 ; idxP < DevP ; ++idxP) printf("%12d             |",idxP);
    printf("\n");
    for(int idxH = 2 ; idxH < NbLevels ; ++idxH){
        printf("[%d]\t|", idxH);
        for(int idxP = 2 ; idxP < DevP ; ++idxP){
            printf("%10e/%10e|", diffs[idxH*DevP+idxP].getL2Norm(),
                   diffs[idxH*DevP+idxP].getInfNorm());
        }
        printf("\n");
    }
}


/** This program show an example of use of
  * the fmm basic algo
  * it also check that each particles is little or longer
  * related that each other
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

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef IndexedParticle                ParticleClass;
    typedef FSphericalCell                 CellClass;
    typedef FVector<ParticleClass>         ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >     KernelClass;

    typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    typedef FFmmAlgorithmThread<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassThread;
    typedef FFmmAlgorithmTask<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassTask;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    std::cout << ">> You can pass -sequential or -task (thread by default).\n";
    //////////////////////////////////////////////////////////////
    const int DevP = FParameters::getValue(argc,argv,"-p", 8);
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 5);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");

    FMath::FAccurater* allPotentialDiff= new FMath::FAccurater[NbLevels*DevP];
    FMath::FAccurater* allXDiff= new FMath::FAccurater[NbLevels*DevP];
    FMath::FAccurater* allYDiff= new FMath::FAccurater[NbLevels*DevP];
    FMath::FAccurater* allZDiff= new FMath::FAccurater[NbLevels*DevP];


    ParticleClass* particles = 0;
    {
        FFmaLoader<ParticleClass> loader(filename);
        particles = new ParticleClass[loader.getNumberOfParticles()];
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(particles[idxPart]);
            particles[idxPart].setIndex( idxPart );
        }

        // Compute direct
        printf("Compute direct!\n");
        KernelClass kernels(DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
            for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                kernels.directInteractionMutual(&particles[idxTarget], &particles[idxOther]);
            }
        }
    }

    for(int idxH = 2 ; idxH < NbLevels ; ++idxH){
        for(int idxP = 2 ; idxP < DevP ; ++idxP){
            std::cout << "Running!!! H = " << idxH << " P = " << idxP << std::endl;

            std::cout << "Opening : " << filename << "\n";
            FFmaLoader<ParticleClass> loader(filename);
            if(!loader.isOpen()){
                std::cout << "Loader Error, " << filename << " is missing\n";
                return 1;
            }

            // -----------------------------------------------------
            CellClass::Init(idxP);
            OctreeClass tree(idxH, SizeSubLevels>=idxH?idxH-1:SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

            // -----------------------------------------------------

            std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
            std::cout << "\tHeight : " << idxH << " \t sub-height : " << SizeSubLevels << std::endl;
            counter.tic();

            for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                ParticleClass part;
                loader.fillParticle(part);
                part.setIndex( idxPart );
                tree.insert(part);
            }

            counter.tac();
            std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

            // -----------------------------------------------------

            std::cout << "Create kernel ..." << std::endl;
            counter.tic();

            KernelClass kernels(idxP, idxH, loader.getBoxWidth(), loader.getCenterOfBox());

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

            // -----------------------------------------------------

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

            // Print for information
            printf("Potential diff is = %e \t %e\n",potentialDiff.getL2Norm(),potentialDiff.getInfNorm());
            printf("Fx diff is = %e \t %e\n",fx.getL2Norm(),fx.getInfNorm());
            printf("Fy diff is = %e \t %e\n",fy.getL2Norm(),fy.getInfNorm());
            printf("Fz diff is = %e \t %e\n",fz.getL2Norm(),fz.getInfNorm());

            allPotentialDiff[idxH*DevP+idxP] = potentialDiff;
            allXDiff[idxH*DevP+idxP] = fx;
            allYDiff[idxH*DevP+idxP] = fy;
            allZDiff[idxH*DevP+idxP] = fz;
        }
    }

    delete[] particles;

    printDetails("potential", NbLevels, DevP, allPotentialDiff);
    printDetails("force x", NbLevels, DevP, allXDiff);
    printDetails("force y", NbLevels, DevP, allYDiff);
    printDetails("force z", NbLevels, DevP, allZDiff);

    delete[] allPotentialDiff;
    delete[] allXDiff;
    delete[] allYDiff;
    delete[] allZDiff;

    return 0;
}



