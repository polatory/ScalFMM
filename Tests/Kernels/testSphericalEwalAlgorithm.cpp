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

#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cstdlib>

#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithmPeriodic.hpp"
#include "../../Src/Core/FFmmAlgorithm.hpp"

#include "../../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../../Src/Kernels/Spherical/FSphericalCell.hpp"

#include "../../Src/Files/FEwalLoader.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

#include "../../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"

/** Ewal particle is used in the gadget program
  * here we try to make the same simulation
  */

struct EwalParticle {
    FPoint position;
    FReal forces[3];
    FReal physicalValue;
    FReal potential;
    int index;
};

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FSphericalCell          CellClass;
    typedef FP2PParticleContainerIndexed  ContainerClass;

    typedef FSimpleLeaf< ContainerClass >                     LeafClass;
    typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel< CellClass, ContainerClass >   KernelClass;

    typedef FFmmAlgorithmPeriodic<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    typedef FFmmAlgorithm<OctreeClass,  CellClass, ContainerClass, KernelClass, LeafClass > FmmClassNoPer;

    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    std::cout << ">> options are -h H -sh SH -P p -per PER -f FILE -noper -verbose -gen -direct\n";
    std::cout << ">> Recommanded files : ../Data/EwalTest_Periodic.run ../Data/EwalTest_NoPeriodic.run\n";
    //////////////////////////////////////////////////////////////

    const int NbLevels      = FParameters::getValue(argc,argv,"-h", 4);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 2);
    const int DevP          = FParameters::getValue(argc,argv,"-P", 9);
    const int PeriodicDeep      = FParameters::getValue(argc,argv,"-per", 2);
    const char* const filename = FParameters::getStr(argc,argv,"-f", "../Data/EwalTest_Periodic.run");
    // recommenda

    FTic counter;
    const FReal coeff_MD  = FReal(138935.4835 / 418.4);
    const FReal coeff_MD1 = FReal(138935.4835);

    // -----------------------------------------------------

    std::cout << "Opening : " << filename << "\n";
    //FEwalLoader loader(filename);
    FEwalBinLoader loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    CellClass::Init(DevP);
    OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tWidth : " << loader.getBoxWidth() << " \t center x : " << loader.getCenterOfBox().getX()
              << " y : " << loader.getCenterOfBox().getY() << " z : " << loader.getCenterOfBox().getZ() << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;

    counter.tic();

    EwalParticle * const particles = new EwalParticle[loader.getNumberOfParticles()];
    memset(particles, 0, sizeof(EwalParticle) * loader.getNumberOfParticles());

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(&particles[idxPart].position, particles[idxPart].forces,
                            &particles[idxPart].physicalValue,&particles[idxPart].index);
        // reset forces and insert in the tree
        tree.insert(particles[idxPart].position, idxPart, particles[idxPart].physicalValue);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Create kernel & run simu ..." << std::endl;
    counter.tic();

    FTreeCoordinate min, max;

    if( FParameters::existParameter(argc, argv, "-noper") ){
        KernelClass kernels( DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        FmmClassNoPer algo(&tree,&kernels);
        algo.execute();
    }
    else{
        FmmClass algo(&tree,PeriodicDeep);
        std::cout << "The simulation box is repeated " << algo.repetitions() << " in X/Y/Z" << std::endl;
        KernelClass kernels( DevP, algo.extendedTreeHeight(), algo.extendedBoxWidth(),algo.extendedBoxCenter());
        algo.setKernel(&kernels);
        algo.execute();
        algo.repetitionsIntervals(&min, &max);
    }

    counter.tac();

    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    EwalParticle* particlesDirect = 0;

    // Do direct
    if(FParameters::existParameter(argc, argv, "-direct")){
        printf("Compute direct:\n");
        printf("Box [%d;%d][%d;%d][%d;%d]\n", min.getX(), max.getX(), min.getY(),
               max.getY(), min.getZ(), max.getZ());

        particlesDirect = new EwalParticle[loader.getNumberOfParticles()];

        FReal dpotential = 0.0;
        FMath::FAccurater dfx, dfy, dfz;

        for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
            EwalParticle part = particles[idxTarget];
            part.forces[0] = part.forces[1] = part.forces[2] = 0.0;
            part.potential = 0.0;
            // compute with all other except itself
            for(int idxOther = 0; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                if( idxOther != idxTarget ){
                    FP2P::NonMutualParticles(
                                particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                particles[idxOther].position.getZ(),particles[idxOther].physicalValue,
                                part.position.getX(), part.position.getY(),
                                part.position.getZ(),part.physicalValue,
                                &part.forces[0],&part.forces[1],
                                &part.forces[2],&part.potential);
                }
            }
            for(int idxX = min.getX() ; idxX <= max.getX() ; ++idxX){
                for(int idxY = min.getY() ; idxY <= max.getY() ; ++idxY){
                    for(int idxZ = min.getZ() ; idxZ <= max.getZ() ; ++idxZ){
                        if(idxX == 0 && idxY == 0 && idxZ == 0) continue;

                        const FPoint offset(loader.getBoxWidth() * FReal(idxX),
                                            loader.getBoxWidth() * FReal(idxY),
                                            loader.getBoxWidth() * FReal(idxZ));

                        for(int idxSource = 0 ; idxSource < loader.getNumberOfParticles() ; ++idxSource){
                            EwalParticle source = particles[idxSource];
                            source.position += offset;
                            FP2P::NonMutualParticles(
                                        source.position.getX(), source.position.getY(),
                                        source.position.getZ(),source.physicalValue,
                                        part.position.getX(), part.position.getY(),
                                                  part.position.getZ(),part.physicalValue,
                                                  &part.forces[0],&part.forces[1],
                                                  &part.forces[2],&part.potential
                                            );
                        }
                    }
                }
            }

            dfx.add(particles[idxTarget].forces[0],-part.forces[0]*coeff_MD1);
            dfy.add(particles[idxTarget].forces[1],-part.forces[1]*coeff_MD1);
            dfz.add(particles[idxTarget].forces[2],-part.forces[2]*coeff_MD1);

            dpotential += part.potential * part.physicalValue;

            particlesDirect[idxTarget] = part;
        }

        printf("Difference between direct and poly:\n");
        printf("Fx diff is = \n");
        printf("%e\n",dfx.getL2Norm());
        printf("%e\n",dfx.getInfNorm());
        printf("Fy diff is = \n");
        printf("%e\n",dfy.getL2Norm());
        printf("%e\n",dfy.getInfNorm());
        printf("Fz diff is = \n");
        printf("%e\n",dfz.getL2Norm());
        printf("%e\n",dfz.getInfNorm());
        printf("Direct Potential= %e\n",dpotential*coeff_MD/2);
    }

    // -----------------------------------------------------
    {
        FReal  potential = 0;
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;


        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
            const FReal*const positionsX = leaf->getTargets()->getPositions()[0];
            const FReal*const positionsY = leaf->getTargets()->getPositions()[1];
            const FReal*const positionsZ = leaf->getTargets()->getPositions()[2];
            const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
            const FVector<int>& indexes = leaf->getTargets()->getIndexes();

            for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                const int indexPartOrig = indexes[idxPart];
                potentialDiff.add(particles[indexPartOrig].potential,potentials[idxPart]);
                fx.add(particles[indexPartOrig].forces[0],forcesX[idxPart]);
                fy.add(particles[indexPartOrig].forces[1],forcesY[idxPart]);
                fz.add(particles[indexPartOrig].forces[2],forcesZ[idxPart]);

                if(FParameters::existParameter(argc, argv, "-verbose")){
                    std::cout << ">> index " << particles[indexPartOrig].index << std::endl;
                    std::cout << "Good x " << particles[indexPartOrig].position.getX() << " y " << particles[indexPartOrig].position.getY() << " z " << particles[indexPartOrig].position.getZ() << std::endl;
                    std::cout << "FMM  x " << positionsX[idxPart] << " y " << positionsY[idxPart] << " z " << positionsZ[idxPart] << std::endl;
                    std::cout << "Good fx " <<particles[indexPartOrig].forces[0] << " fy " << particles[indexPartOrig].forces[1] << " fz " << particles[indexPartOrig].forces[2] << std::endl;
                    std::cout << "FMM  fx " << forcesX[idxPart]*coeff_MD1 << " fy " << forcesY[idxPart]*coeff_MD1 << " fz " << forcesZ[idxPart]*coeff_MD1 << std::endl;
                    std::cout << "GOOD physical value " << particles[indexPartOrig].physicalValue << " potential " << particles[indexPartOrig].potential << std::endl;
                    std::cout << "FMM  physical value " << physicalValues[idxPart] << " potential " << potentials[idxPart] << std::endl;

                    std::cout << "\n";
                }
            }
        });
	
        printf("Difference between FMM and poly:\n");
        printf("Potential diff is = \n");
        printf("%e\n",potentialDiff.getL2Norm());
        printf("%e\n",potentialDiff.getInfNorm());
        printf("Fx diff is = \n");
        printf("%e\n",fx.getL2Norm());
        printf("%e\n",fx.getInfNorm());
        printf("Fy diff is = \n");
        printf("%e\n",fy.getL2Norm());
        printf("%e\n",fy.getInfNorm());
        printf("Fz diff is = \n");
        printf("%e\n",fz.getL2Norm());
        printf("%e\n",fz.getInfNorm());
        std::cout << std::endl<< std::endl<< "Potential= " << std::setprecision(8) << potential*coeff_MD/2 << std::endl;

    }

    // -----------------------------------------------------

    { // get sum forces&potential
        FReal  potential = 0;
        FReal fx = 0.0, fy = 0.0, fz = 0.0;

        tree.forEachLeaf([&](LeafClass* leaf){
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const FReal*const forcesX = leaf->getTargets()->getForcesX();
            const FReal*const forcesY = leaf->getTargets()->getForcesY();
            const FReal*const forcesZ = leaf->getTargets()->getForcesZ();
            const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();

            for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                potential += potentials[idxPart];
                fx += forcesX[idxPart];
                fy += forcesY[idxPart];
                fz += forcesZ[idxPart];
            }
        });

        std::cout << "Foces Sum  x = " << fx << " y = " << fy << " z = " << fz << std::endl;
        std::cout << "Potential = " << std::setprecision(8) << potential*coeff_MD/2 << std::endl;
        std::cout << "Energy from file = " << std::setprecision(8) << loader.getEnergy() << std::endl;
	//        std::cout << "Constante DL_POLY: " << coeff_MD << std::endl;
    }
    // -----------------------------------------------------

    // ReGenerate file
    if( FParameters::existParameter(argc, argv, "-gen") ){
        std::cout << "Generate ewal.out from input file" << std::endl;
        std::ofstream fileout("ewal.out",std::ifstream::out);
        std::ifstream file(filename,std::ifstream::in);
        if(file.is_open()){
            const int bufferSize = 512;
            char buffer[bufferSize];
            file.getline(buffer, bufferSize);
            fileout << buffer << '\n';
            file.getline(buffer, bufferSize);
            fileout << buffer << '\n';
            if( !FParameters::existParameter(argc, argv, "-noper") ){
                file.getline(buffer, bufferSize);
                fileout << buffer << '\n';
                file.getline(buffer, bufferSize);
                fileout << buffer << '\n';
                file.getline(buffer, bufferSize);
                fileout << buffer << '\n';
            }
            for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
                file.getline(buffer, bufferSize);
                fileout << buffer << '\n';
                file.getline(buffer, bufferSize);
                fileout << buffer << '\n';
                file.getline(buffer, bufferSize);
                file.getline(buffer, bufferSize);
            }
        }
    }
    // end generate

    // -----------------------------------------------------

    delete[] particles;
    delete[] particlesDirect;

    return 0;
}



