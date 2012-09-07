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
#include "../../Src/Kernels/Spherical/FSphericalParticle.hpp"

#include "../../Src/Files/FEwalLoader.hpp"
#include "../../Src/Components/FSimpleLeaf.hpp"

/** Ewal particle is used in the gadget program
  * here we try to make the same simulation
  */



// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FEwalParticle<FSphericalParticle>    ParticleClass;
    typedef FSphericalCell          CellClass;
    typedef FVector<ParticleClass>  ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >   KernelClass;

    typedef FFmmAlgorithmPeriodic<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassNoPer;

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
    FEwalLoader<ParticleClass> loader(filename);
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

    ParticleClass * const particles = new ParticleClass[loader.getNumberOfParticles()];

    for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
        loader.fillParticle(particles[idxPart]);
        // reset forces and insert in the tree
        particles[idxPart].setIndex(idxPart);
        ParticleClass part = particles[idxPart];
        part.setForces(0,0,0);
        part.setPotential(0);
        tree.insert(part);
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Create kernel & run simu ..." << std::endl;
    counter.tic();

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
    }

    counter.tac();

    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;


    // -----------------------------------------------------

    // Do direct
    if(FParameters::existParameter(argc, argv, "-direct")){
        printf("Compute direct:\n");
        FMath::FAccurater fx, fy, fz;
        KernelClass kernels( DevP, NbLevels, loader.getBoxWidth(), loader.getCenterOfBox());

        for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
            ParticleClass part = particles[idxTarget];
            part.setForces(0,0,0);
            part.setPotential(0);
            // compute with all other except itself
            for(int idxOther = 0; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                if( idxOther != idxTarget ){
                    kernels.directInteraction(&part, particles[idxOther]);
                }
            }
            fx.add(particles[idxTarget].getForces().getX(),-part.getForces().getX()*coeff_MD1);
            fy.add(particles[idxTarget].getForces().getY(),-part.getForces().getY()*coeff_MD1);
            fz.add(particles[idxTarget].getForces().getZ(),-part.getForces().getZ()*coeff_MD1);

            // print result
            if(FParameters::existParameter(argc, argv, "-verbose")){
                std::cout << ">> index " << part.getIndex() << " type " << part.getType() << std::endl;
                std::cout << "Good x " << particles[idxTarget].getPosition().getX() << " y " << particles[idxTarget].getPosition().getY() << " z " << particles[idxTarget].getPosition().getZ() << std::endl;
                std::cout << "DIRECT  x " << part.getPosition().getX() << " y " << part.getPosition().getY() << " z " << part.getPosition().getZ() << std::endl;
                std::cout << "Good fx " <<particles[idxTarget].getForces().getX() << " fy " << particles[idxTarget].getForces().getY() << " fz " << particles[idxTarget].getForces().getZ() << std::endl;
                std::cout << "DIRECT  fx " << part.getForces().getX()*coeff_MD1 << " fy " << part.getForces().getY()*coeff_MD1 << " fz " << part.getForces().getZ()*coeff_MD1 << std::endl;
                std::cout << "GOOD physical value " << particles[idxTarget].getPhysicalValue() << " potential " << particles[idxTarget].getPotential() << std::endl;
                std::cout << "DIRECT  physical value " << part.getPhysicalValue() << " potential " << part.getPotential() << std::endl;
                std::cout << "\n";
            }
        }

        printf("Difference between direct and poly:\n");
        printf("Fx diff is = \n");
        printf("%e\n",fx.getL2Norm());
        printf("%e\n",fx.getInfNorm());
        printf("Fy diff is = \n");
        printf("%e\n",fy.getL2Norm());
        printf("%e\n",fy.getInfNorm());
        printf("Fz diff is = \n");
        printf("%e\n",fz.getL2Norm());
        printf("%e\n",fz.getInfNorm());
    }

    // -----------------------------------------------------
    {
        FReal  potential = 0;
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;

        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){

                const ParticleClass& part = particles[iter.data().getIndex()];

                if(FParameters::existParameter(argc, argv, "-verbose")){
                    std::cout << ">> index " << iter.data().getIndex() << " type " << iter.data().getType() << std::endl;
                    std::cout << "Good x " << part.getPosition().getX() << " y " << part.getPosition().getY() << " z " << part.getPosition().getZ() << std::endl;
                    std::cout << "FMM  x " << iter.data().getPosition().getX() << " y " << iter.data().getPosition().getY() << " z " << iter.data().getPosition().getZ() << std::endl;
                    std::cout << "Good fx " <<part.getForces().getX() << " fy " << part.getForces().getY() << " fz " << part.getForces().getZ() << std::endl;
                    std::cout << "FMM  fx " << iter.data().getForces().getX()*coeff_MD1 << " fy " << iter.data().getForces().getY()*coeff_MD1 << " fz " << iter.data().getForces().getZ()*coeff_MD1 << std::endl;
                    std::cout << "GOOD physical value " << part.getPhysicalValue() << " potential " << part.getPotential() << std::endl;
                    std::cout << "FMM  physical value " << iter.data().getPhysicalValue() << " potential " << iter.data().getPotential() << std::endl;
                    std::cout << "\n";
                }

                potential += iter.data().getPotential() * iter.data().getPhysicalValue();


                potentialDiff.add(part.getPotential(),-iter.data().getPotential());
                fx.add(part.getForces().getX(),-iter.data().getForces().getX()*coeff_MD1);
                fy.add(part.getForces().getY(),-iter.data().getForces().getY()*coeff_MD1);
                fz.add(part.getForces().getZ(),-iter.data().getForces().getZ()*coeff_MD1);

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());
	
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
        FPoint forces;
        OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished() ){
                potential += iter.data().getPotential() * iter.data().getPhysicalValue();
                forces += iter.data().getForces();

                if(FParameters::existParameter(argc, argv, "-verbose")){
                    std::cout << " " << iter.data().getIndex()+1 << " \t "<< iter.data().getType() << "  \t " <<
                                 std::setprecision(5)<< iter.data().getPosition().getX() << "  \t" <<
                                 iter.data().getPosition().getY() << "  \t" <<
                                 iter.data().getPosition().getZ() << "   Forces: \t"<<
                                 std::setprecision(8) << iter.data().getForces().getX()*coeff_MD1 << "  \t " <<
                                 iter.data().getForces().getY()*coeff_MD1 << "  \t " <<
                                 iter.data().getForces().getZ()*coeff_MD1 << std::endl;
                }

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());

        std::cout << "Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
        std::cout << "Potential = " << std::setprecision(8) << potential*coeff_MD/2 << std::endl;
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

    return 0;
}



