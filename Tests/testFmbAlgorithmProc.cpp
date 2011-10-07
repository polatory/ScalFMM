// /!\ Please, you must read the license at the bottom of this page


#include "../Src/Utils/FTic.hpp"
#include "../Src/Utils/FMpi.hpp"
#include "../Src/Utils/FAbstractSendable.hpp"
#include "../Src/Utils/FParameters.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Components/FFmaParticle.hpp"
#include "../Src/Extensions/FExtendForces.hpp"
#include "../Src/Extensions/FExtendPotential.hpp"

#include "../Src/Components/FBasicCell.hpp"
#include "../Src/Fmb/FExtendFmbCell.hpp"

#include "../Src/Core/FFmmAlgorithmThreadProc.hpp"
#include "../Src/Core/FFmmAlgorithmThread.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Fmb/FFmbKernels.hpp"

#include "../Src/Files/FMpiFmaLoader.hpp"
#include "../Src/Files/FMpiTreeBuilder.hpp"
#include "../Src/Files/FFmaBinLoader.hpp"

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

//#define VALIDATE_FMM

// With openmp : mpicxx -g testFmbAlgorithmProc.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp ../Src/Utils/FMath.cpp ../Src/Utils/F3DPosition.cpp -lgomp -fopenmp -O2 -o testFmbAlgorithmProc.exe
// mpicxx -openmp testFmbAlgorithmProc.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp ../Src/Utils/FMath.cpp ../Src/Utils/F3DPosition.cpp -O2 -o testFmbAlgorithmProc.exe
// mpirun -perhost 1 -trace -l -np 8 ../testFmbAlgorithmProc.exe /lustre/bramas/200kk.fma.bin -h 9

// mpirun -np 3 `eztrace -e ./Tests/Debug/testFmbAlgorithmProc ../Data/testFMAlgorithm.fma.tmp `
// eztrace_convert -o my_paje /tmp/berenger_eztrace_log_rank_0 /tmp/berenger_eztrace_log_rank_1 /tmp/berenger_eztrace_log_rank_2

//mpicxx -I$VT_ROOT/include -trace -openmp testFmbAlgorithmProc.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp ../Src/Utils/FMath.cpp ../Src/Utils/F3DPosition.cpp -O2 -o testFmbAlgorithmProc.exe

// mpirun -np 8 -npernode 1 -output-filename tempmpi/out.mpi.temp `eztrace -e ./testFmbAlgorithmProc.exe /lustre/bramas/200kk.fma.bin -h 9`
// eztrace_convert -o my_paje /tmp/berenger_eztrace_log_rank_0 /tmp/berenger_eztrace_log_rank_1 /tmp/berenger_eztrace_log_rank_2 /tmp/berenger_eztrace_log_rank_3 /tmp/berenger_eztrace_log_rank_4 /tmp/berenger_eztrace_log_rank_5 /tmp/berenger_eztrace_log_rank_6 /tmp/berenger_eztrace_log_rank_7

/** This program show an example of use of
  * the fmm basic algo
  * it also check that eachh particles is little or longer
  * related that each other
  */


/** Fmb class has to extend {FExtendForces,FExtendPotential,FExtendPhysicalValue}
  * Because we use fma loader it needs {FFmaParticle}
  */
class FmbParticle : public FFmaParticle, public FExtendForces, public FExtendPotential {
public:
};

/** Custom cell
  *
  */
class FmbCell : public FBasicCell, public FExtendFmbCell , public FAbstractSendable {
public:
    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    static const int SerializedSizeUp = sizeof(FComplexe)*MultipoleSize + sizeof(FBasicCell);
    void serializeUp(void* const buffer) const {
        memcpy(buffer, (FBasicCell*)this, sizeof(FBasicCell));
        memcpy((char*)(buffer) + sizeof(FBasicCell), multipole_exp, sizeof(FComplexe)*MultipoleSize );
    }
    void deserializeUp(const void* const buffer){
        memcpy((FBasicCell*)this, buffer, sizeof(FBasicCell));
        memcpy(multipole_exp, (char*)(buffer) + sizeof(FBasicCell), sizeof(FComplexe)*MultipoleSize );
    }

    static const int SerializedSizeDown = sizeof(FComplexe)*MultipoleSize + sizeof(FBasicCell);
    void serializeDown(void* const buffer) const {
        memcpy(buffer, (FBasicCell*)this, sizeof(FBasicCell));
        memcpy((char*)(buffer) + sizeof(FBasicCell), local_exp, sizeof(FComplexe)*MultipoleSize );
    }
    void deserializeDown(const void* const buffer){
        memcpy((FBasicCell*)this, buffer, sizeof(FBasicCell));
        memcpy(local_exp, (char*)(buffer) + sizeof(FBasicCell), sizeof(FComplexe)*MultipoleSize );
    }

    ///////////////////////////////////////////////////////
    // to test equality between good and potentialy bad solution
    ///////////////////////////////////////////////////////
    /** To compare data */
    bool isEqualPole(const FmbCell& other, FReal*const cumul){
        //return memcmp(multipole_exp, other.multipole_exp, sizeof(FComplexe)*MultipoleSize) == 0 &&
        //        memcmp(local_exp, other.local_exp, sizeof(FComplexe)*MultipoleSize) == 0;
        *cumul = 0.0;
        for(int idx = 0; idx < MultipoleSize; ++idx){
            *cumul += FMath::Abs( multipole_exp[idx].getImag() - other.multipole_exp[idx].getImag() );
            *cumul += FMath::Abs( multipole_exp[idx].getReal() - other.multipole_exp[idx].getReal() );
        }

        return *cumul < 0.0001;//FMath::LookEqual(cumul,FReal(0.0));
    }

    /** To compare data */
    bool isEqualLocal(const FmbCell& other, FReal*const cumul){
        //return memcmp(multipole_exp, other.multipole_exp, sizeof(FComplexe)*MultipoleSize) == 0 &&
        //        memcmp(local_exp, other.local_exp, sizeof(FComplexe)*MultipoleSize) == 0;
        *cumul = 0.0;
        for(int idx = 0; idx < MultipoleSize; ++idx){
            *cumul += FMath::Abs( local_exp[idx].getImag() - other.local_exp[idx].getImag() );
            *cumul += FMath::Abs( local_exp[idx].getReal() - other.local_exp[idx].getReal() );
        }

        return *cumul < 0.0001;//FMath::LookEqual(cumul,FReal(0.0));
    }
};

#ifdef VALIDATE_FMM
template<class OctreeClass, class ContainerClass>
void ValidateFMMAlgoProc(OctreeClass* const badTree,
                         OctreeClass* const valideTree){
    std::cout << "Check Result\n";
    {
        const int OctreeHeight = valideTree->getHeight();
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level >= 1 ; --level){
            while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
                octreeIteratorValide.moveRight();
            }

            do {
                if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << "Error index are not equal!" << std::endl;
                }
                else{
                    FReal cumul;
                    if( !octreeIterator.getCurrentCell()->isEqualPole(*octreeIteratorValide.getCurrentCell(),&cumul) ){
                        std::cout << "Pole Data are different." << " Cumul " << cumul << std::endl;
                    }
                    if( !octreeIterator.getCurrentCell()->isEqualLocal(*octreeIteratorValide.getCurrentCell(),&cumul) ){
                        std::cout << "Local Data are different." << " Cumul " << cumul << std::endl;
                    }
                }

            } while(octreeIterator.moveRight() && octreeIteratorValide.moveRight());

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            octreeIteratorValide.moveUp();
            octreeIteratorValide.gotoLeft();
        }
    }
    {
        // Check that each particle has been summed with all other
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        while(octreeIteratorValide.getCurrentGlobalIndex() != octreeIterator.getCurrentGlobalIndex()){
            octreeIteratorValide.moveRight();
        }

        do {
            typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());

            typename ContainerClass::BasicIterator iterValide(*octreeIteratorValide.getCurrentListTargets());

            if( octreeIterator.getCurrentListSrc()->getSize() != octreeIteratorValide.getCurrentListSrc()->getSize()){
                std::cout << " Particules numbers is different " << std::endl;
            }
            if( octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                std::cout << " Index are differents " << std::endl;
            }

            while( iter.hasNotFinished() && iterValide.hasNotFinished() ){
                // If a particles has been impacted by less than NbPart - 1 (the current particle)
                // there is a problem

                if( !FMath::LookEqual(iterValide.data().getPotential() , iter.data().getPotential()) ){
                    std::cout << " Potential error : " << iterValide.data().getPotential()  << " " << iter.data().getPotential() << "\n";
                }
                if( !FMath::LookEqual(iterValide.data().getForces().getX(),iter.data().getForces().getX())
                        || !FMath::LookEqual(iterValide.data().getForces().getY(),iter.data().getForces().getY())
                        || !FMath::LookEqual(iterValide.data().getForces().getZ(),iter.data().getForces().getZ()) ){
                    /*std::cout << idx << " Forces error : " << (iterValide.data().getForces().getX() - iter.data().getForces().getX())
                              << " " << (iterValide.data().getForces().getY() - iter.data().getForces().getY())
                              << " " << (iterValide.data().getForces().getZ() - iter.data().getForces().getZ()) << "\n";*/
                    std::cout << " Forces error : x " << iterValide.data().getForces().getX() << " " << iter.data().getForces().getX()
                              << " y " << iterValide.data().getForces().getY()  << " " << iter.data().getForces().getY()
                              << " z " << iterValide.data().getForces().getZ()  << " " << iter.data().getForces().getZ() << "\n";
                }
                iter.gotoNext();
                iterValide.gotoNext();
            }

        } while(octreeIterator.moveRight() && octreeIteratorValide.moveRight());
    }

    std::cout << "Done\n";
}
#endif


// Simply create particles and try the kernels
int main(int argc, char ** argv){
    typedef FmbParticle             ParticleClass;
    typedef FmbCell                 CellClass;
    typedef FVector<ParticleClass>  ContainerClass;

    typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
    typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
    typedef FFmbKernels<ParticleClass, CellClass, ContainerClass >          KernelClass;

    typedef FFmmAlgorithmThreadProc<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
    typedef FFmmAlgorithmThread<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassNoProc;
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test fmb algorithm.\n";
    //////////////////////////////////////////////////////////////

    FMpi app( argc, argv);

    const int NbLevels = FParameters::getValue(argc,argv,"-h", 9);
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);
    FTic counter;
    const char* const defaultFilename = "testLoaderFMA.fma"; //../../Data/ "testLoaderFMA.fma" "testFMAlgorithm.fma" Sphere.fma
    const char* filename;

    if(argc == 1){
        std::cout << "You have to give a .fma file in argument.\n";
        std::cout << "The program will try a default file : " << defaultFilename << "\n";
        filename = defaultFilename;
    }
    else{
        filename = argv[1];
        std::cout << "Opening : " << filename << "\n";
    }

    FMpiFmaLoader<ParticleClass> loader(filename, app);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());

    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    if( app.processCount() != 1){
        //////////////////////////////////////////////////////////////////////////////////
        // We sort our particles
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Create intervals ..." << std::endl;
        counter.tic();

        FMpiTreeBuilder<ParticleClass> builder;

        builder.splitAndSortFile(loader, NbLevels);

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
        // Now we can build the real tree
        //////////////////////////////////////////////////////////////////////////////////
        std::cout << "Create real tree ..." << std::endl;
        counter.tic();

        builder.intervalsToTree(tree);

        counter.tac();
        std::cout << "Done  " << "(" << counter.elapsed() << "s)." << std::endl;

        //////////////////////////////////////////////////////////////////////////////////
    }
    else{
        FFmaBinLoader<ParticleClass> loaderSeq(filename);
        ParticleClass partToInsert;
        for(FSize idxPart = 0 ; idxPart < loaderSeq.getNumberOfParticles() ; ++idxPart){
            loaderSeq.fillParticle(partToInsert);
            tree.insert(partToInsert);
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels(NbLevels,loader.getBoxWidth());
    FmmClass algo(app,&tree,&kernels);
    algo.execute();

    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

    { // get sum forces&potential
        FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Sum Result" , __FILE__ , __LINE__) );

        FReal potential = 0;
        F3DPosition forces;

        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
        do{
            typename ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
            while( iter.hasNotFinished()){
                potential += iter.data().getPotential() * iter.data().getPhysicalValue();
                forces += iter.data().getForces();

                iter.gotoNext();
            }
        } while(octreeIterator.moveRight());

        std::cout << "My potential is " << potential << std::endl;

        potential = app.reduceSum(potential);
        forces.setX(app.reduceSum(forces.getX()));
        forces.setY(app.reduceSum(forces.getY()));
        forces.setZ(app.reduceSum(forces.getZ()));


        if(app.isMaster()){
            std::cout << "Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
            std::cout << "Potential Sum = " << potential << std::endl;
        }
    }

#ifdef VALIDATE_FMM
    {
        OctreeClass treeValide(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
        {
            FFmaBinLoader<ParticleClass> loaderSeq(filename);
            ParticleClass partToInsert;
            for(FSize idxPart = 0 ; idxPart < loaderSeq.getNumberOfParticles() ; ++idxPart){
                loaderSeq.fillParticle(partToInsert);
                treeValide.insert(partToInsert);
            }
        }

        std::cout << "Working on particles ..." << std::endl;
        counter.tic();
        FmmClassNoProc algoValide(&treeValide,&kernels);
        algoValide.execute();
        counter.tac();
        std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

        FReal potentialValide = 0;
        F3DPosition forcesValide;

        typename OctreeClass::Iterator octreeIteratorValide(&treeValide);
        octreeIteratorValide.gotoBottomLeft();

        do{
            typename ContainerClass::ConstBasicIterator iterValide(*octreeIteratorValide.getCurrentListTargets());
            while( iterValide.hasNotFinished()){
                potentialValide += iterValide.data().getPotential() * iterValide.data().getPhysicalValue();
                forcesValide += iterValide.data().getForces();

                iterValide.gotoNext();
            }

        } while(octreeIteratorValide.moveRight());

        std::cout << "Valide Foces Sum  x = " << forcesValide.getX() << " y = " << forcesValide.getY() << " z = " << forcesValide.getZ() << std::endl;
        std::cout << "Valide Potential = " << potentialValide << std::endl;

        ValidateFMMAlgoProc<OctreeClass,ContainerClass>(&tree,&treeValide);
    }
#endif


    // -----------------------------------------------------

    return 0;
}


// [--LICENSE--]
