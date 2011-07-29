// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include <stdio.h>
#include <stdlib.h>

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
#include "../Src/Core/FFmmAlgorithm.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"

#include "../Src/Fmb/FFmbKernels.hpp"

#include "../Src/Files/FFmaLoader.hpp"



//#define VALIDATE_FMM

// With openmp : mpicxx -g testFmbAlgorithmProc.cpp ../Src/Utils/FDebug.cpp ../Src/Utils/FTrace.cpp -lgomp -fopenmp -O2 -o testFmbAlgorithmProc.exe
// icpc -openmp -openmp-lib=compat testFmbAlgorithm.cpp ../Src/Utils/FAssertable.cpp ../Src/Utils/FDebug.cpp -O2 -o testFmbAlgorithm.exe

// mpirun -np 3 `eztrace -e ./Tests/Debug/testFmbAlgorithmProc ../Data/testFMAlgorithm.fma.tmp `
// eztrace_convert -o my_paje /tmp/berenger_eztrace_log_rank_0 /tmp/berenger_eztrace_log_rank_1 /tmp/berenger_eztrace_log_rank_2

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
class FmbCell : public FBasicCell, public FExtendFmbCell , public FAbstractSendable{
public:
    ///////////////////////////////////////////////////////
    // to extend FAbstractSendable
    ///////////////////////////////////////////////////////
    int bytesToSendUp() const{
        return sizeof(FComplexe)*MultipoleSize;
    }
    int writeUp(void* const buffer, const int) const {
        memcpy(buffer,multipole_exp,bytesToSendUp());
        return bytesToSendUp();
    }
    int bytesToReceiveUp() const{
        return sizeof(FComplexe)*MultipoleSize;
    }
    int readUp(void* const buffer, const int) {
        memcpy(multipole_exp,buffer,bytesToSendUp());
        return bytesToReceiveUp();
    }

    int bytesToSendDown() const{
        return sizeof(FComplexe)*MultipoleSize;
    }
    int writeDown(void* const buffer, const int) const {
        memcpy(buffer,local_exp,bytesToSendDown());
        return bytesToSendDown();
    }
    int bytesToReceiveDown() const{
        return sizeof(FComplexe)*MultipoleSize;
    }
    int readDown(void* const buffer, const int) {
        FComplexe*const otherLocal = static_cast<FComplexe*>(buffer);
        for(int idx = 0 ; idx < MultipoleSize ; ++idx){
            local_exp[idx] += otherLocal[idx];
        }
        return bytesToReceiveDown();
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

        return *cumul == 0.0;//FMath::LookEqual(cumul,FReal(0.0));
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
template<class OctreeClass, class ContainerClass, class FmmClass>
void ValidateFMMAlgoProc(OctreeClass* const badTree,
                         OctreeClass* const valideTree,
                         FmmClass*const fmm){
    const int OctreeHeight = valideTree->getHeight();
    std::cout << "Check Result\n";
    {
        typename OctreeClass::Iterator octreeIterator(badTree);
        octreeIterator.gotoBottomLeft();

        typename OctreeClass::Iterator octreeIteratorValide(valideTree);
        octreeIteratorValide.gotoBottomLeft();

        for(int level = OctreeHeight - 1 ; level >= 1 ; --level){
            int NbLeafs = 0;
            do{
                ++NbLeafs;
            } while(octreeIterator.moveRight());
            octreeIterator.gotoLeft();

            const int startIdx = fmm->getLeft(NbLeafs);
            const int endIdx = fmm->getRight(NbLeafs);
            // Check that each particle has been summed with all other

            for(int idx = 0 ; idx < startIdx ; ++idx){
                octreeIterator.moveRight();
                octreeIteratorValide.moveRight();
            }

            for(int idx = startIdx ; idx < endIdx ; ++idx){
                if(octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << "Error index are not equal!" << std::endl;
                }
                else{
                    FReal cumul;
                    if( !octreeIterator.getCurrentCell()->isEqualPole(*octreeIteratorValide.getCurrentCell(),&cumul) ){
                        std::cout << "Pole Data are different." << idx << " Cumul " << cumul << std::endl;
                    }
                    if( !octreeIterator.getCurrentCell()->isEqualLocal(*octreeIteratorValide.getCurrentCell(),&cumul) ){
                        std::cout << "Local Data are different." << idx << " Cumul " << cumul << std::endl;
                    }
                }

                octreeIterator.moveRight();
                octreeIteratorValide.moveRight();
            }

            octreeIterator.moveUp();
            octreeIterator.gotoLeft();

            octreeIteratorValide.moveUp();
            octreeIteratorValide.gotoLeft();
        }
    }
    {
        int NbLeafs = 0;
        { // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();
            do{
                ++NbLeafs;
            } while(octreeIterator.moveRight());
            std::cout << "There is " << NbLeafs << " Leafs" << std::endl;
        }
        {
            const int startIdx = fmm->getLeft(NbLeafs);
            const int endIdx = fmm->getRight(NbLeafs);
            // Check that each particle has been summed with all other
            typename OctreeClass::Iterator octreeIterator(badTree);
            octreeIterator.gotoBottomLeft();

            typename OctreeClass::Iterator octreeIteratorValide(valideTree);
            octreeIteratorValide.gotoBottomLeft();

            for(int idx = 0 ; idx < startIdx ; ++idx){
                octreeIterator.moveRight();
                octreeIteratorValide.moveRight();
            }

            for(int idx = startIdx ; idx < endIdx ; ++idx){
                typename ContainerClass::BasicIterator iter(*octreeIterator.getCurrentListTargets());

                typename ContainerClass::BasicIterator iterValide(*octreeIteratorValide.getCurrentListTargets());

                if( octreeIterator.getCurrentListSrc()->getSize() != octreeIteratorValide.getCurrentListSrc()->getSize()){
                    std::cout << idx << " Particules numbers is different " << std::endl;
                }
                if( octreeIterator.getCurrentGlobalIndex() != octreeIteratorValide.getCurrentGlobalIndex()){
                    std::cout << idx << " Index are differents " << std::endl;
                }

                while( iter.hasNotFinished() && iterValide.hasNotFinished() ){
                    // If a particles has been impacted by less than NbPart - 1 (the current particle)
                    // there is a problem

                    if( iterValide.data().getPotential() != iter.data().getPotential() ){
                        std::cout << idx << " Potential error : " << iterValide.data().getPotential()  << " " << iter.data().getPotential() << "\n";
                    }
                    if( !FMath::LookEqual(iterValide.data().getForces().getX(),iter.data().getForces().getX())
                            || !FMath::LookEqual(iterValide.data().getForces().getY(),iter.data().getForces().getY())
                            || !FMath::LookEqual(iterValide.data().getForces().getZ(),iter.data().getForces().getZ()) ){
                        /*std::cout << idx << " Forces error : " << (iterValide.data().getForces().getX() - iter.data().getForces().getX())
                                  << " " << (iterValide.data().getForces().getY() - iter.data().getForces().getY())
                                  << " " << (iterValide.data().getForces().getZ() - iter.data().getForces().getZ()) << "\n";*/
                        std::cout << idx << " Forces error : x " << iterValide.data().getForces().getX() << " " << iter.data().getForces().getX()
                                  << " y " << iterValide.data().getForces().getY()  << " " << iter.data().getForces().getY()
                                  << " z " << iterValide.data().getForces().getZ()  << " " << iter.data().getForces().getZ() << "\n";
                    }
                    iter.gotoNext();
                    iterValide.gotoNext();
                }

                octreeIterator.moveRight();
                octreeIteratorValide.moveRight();
            }
        }
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
    typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassNoProc;
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

    FFmaLoader<ParticleClass> loader(filename);
    if(!loader.isOpen()){
        std::cout << "Loader Error, " << filename << " is missing\n";
        return 1;
    }

    // -----------------------------------------------------

    OctreeClass tree(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
#ifdef VALIDATE_FMM
    OctreeClass treeValide(NbLevels, SizeSubLevels,loader.getBoxWidth(),loader.getCenterOfBox());
#endif
    // -----------------------------------------------------

    std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    counter.tic();

    {
        ParticleClass particleToFill;
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(particleToFill);
            tree.insert(particleToFill);
            #ifdef VALIDATE_FMM
            treeValide.insert(particleToFill);
            #endif
        }
    }

    counter.tac();
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl;

    // -----------------------------------------------------

    std::cout << "Working on particles ..." << std::endl;
    counter.tic();

    KernelClass kernels(NbLevels,loader.getBoxWidth());

    //FmmClass algo(app,&tree,&kernels);
    //algo.execute();
#ifdef VALIDATE_FMM
    FmmClassNoProc algoValide(&treeValide,&kernels);
    algoValide.execute();
#endif
    counter.tac();
    std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;


    { // get sum forces&potential
        FReal potential = 0;
        F3DPosition forces;
#ifdef VALIDATE_FMM
        FReal potentialValide = 0;
        F3DPosition forcesValide;
#endif
        typename OctreeClass::Iterator octreeIterator(&tree);
        octreeIterator.gotoBottomLeft();
#ifdef VALIDATE_FMM
        typename OctreeClass::Iterator octreeIteratorValide(&treeValide);
        octreeIteratorValide.gotoBottomLeft();
#endif
        typename OctreeClass::Iterator countLeafsIterator(octreeIterator);
        int NbLeafs = 0;
        do{
            ++NbLeafs;
        } while(countLeafsIterator.moveRight());

        const int startIdx = 0;//algo.getLeft(NbLeafs);
        const int endIdx = 0;//algo.getRight(NbLeafs);

        std::cout <<"From " << startIdx << " to " << endIdx << "  NbLeafs is " << NbLeafs << std::endl;

        for(int idxLeaf = 0 ; idxLeaf < startIdx ; ++idxLeaf){
            octreeIterator.moveRight();
#ifdef VALIDATE_FMM
            octreeIteratorValide.moveRight();
#endif
        }

        for(int idxLeaf = startIdx ; idxLeaf < endIdx ; ++idxLeaf){
            typename ContainerClass::ConstBasicIterator iter(*octreeIterator.getCurrentListTargets());
#ifdef VALIDATE_FMM
            typename ContainerClass::ConstBasicIterator iterValide(*octreeIteratorValide.getCurrentListTargets());
#endif
            while( iter.hasNotFinished()
#ifdef VALIDATE_FMM
                  && iterValide.hasNotFinished()
#endif
                  ){
                potential += iter.data().getPotential() * iter.data().getPhysicalValue();
                forces += iter.data().getForces();
#ifdef VALIDATE_FMM
                potentialValide += iterValide.data().getPotential() * iterValide.data().getPhysicalValue();
                forcesValide += iterValide.data().getForces();
                iterValide.gotoNext();
#endif
                iter.gotoNext();
            }

            octreeIterator.moveRight();
#ifdef VALIDATE_FMM
            octreeIteratorValide.moveRight();
#endif
        }


#ifdef VALIDATE_FMM
        std::cout << "MPI Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
        std::cout << "Valide Foces Sum  x = " << forcesValide.getX() << " y = " << forcesValide.getY() << " z = " << forcesValide.getZ() << std::endl;
        std::cout << "MPI Potential = " << potential << std::endl;
        std::cout << "Valide Potential = " << potentialValide << std::endl;
#endif
        potential = app.reduceSum(potential);
        forces.setX(app.reduceSum(forces.getX()));
        forces.setY(app.reduceSum(forces.getY()));
        forces.setZ(app.reduceSum(forces.getZ()));
        if(app.isMaster()){
            std::cout << "Foces Sum  x = " << forces.getX() << " y = " << forces.getY() << " z = " << forces.getZ() << std::endl;
            std::cout << "Potential = " << potential << std::endl;
        }
    }

#ifdef VALIDATE_FMM
    ValidateFMMAlgoProc<OctreeClass,ContainerClass,FmmClass>(&tree,&treeValide,&algo);
#endif
    // -----------------------------------------------------

    return 0;
}


// [--LICENSE--]
