
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
#include "../../Src/Files/FRandomLoader.hpp"



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


int restultIndex(const int idxP, const int idxH, const int minP,
                 const int maxP, const int minH){
    return (idxH-minH)*(maxP+1-minP)+(idxP-minP);
}


typedef IndexedParticle                ParticleClass;
typedef FSphericalCell                 CellClass;
typedef FVector<ParticleClass>         ContainerClass;

typedef FSimpleLeaf<ParticleClass, ContainerClass >                     LeafClass;
typedef FOctree<ParticleClass, CellClass, ContainerClass , LeafClass >  OctreeClass;
typedef FSphericalKernel<ParticleClass, CellClass, ContainerClass >     KernelClass;

typedef FFmmAlgorithm<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;
typedef FFmmAlgorithmThread<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassThread;
typedef FFmmAlgorithmTask<OctreeClass, ParticleClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClassTask;


void doATest(const int NbParticles, const int minP, const int maxP, const int minH, const int maxH,
             const FReal physicalValue, const bool neutral ,
             FMath::FAccurater* allPotentialDiff, FReal* allAbsoluteDiff, double* timing,
             const int SizeSubLevels = 3, FReal* timeForDirect = 0){
    FTic counter;
    FRandomLoader<ParticleClass> loader(NbParticles);
    ParticleClass*const particles = new ParticleClass[loader.getNumberOfParticles()];

    const bool computeDirectAndDiff = timeForDirect || allAbsoluteDiff || allPotentialDiff;
    {
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            loader.fillParticle(particles[idxPart]);
            particles[idxPart].setIndex( idxPart );
            if((idxPart & 1) && neutral){
                particles[idxPart].setPhysicalValue( -physicalValue );
            }
            else particles[idxPart].setPhysicalValue( physicalValue );
        }

        // Compute direct
        if(computeDirectAndDiff){
            printf("Compute direct!\n");
            KernelClass kernels(minP, minH, loader.getBoxWidth(), loader.getCenterOfBox());
            counter.tic();
            for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
                for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                    kernels.directInteractionMutual(&particles[idxTarget], &particles[idxOther]);
                }
            }
            if(timeForDirect) *timeForDirect = counter.tacAndElapsed();
        }
    }

    for(int idxH = minH ; idxH <= maxH ; ++idxH){
        for(int idxP = minP ; idxP <= maxP ; ++idxP){
            std::cout << "Running!!! H = " << idxH << " P = " << idxP << std::endl;
            const int outputIndex = restultIndex(idxP, idxH, minP, maxP, minH);

            // -----------------------------------------------------
            CellClass::Init(idxP);
            OctreeClass tree(idxH, SizeSubLevels>=idxH?idxH-1:SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());

            // -----------------------------------------------------

            std::cout << "Creating & Inserting " << loader.getNumberOfParticles() << " particles ..." << std::endl;
            std::cout << "\tHeight : " << idxH << " \t sub-height : " << SizeSubLevels << std::endl;
            counter.tic();

            for(int idxPart = 0 ; idxPart < NbParticles ; ++idxPart){
                ParticleClass part = particles[idxPart];
                part.setPotential(0.0);
                part.setForces(0.0,0.0,0.0);
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

            FmmClass algo(&tree,&kernels);
            counter.tic();
            algo.execute();
            counter.tac();


            if(timing) timing[outputIndex] = counter.elapsed();
            std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl;

            // -----------------------------------------------------
            if(computeDirectAndDiff){
                FMath::FAccurater potentialDiff;
                FMath::FAccurater fx, fy, fz;
                FReal absoluteDiff = FReal(0.0);
                { // Check that each particle has been summed with all other
                    typename OctreeClass::Iterator octreeIterator(&tree);
                    octreeIterator.gotoBottomLeft();

                    do{
                        typename ContainerClass::BasicIterator leafIter(*octreeIterator.getCurrentListTargets());

                        while( leafIter.hasNotFinished() ){
                            const ParticleClass& other = particles[leafIter.data().getIndex()];

                            potentialDiff.add(other.getPotential(),leafIter.data().getPotential());
                            absoluteDiff = FMath::Max(absoluteDiff,FMath::Abs(other.getPotential()-leafIter.data().getPotential()));

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

                if(allPotentialDiff) allPotentialDiff[outputIndex] = potentialDiff;
                if(allAbsoluteDiff) allAbsoluteDiff[outputIndex] = absoluteDiff;
            }
        }
    }

    delete[] particles;
}

// Simply create particles and try the kernels
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable has to be used to test Spherical algorithm.\n";
    std::cout << ">> You can pass -help to know more\n";
    //////////////////////////////////////////////////////////////
    if( FParameters::existParameter(argc,argv,"-help") ){
        std::cout << ">> -test-hp [-p P] [-h H] [-sh SH] [-f file]\n";
        std::cout << "To test all P from 1 to P and all H from 2 to H\n";
        std::cout << "It creates the Potential errors and output a map.\n\n";

        std::cout << ">> -test-axe [-p P] [-sp start position] [-ep end position]\n";
        std::cout << "      [-h H] [-pv physical value]\n";
        std::cout << "To test between two particles from star position to end position\n";
        std::cout << "It check the accuracy.\n\n";

        std::cout << ">> -test-h [-h H] [-nb NB PART] [-pv physical value]\n";
        std::cout << "To test between for NB PART different H for 3 accuracies\n\n";

        std::cout << ">> -test-p [-p P] [-h H] [-nb NB PART] [-pv PHYSICAL VALUE]\n";
        std::cout << "To test the time for a run of different nb particlues\n\n";
    }

    const int NbParticles = FParameters::getValue(argc,argv,"-nb", 60000);
    const int NbLevels = FParameters::getValue(argc,argv,"-h", 6);
    const int DevP = FParameters::getValue(argc,argv,"-p", 30);
    const FReal physicalValue = FParameters::getValue(argc,argv,"-pv", 1.0);
    const bool neutral = FParameters::existParameter(argc,argv,"-neutral");
    const int SizeSubLevels = FParameters::getValue(argc,argv,"-sh", 3);

    if( FParameters::existParameter(argc,argv,"-test-hp") ){
        std::cout << "Execute : test-hp\n";

        FMath::FAccurater allPotentialDiff[(NbLevels+1-2)*(DevP+1-2)];
        printf("Size array %d \n", (NbLevels+1-2)*(DevP+1-2));

        doATest(NbParticles,2,DevP,2,NbLevels,
                physicalValue, neutral, allPotentialDiff, 0,
                0, SizeSubLevels);

        {
            FILE* fres = fopen("test-hp.res", "w");
            fprintf(fres, "# H P Res\n");

            for(int idxH = 2 ; idxH < NbLevels ; ++idxH){
                for(int idxP = 2 ; idxP < DevP ; ++idxP){
                    const int index = restultIndex(idxP, idxH, 2, DevP, 2);
                    fprintf(fres, "%d\t%d\t%e\t%e\n", idxH, idxP,
                            allPotentialDiff[index].getL2Norm(),
                            allPotentialDiff[index].getInfNorm());
                }
                fprintf(fres,"\n");
            }
            fclose(fres);
        }
        {
            FILE* fplot = fopen("test-hp.plot", "w");
            fprintf(fplot, "# H and P test\n");
            fprintf(fplot, "set terminal svg\n");
            fprintf(fplot, "set output \"test-hp.svg\"\n");
            fprintf(fplot, "set title \"Precision for different H and P (NbParticles = %d)\n", NbParticles);
            fprintf(fplot, "set xlabel \"H\"\n");
            fprintf(fplot, "set ylabel \"P\"\n");
            fprintf(fplot, "set zlabel \"Accuracy\"\n");
            fprintf(fplot, "set logscale z\n");
            fprintf(fplot, "unset surface\n");
            fprintf(fplot, "set pm3d\n");
            fprintf(fplot, "set view map\n");
            fprintf(fplot, "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"#ffffff\" behind\n");
            fprintf(fplot, "set macros\n");
            fprintf(fplot, "splot \\\n");
            fprintf(fplot, "\"./test-hp.res\" u 1:2:3 t \"Error L2\" ,\\\n");
            fprintf(fplot, "\"./test-hp.res\" u 1:3:4 t \"Error Inf\" ;\n");
        }
    }
    if( FParameters::existParameter(argc,argv,"-test-axe") ){
        std::cout << "Execute : test-axe\n";
        const FReal startPosition = FParameters::getValue(argc,argv,"-sp", 1.5);
        const FReal endPosition = FParameters::getValue(argc,argv,"-ep", 4.5);
        const FReal stepValue = FParameters::getValue(argc,argv,"-sep", 0.1);
        const int nbStep = int((endPosition-startPosition)/stepValue);

        const FReal boxWidth = FParameters::getValue(argc,argv,"-bw", 10.0);
        const FPoint boxCenter(0,0,0);


        FMath::FAccurater* potentialDiff = new FMath::FAccurater[DevP+1];

        for(int idxP = 1 ; idxP <= DevP ; ++idxP){
            for(int idxStep = 0 ; idxStep <= nbStep ; ++idxStep){
                ParticleClass centeredParticle;
                centeredParticle.setPosition(0.0,0.0,0.5);
                centeredParticle.setPhysicalValue(physicalValue);
                centeredParticle.setIndex(0);

                ParticleClass otherParticle;
                otherParticle.setPosition(0.0,0.0,startPosition + (idxStep*stepValue));
                otherParticle.setPhysicalValue(physicalValue);
                otherParticle.setIndex(1);

                CellClass::Init(idxP);
                OctreeClass tree(NbLevels, 2, boxWidth, boxCenter);

                tree.insert(centeredParticle);
                tree.insert(otherParticle);

                KernelClass kernels(idxP, NbLevels, boxWidth, boxCenter);
                FmmClass algo(&tree,&kernels);
                algo.execute();

                kernels.directInteractionMutual(&centeredParticle, &otherParticle);

                { // Check that each particle has been summed with all other
                    typename OctreeClass::Iterator octreeIterator(&tree);
                    octreeIterator.gotoBottomLeft();

                    do{
                        typename ContainerClass::BasicIterator leafIter(*octreeIterator.getCurrentListTargets());

                        while( leafIter.hasNotFinished() ){
                            const ParticleClass& other = (leafIter.data().getIndex()==0?centeredParticle:otherParticle);

                            potentialDiff[idxP].add(other.getPotential(),leafIter.data().getPotential());

                            leafIter.gotoNext();
                        }
                    } while(octreeIterator.moveRight());
                }

                // Print for information
                printf("For P %d, particle pos %lf\n", idxP, otherParticle.getPosition().getZ());
                printf("Potential diff is = %e \t %e\n",
                       potentialDiff[idxP].getL2Norm(),
                       potentialDiff[idxP].getInfNorm());
            }
        }


        {
            FILE* fres = fopen("test-axe.res", "w");
            fprintf(fres, "# P\tL2\tInf\n");

            for(int idxP = 1 ; idxP <= DevP ; ++idxP){
                fprintf(fres, "%d\t%e\t%e\n", idxP, potentialDiff[idxP].getL2Norm(),
                        potentialDiff[idxP].getInfNorm());
            }
            fclose(fres);
        }
        {
            FILE* fplot = fopen("test-axe.plot", "w");
            fprintf(fplot, "# Particles on Z test\n");
            fprintf(fplot, "set terminal svg\n");
            fprintf(fplot, "set output \"test-axe.svg\"\n");
            fprintf(fplot, "set title \"Precision for particles on the Z axis [%4.2f/%4.2f] H = %d, Box width %4.2f\"\n",
                    startPosition, endPosition, NbLevels, boxWidth);
            fprintf(fplot, "set xlabel \"P\"\n");
            fprintf(fplot, "set ylabel \"Precisions\"\n");
            fprintf(fplot, "set logscale y\n");
            fprintf(fplot, "set size ratio 0.5\n");
            fprintf(fplot, "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"#ffffff\" behind\n");
            fprintf(fplot, "set macros\n");
            fprintf(fplot, "mystyle = \"with lines\"\n");
            fprintf(fplot, "plot \\\n");
            fprintf(fplot, "\"./test-axe.res\" u 1:2 t \"Error L2\" @mystyle,\\\n");
            fprintf(fplot, "\"./test-axe.res\" u 1:3 t \"Error Inf\" @mystyle;\n");

        }

        delete[] potentialDiff;
    }
    if( FParameters::existParameter(argc,argv,"-test-h") ){
        std::cout << "Execute : test-h\n";
        const int DevsP[3] = {3, 8, 12};
        FReal timeCounter[3][NbLevels+1];
        FReal directTime;

        for(int idxP = 0 ; idxP < 3 ; ++idxP){
            doATest(NbParticles,DevsP[idxP],DevsP[idxP],2,NbLevels,
                    physicalValue, neutral, 0, 0,
                    timeCounter[idxP], SizeSubLevels, &directTime);
        }

        {
            FILE* fres = fopen("test-h.res", "w");
            fprintf(fres, "# H\tLow\tMoy\tHigh\n");
            fprintf(fres, "%d\t%e\t%e\t%e\n", 1, directTime, directTime, directTime);
            for(int idxH = 2 ; idxH <= NbLevels ; ++idxH){
                fprintf(fres, "%d\t%e\t%e\t%e\n", idxH, timeCounter[0][idxH-2],
                        timeCounter[1][idxH-2],timeCounter[2][idxH-2]);
            }
            fclose(fres);
        }
        {
            FILE* fplot = fopen("test-h.plot", "w");
            fprintf(fplot, "# several height test\n");
            fprintf(fplot, "set terminal svg\n");
            fprintf(fplot, "set output \"test-h.svg\"\n");
            fprintf(fplot, "set title \"Time for 3 different P and different H (NbParticles = %d)\"\n", NbParticles);
            fprintf(fplot, "set xlabel \"Height of the tree\"\n");
            fprintf(fplot, "set ylabel \"Time\"\n");
            fprintf(fplot, "set size ratio 0.5\n");
            fprintf(fplot, "set xtics autofreq 1\n");
            fprintf(fplot, "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"#ffffff\" behind\n");
            fprintf(fplot, "set macros\n");
            fprintf(fplot, "mystyle = \"with lines\"\n");
            fprintf(fplot, "plot \\\n");

            fprintf(fplot, "\"./test-h.res\" u 1:2 t \"Low Precision (P=%d)\" @mystyle,\\\n", DevsP[0]);
            fprintf(fplot, "\"./test-h.res\" u 1:3 t \"Medium Precision (P=%d)\" @mystyle,\\\n", DevsP[1]);
            fprintf(fplot, "\"./test-h.res\" u 1:4 t \"High Precision (P=%d)\" @mystyle;\n", DevsP[2]);
        }
    }
    if( FParameters::existParameter(argc,argv,"-test-time") ){
        std::cout << "Execute : test-time\n";
        const int DevsP[3] = {3, 8, 12};
        const int NbSteps = 4;
        const int ParticlesNumbers[NbSteps] = {100,1000,10000,100000};
        const int AllLevels[NbSteps] = {2,2,3,4};

        FReal timeCounter[3][NbSteps];
        FMath::FAccurater potentialDiff[3][NbSteps];
        FReal potentialAbsoluteDiff[3][NbSteps];

        for(int idxPart = 0 ; idxPart < NbSteps ; ++idxPart){
             for(int idxP = 0 ; idxP < 3 ; ++idxP){
                 doATest(ParticlesNumbers[idxPart],DevsP[idxP],DevsP[idxP],AllLevels[idxPart],AllLevels[idxPart],
                         physicalValue, neutral, potentialDiff[idxP]+idxPart, potentialAbsoluteDiff[idxP]+idxPart,
                         timeCounter[idxP]+idxPart, SizeSubLevels);
             }
        }

        {
            FILE* fres = fopen("test-np.res", "w");
            fprintf(fres, "# Particles");
            for(int idxP = 0 ; idxP < 4 ; ++idxP){
                fprintf(fres, "\t%d\t%d\t%d", DevsP[0],DevsP[1],DevsP[2]);
            }
            fprintf(fres, "\n");
            for(int idxPart = 0 ; idxPart < NbSteps ; ++idxPart){
                fprintf(fres, "%d", ParticlesNumbers[idxPart]);
                for(int idxP = 0 ; idxP < 3 ; ++idxP){
                    fprintf(fres, "\t%e", timeCounter[idxP][idxPart]);
                }
                for(int idxP = 0 ; idxP < 3 ; ++idxP){
                    fprintf(fres, "\t%e", potentialDiff[idxP][idxPart].getL2Norm());
                }
                for(int idxP = 0 ; idxP < 3 ; ++idxP){
                    fprintf(fres, "\t%e", potentialDiff[idxP][idxPart].getInfNorm());
                }
                for(int idxP = 0 ; idxP < 3 ; ++idxP){
                    fprintf(fres, "\t%e", potentialAbsoluteDiff[idxP][idxPart]);
                }
                fprintf(fres, "\n");
            }
            fclose(fres);
        }
        {
            FILE* fplot = fopen("test-np.plot", "w");
            fprintf(fplot, "# several height test\n");
            fprintf(fplot, "set terminal svg\n");
            fprintf(fplot, "set output \"test-np.svg\"\n");
            fprintf(fplot, "set title \"Time for three different P and different simulation\"\n");
            fprintf(fplot, "set xlabel \"Number of particles\"\n");
            fprintf(fplot, "set ylabel \"Time\"\n");
            fprintf(fplot, "set size ratio 0.5\n");
            fprintf(fplot, "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"#ffffff\" behind\n");
            fprintf(fplot, "set macros\n");
            fprintf(fplot, "mystyle = \"with lines\"\n");
            fprintf(fplot, "plot \\\n");

            fprintf(fplot, "\"./test-np.res\" u 1:2 t \"Time for low precision (P = %d)\" @mystyle,\\\n",DevsP[0]);
            fprintf(fplot, "\"./test-np.res\" u 1:3 t \"Time for medium precision (P = %d)\" @mystyle,\\\n",DevsP[1]);
            fprintf(fplot, "\"./test-np.res\" u 1:4 t \"Time for high precision (P = %d)\" @mystyle;\n",DevsP[2]);

            fprintf(fplot, "set terminal svg\n");
            fprintf(fplot, "set output \"test-np-error.svg\"\n");
            fprintf(fplot, "set title \"Error for three different P\"\n");
            fprintf(fplot, "set xlabel \"Number of particles\"\n");
            fprintf(fplot, "set ylabel \"Accuracy\"\n");
            fprintf(fplot, "set y2label \"Absolute Error\"\n");
            fprintf(fplot, "set size ratio 0.5\n");
            fprintf(fplot, "set logscale y\n");
            fprintf(fplot, "set logscale y2\n");
            fprintf(fplot, "set ytics nomirror\n");
            fprintf(fplot, "set y2tics nomirror\n");
            fprintf(fplot, "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"#ffffff\" behind\n");
            fprintf(fplot, "set macros\n");
            fprintf(fplot, "mystyle = \"with lines\"\n");
            fprintf(fplot, "plot \\\n");

            fprintf(fplot, "\"./test-np.res\" u 1:5 t \"L2 for low precision (P = %d)\" @mystyle axes x1y1,\\\n",DevsP[0]);
            fprintf(fplot, "\"./test-np.res\" u 1:6 t \"L2 for medium precision (P = %d)\" @mystyle axes x1y1,\\\n",DevsP[1]);
            fprintf(fplot, "\"./test-np.res\" u 1:7 t \"L2 for high precision (P = %d)\" @mystyle axes x1y1,\\\n",DevsP[2]);

            fprintf(fplot, "\"./test-np.res\" u 1:8 t \"Inf for low precision (P = %d)\" @mystyle axes x1y1,\\\n",DevsP[0]);
            fprintf(fplot, "\"./test-np.res\" u 1:9 t \"Inf for medium precision (P = %d)\" @mystyle axes x1y1,\\\n",DevsP[1]);
            fprintf(fplot, "\"./test-np.res\" u 1:10 t \"Inf for high precision (P = %d)\" @mystyle axes x1y1,\\\n",DevsP[2]);

            fprintf(fplot, "\"./test-np.res\" u 1:11 t \"Erreur for low precision (P = %d)\" @mystyle axes x1y2,\\\n",DevsP[0]);
            fprintf(fplot, "\"./test-np.res\" u 1:12 t \"Erreur for medium precision (P = %d)\" @mystyle axes x1y2,\\\n",DevsP[1]);
            fprintf(fplot, "\"./test-np.res\" u 1:13 t \"Erreur for high precision (P = %d)\" @mystyle axes x1y2;\n",DevsP[2]);

            fclose(fplot);
        }
    }
    if( FParameters::existParameter(argc,argv,"-test-p") ){
        std::cout << "Execute : test-p\n";
        FMath::FAccurater potentialDiff[DevP+1];
        FReal potentialAbsoluteDiff[DevP+1];

        doATest(NbParticles,1,DevP,NbLevels,NbLevels,
                physicalValue, neutral, potentialDiff, potentialAbsoluteDiff,
                0, SizeSubLevels);

        {
            FILE* fres = fopen("test-p.res", "w");
            fprintf(fres, "# P\tL2\tInf\n");
            for(int idxP = 1 ; idxP <= DevP ; ++idxP){
                fprintf(fres, "%d\t%e\t%e\t%e\n", idxP,
                        potentialDiff[idxP-1].getL2Norm(),
                        potentialDiff[idxP-1].getInfNorm(),
                        potentialAbsoluteDiff[idxP-1]);
            }
            fclose(fres);
        }
        {
            FILE* fplot = fopen("test-p.plot", "w");
            fprintf(fplot, "# several height test\n");
            fprintf(fplot, "set terminal svg\n");
            fprintf(fplot, "set output \"test-p.svg\"\n");
            fprintf(fplot, "set title \"Accuracy for different P (NbParticles = %d)\"\n", NbParticles);
            fprintf(fplot, "set xlabel \"P\"\n");
            fprintf(fplot, "set ylabel \"Accuracy\"\n");
            fprintf(fplot, "set y2label \"Absolute Error\"\n");
            fprintf(fplot, "set size ratio 0.5\n");
            fprintf(fplot, "set logscale y\n");
            fprintf(fplot, "set logscale y2\n");
            fprintf(fplot, "set ytics nomirror\n");
            fprintf(fplot, "set y2tics nomirror\n");
            fprintf(fplot, "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"#ffffff\" behind\n");
            fprintf(fplot, "set macros\n");
            fprintf(fplot, "mystyle = \"with lines\"\n");
            fprintf(fplot, "plot \\\n");

            fprintf(fplot, "\"./test-p.res\" u 1:2 t \"L2\" @mystyle axes x1y1,\\\n");
            fprintf(fplot, "\"./test-p.res\" u 1:3 t \"Inf\" @mystyle axes x1y1,\\\n");
            fprintf(fplot, "\"./test-p.res\" u 1:4 t \"Absolute\" @mystyle axes x1y2;\n");
        }
    }

    return 0;
}



