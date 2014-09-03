// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Berenger Bramas, Matthias Messner
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
#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Kernels/Spherical/FSphericalCell.hpp"
#include "../Src/Kernels/Spherical/FSphericalKernel.hpp"
#include "../Src/Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "../Src/Kernels/Interpolation/FInterpMatrixKernel.hpp"

#include "../Src/Components/FSimpleLeaf.hpp"
#include "../Src/Core/FFmmAlgorithmPeriodic.hpp"
#include "../Src/Files/FRandomLoader.hpp"

#include "FUTester.hpp"

#include "../Src/Components/FTestCell.hpp"
#include "../Src/Components/FTestKernels.hpp"



/** The class to run the test */
class TestSphericalDirectPeriodic : public FUTester<TestSphericalDirectPeriodic> {
    /** Here we test only the P2P */
    void TestPeriodicFmm(){
        typedef FSphericalCell            CellClass;
        typedef FP2PParticleContainerIndexed<>  ContainerClass;

        typedef FSphericalKernel<CellClass, ContainerClass >   KernelClass;

        typedef FSimpleLeaf< ContainerClass >                     LeafClass;
        typedef FOctree< CellClass, ContainerClass , LeafClass >  OctreeClass;

        typedef FFmmAlgorithmPeriodic<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

        // Parameters
        const int NbLevels        = 4;
        const int SizeSubLevels = 2;
        const int PeriodicDeep  = 2;
        const int DevP              = 14;
        const int NbParticles     = 100;

        FSphericalCell::Init(DevP);

        FRandomLoader loader(NbParticles);
        OctreeClass tree(NbLevels, SizeSubLevels, loader.getBoxWidth(), loader.getCenterOfBox());
        struct TestParticle{
            FPoint position;
            FReal forces[3];
            FReal physicalValue;
            FReal potential;
        };
        FReal coeff = -1.0, value = 0.10, sum = 0.0;
       TestParticle* const particles = new TestParticle[loader.getNumberOfParticles()];
        for(int idxPart = 0 ; idxPart < loader.getNumberOfParticles() ; ++idxPart){
            FPoint position;
            loader.fillParticle(&position);
            value *= coeff ;
            sum += value ;
            // put in tree
            tree.insert(position, idxPart, value);
            // get copy
            particles[idxPart].position         = position;
            particles[idxPart].physicalValue = value;
            particles[idxPart].potential = 0.0;
            particles[idxPart].forces[0] = 0.0;
            particles[idxPart].forces[1] = 0.0;
            particles[idxPart].forces[2] = 0.0;
        }

        // Run FMM
        Print("Fmm...");
        FmmClass algo(&tree,PeriodicDeep);
        KernelClass kernels( DevP, algo.extendedTreeHeight(), algo.extendedBoxWidth(), algo.extendedBoxCenter());
        algo.setKernel(&kernels);
        algo.execute();

        // Run Direct
        Print("Run direct...");
        FTreeCoordinate min, max;
        algo.repetitionsIntervals(&min, &max);
        for(int idxTarget = 0 ; idxTarget < loader.getNumberOfParticles() ; ++idxTarget){
            for(int idxOther = idxTarget + 1 ; idxOther < loader.getNumberOfParticles() ; ++idxOther){
                FP2PR::MutualParticles(particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                      particles[idxTarget].position.getZ(),particles[idxTarget].physicalValue,
                                      &particles[idxTarget].forces[0],&particles[idxTarget].forces[1],
                                      &particles[idxTarget].forces[2],&particles[idxTarget].potential,
                                particles[idxOther].position.getX(), particles[idxOther].position.getY(),
                                particles[idxOther].position.getZ(),particles[idxOther].physicalValue,
                                &particles[idxOther].forces[0],&particles[idxOther].forces[1],
                                      &particles[idxOther].forces[2],&particles[idxOther].potential);
            }
            for(int idxX = min.getX() ; idxX <= max.getX() ; ++idxX){
                for(int idxY = min.getY() ; idxY <= max.getY() ; ++idxY){
                    for(int idxZ = min.getZ() ; idxZ <= max.getZ() ; ++idxZ){
                        if(idxX ==0 && idxY == 0 && idxZ == 0) continue;
                        // next lines for test

                        const FPoint offset(loader.getBoxWidth() * FReal(idxX),
                                            loader.getBoxWidth() * FReal(idxY),
                                            loader.getBoxWidth() * FReal(idxZ));

                        for(int idxSource = 0 ; idxSource < NbParticles ; ++idxSource){
                            TestParticle source = particles[idxSource];
                            source.position += offset;

                            FP2PR::NonMutualParticles(
                                        source.position.getX(), source.position.getY(),
                                        source.position.getZ(),source.physicalValue,
                                        particles[idxTarget].position.getX(), particles[idxTarget].position.getY(),
                                          particles[idxTarget].position.getZ(),particles[idxTarget].physicalValue,
                                          &particles[idxTarget].forces[0],&particles[idxTarget].forces[1],
                                        &particles[idxTarget].forces[2],&particles[idxTarget].potential);
                        }
                    }
                }
            }
        }

        // Compare
        Print("Compute Diff...");
        FMath::FAccurater potentialDiff;
        FMath::FAccurater fx, fy, fz;
        FReal energy= 0.0 , energyD = 0.0 ;
        { // Check that each particle has been summed with all other

            tree.forEachLeaf([&](LeafClass* leaf){
                const FReal*const potentials = leaf->getTargets()->getPotentials();
                const FReal*const physicalValues = leaf->getTargets()->getPhysicalValues();
                const FReal*const forcesX     = leaf->getTargets()->getForcesX();
                const FReal*const forcesY     = leaf->getTargets()->getForcesY();
                const FReal*const forcesZ     = leaf->getTargets()->getForcesZ();
                const int nbParticlesInLeaf    = leaf->getTargets()->getNbParticles();
                const FVector<int>& indexes = leaf->getTargets()->getIndexes();

                for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                    const int indexPartOrig = indexes[idxPart];
                    potentialDiff.add(particles[indexPartOrig].potential,potentials[idxPart]);
                    fx.add(particles[indexPartOrig].forces[0],forcesX[idxPart]);
                    fy.add(particles[indexPartOrig].forces[1],forcesY[idxPart]);
                    fz.add(particles[indexPartOrig].forces[2],forcesZ[idxPart]);
                    energy   += potentials[idxPart]*physicalValues[idxPart];
                    energyD +=particles[indexPartOrig].potential*particles[indexPartOrig].potential;

                }
            });
        }

        Print("Potential diff is = ");
        printf("         L2Norm   %e\n",potentialDiff.getRelativeL2Norm());
		printf("         RMSError %e\n",potentialDiff.getRMSError());
        Print("Fx diff is = ");
		printf("         L2Norm   %e\n",fx.getRelativeL2Norm());
		printf("         RMSError %e\n",fx.getRMSError());
        Print(fx.getRelativeL2Norm());
        Print(fx.getRelativeInfNorm());
        Print("Fy diff is = ");
		printf("        L2Norm   %e\n",fy.getRelativeL2Norm());
		printf("        RMSError %e\n",fy.getRMSError());
        Print("Fz diff is = ");
		printf("        L2Norm   %e\n",fz.getRelativeL2Norm());
		printf("        RMSError %e\n",fz.getRMSError());
        FReal L2error = (fx.getRelativeL2Norm()*fx.getRelativeL2Norm() + fy.getRelativeL2Norm()*fy.getRelativeL2Norm()  + fz.getRelativeL2Norm() *fz.getRelativeL2Norm()  );
		printf(" Total L2 Force Error= %e\n",FMath::Sqrt(L2error)) ;
		printf("  Energy Error =   %.12e\n",FMath::Abs(energy-energyD));

        const FReal MaximumDiff = FReal(0.001);
        uassert(potentialDiff.getRelativeL2Norm() < MaximumDiff);  // 1
        uassert(potentialDiff.getRMSError() < MaximumDiff);  // 2
        uassert(fx.getRelativeL2Norm()  < MaximumDiff);     // 3
        uassert(fx.getRMSError() < MaximumDiff);                // 4
        uassert(fy.getRelativeL2Norm()  < MaximumDiff);     // 5
        uassert(fy.getRMSError() < MaximumDiff);                 // 6
        uassert(fz.getRelativeL2Norm()  < MaximumDiff);    // 7
        uassert(fz.getRMSError() < MaximumDiff);                // 8

        delete[] particles;
    }

    /** After check the memory if needed */
    void After() {
        if( FMemStats::controler.isUsed() ){
            std::cout << "Memory used at the end " << FMemStats::controler.getCurrentAllocated() << " Bytes (" << FMemStats::controler.getCurrentAllocatedMB() << "MB)\n";
            std::cout << "Max memory used " << FMemStats::controler.getMaxAllocated() << " Bytes (" << FMemStats::controler.getMaxAllocatedMB() << "MB)\n";
            std::cout << "Total memory used " << FMemStats::controler.getTotalAllocated() << " Bytes (" << FMemStats::controler.getTotalAllocatedMB() << "MB)\n";
        }
    }

    // set test
    void SetTests(){        
        AddTest(&TestSphericalDirectPeriodic::TestPeriodicFmm,"Test Simu and with direct compare to Test fmm periodic");
    }
};


// You must do this
TestClass(TestSphericalDirectPeriodic)



