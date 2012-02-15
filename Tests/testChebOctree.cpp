// [--License--]

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../Src/Utils/FTic.hpp"

#include "../Src/Containers/FOctree.hpp"
#include "../Src/Containers/FVector.hpp"

#include "../Src/Utils/FAssertable.hpp"
#include "../Src/Utils/F3DPosition.hpp"

#include "../Src/Chebyshev/FChebParticle.hpp"
#include "../Src/Chebyshev/FChebCell.hpp"
#include "../Src/Chebyshev/FChebLeaf.hpp"


/**
* In this file we show how to use octree
*/

int main(int, char **){    

	const int ORDER = 5;

	typedef FChebParticle ParticleClass;
	typedef FVector<FChebParticle> ContainerClass;
	typedef FChebLeaf<ParticleClass,ContainerClass> LeafClass;
	typedef FChebCell<ORDER> CellClass;
	typedef FOctree<ParticleClass,CellClass,ContainerClass,LeafClass> OctreeClass;
	
	///////////////////////What we do/////////////////////////////
	std::cout << ">> This executable is useless to execute.\n";
	std::cout << ">> It is only interesting to understand the code\n";
	std::cout << ">> and how to use the Octree\n";
	//////////////////////////////////////////////////////////////
	
	const long NbPart = 100000;
	const FReal FRandMax = FReal(RAND_MAX);
	FTic counter;
	
	srand( static_cast<unsigned int>(time(NULL)) );
	
	const FReal BoxWidth = 1.;
	const F3DPosition BoxCenter(.5, .5, .5);
	const unsigned int TreeHeight = 10;
	OctreeClass tree(TreeHeight, 3, BoxWidth, BoxCenter);

	// -----------------------------------------------------
	std::cout << "Creating and inserting " << NbPart << " particles ..." << std::endl;
	counter.tic();
	{
		FChebParticle particle;
		for(long idxPart = 0 ; idxPart < NbPart ; ++idxPart){
			particle.setPosition(FReal(rand())/FRandMax,FReal(rand())/FRandMax,FReal(rand())/FRandMax);
			tree.insert(particle);
		}
	}
	counter.tac();
	std::cout << "Done  " << "(" << counter.elapsed() << ")." << std::endl;
	// -----------------------------------------------------


	
	// Check if particles are strictly within its containing leaf cells
	{
		const FReal BoxWidthLeaf = BoxWidth / FReal(FMath::pow(2, TreeHeight-1));
		OctreeClass::Iterator octreeIterator(&tree);
		octreeIterator.gotoBottomLeft();
		do{
			const CellClass *const LeafCell = octreeIterator.getCurrentCell();

			const F3DPosition Origin(BoxCenter - BoxWidth / FReal(2.));
			const F3DPosition LeafCellCenter(Origin.getX() + (FReal(LeafCell->getCoordinate().getX()) + FReal(.5)) * BoxWidthLeaf,
																			 Origin.getY() + (FReal(LeafCell->getCoordinate().getY()) + FReal(.5)) * BoxWidthLeaf,
																			 Origin.getZ() + (FReal(LeafCell->getCoordinate().getZ()) + FReal(.5)) * BoxWidthLeaf);

			const ContainerClass *const Particles = octreeIterator.getCurrentListSrc();
			ContainerClass::ConstBasicIterator particleIterator(*Particles);
			while(particleIterator.hasNotFinished()) {
				const F3DPosition distance(LeafCellCenter-particleIterator.data().getPosition());
				if (std::abs(distance.getX())>BoxWidthLeaf/FReal(2.) ||
						std::abs(distance.getY())>BoxWidthLeaf/FReal(2.) ||
						std::abs(distance.getZ())>BoxWidthLeaf/FReal(2.)) {
					std::cout << "Particle (center - particle = " << distance << " < " << BoxWidthLeaf/FReal(2.) << ") is out of cell. STOP"
										<< std::endl;
					//exit(-1);
				}
				particleIterator.gotoNext();
			}
		} while(octreeIterator.moveRight());
	}

	return 0;
}


// [--END--]
