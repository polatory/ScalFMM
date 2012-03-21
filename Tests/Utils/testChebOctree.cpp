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

#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "../../Src/Utils/FTic.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Utils/FAssertable.hpp"
#include "../../Src/Utils/FPoint.hpp"

#include "../../Src/Kernels/Chebyshev/FChebParticle.hpp"
#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Chebyshev/FChebLeaf.hpp"


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
	const FPoint BoxCenter(.5, .5, .5);
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

			const FPoint Origin(BoxCenter - BoxWidth / FReal(2.));
			const FPoint LeafCellCenter(Origin.getX() + (FReal(LeafCell->getCoordinate().getX()) + FReal(.5)) * BoxWidthLeaf,
																			 Origin.getY() + (FReal(LeafCell->getCoordinate().getY()) + FReal(.5)) * BoxWidthLeaf,
																			 Origin.getZ() + (FReal(LeafCell->getCoordinate().getZ()) + FReal(.5)) * BoxWidthLeaf);

			const ContainerClass *const Particles = octreeIterator.getCurrentListSrc();
			ContainerClass::ConstBasicIterator particleIterator(*Particles);
			while(particleIterator.hasNotFinished()) {
				const FPoint distance(LeafCellCenter-particleIterator.data().getPosition());
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