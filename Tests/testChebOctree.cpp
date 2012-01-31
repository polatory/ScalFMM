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
	
	const long NbPart = 200000;
	const FReal FRandMax = FReal(RAND_MAX);
	FTic counter;
	
	srand( static_cast<unsigned int>(time(NULL)) );
	
	OctreeClass tree(10, 3, 1.0, F3DPosition(0.5,0.5,0.5));

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

	return 0;
}


// [--END--]
