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

// ==== CMAKE =====
// @FUSE_BLAS
// @FUSE_STARPU
// ================

#include <iostream>

#include <cstdio>
#include <cstdlib>


#include "../../Src/Starpu/FCommon.hpp"

#include "../../Src/Files/FFmaScanfLoader.hpp"

#include "../../Src/Kernels/Chebyshev/FChebParticle.hpp"
#include "../../Src/Kernels/Chebyshev/FChebLeaf.hpp"
#include "../../Src/Kernels/Chebyshev/FChebCell.hpp"
#include "../../Src/Kernels/Chebyshev/FChebMatrixKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebSymKernel.hpp"
#include "../../Src/Kernels/Chebyshev/FChebFlopsSymKernel.hpp"

#include "../../Src/Utils/FParameters.hpp"

#include "../../Src/Containers/FOctree.hpp"
#include "../../Src/Containers/FVector.hpp"

#include "../../Src/Core/FFmmAlgorithm.hpp"
#include "../../Src/Core/FFmmAlgorithmThread.hpp"
#include "../../Src/Core/FFmmAlgorithmStarpuGroup.hpp"


#include "../Cuda/eventTimer.h"

class StarpuChebCell : public GChebCell {
public:
	void intialCopy(const GChebCell*const other){
		FMemUtils::copyall(getLocal(), other->getLocal(), GChebCell::DataSize);
		FMemUtils::copyall(getMultipole(), other->getMultipole(), GChebCell::DataSize);
	}
	void copyUp(const GChebCell*const other){
		FMemUtils::copyall(getMultipole(), other->getMultipole(), GChebCell::DataSize);
	}
	void restoreCopy(const GChebCell*const other){
		FMemUtils::copyall(getLocal(), other->getLocal(), GChebCell::DataSize);
		FMemUtils::copyall(getMultipole(), other->getMultipole(), GChebCell::DataSize);
	}
};



FReal computeL2norm(unsigned int N, FReal *const u, FReal *const v)
{
	FReal      dot = FReal(0.);
	FReal diff_dot = FReal(0.);
	for (unsigned int i=0; i<N; ++i) {
		FReal w = v[i] - u[i];
		diff_dot += w    * w;
		dot      += u[i] * u[i];
	}
	return FMath::Sqrt(diff_dot / dot);
}



FReal computeINFnorm(unsigned int N, FReal *const u, FReal *const v)
{
	FReal      max = FReal(0.);
	FReal diff_max = FReal(0.);
	for (unsigned int n=0; n<N; ++n) {
		if (     max<std::abs(u[n]))           max = std::abs(u[n]);
		if (diff_max<std::abs(u[n]-v[n])) diff_max = std::abs(u[n]-v[n]);
	}
	return diff_max / max;
}






// a function using usual starpu codelet prototype
void m2l_copy_const_data_cuda_func(void* constant_buffer)
{
	// copy permutation vectors to global memory
	int *const _pvectors = (int*)constant_buffer;
	cudaMemcpyToSymbol("pvectors", _pvectors, 343 * NNODES * sizeof(int));
	
	// copy permutation indices to constant memory
	int *const _pindices = (int*)constant_buffer + 343 * NNODES;
	cudaMemcpyToSymbol("pindices", _pindices, 343 * sizeof(int));
	
	// copy incremental lowranks to constant memory
	int *const _lowranks = (int*)constant_buffer + 343 * (NNODES + 1);
	cudaMemcpyToSymbol("lowranks", _lowranks, 16 * sizeof(int));
	
	// copy compressed m2l operators (K~UV') to global memory
	FReal *const _K = (FReal*)((int*)constant_buffer + (343 * (NNODES + 1) + 16));
	cudaMemcpyToSymbol("K", _K, _lowranks[15] * NNODES * 2 * sizeof(FReal));

	// copy root box width to constant memory
	FReal *const _width = _K + _lowranks[15] * NNODES * 2;
	cudaMemcpyToSymbol("width", _width, 1 * sizeof(FReal));

	// ask CUDA for the last error to occur (if one exists)
	cudaError_t error = cudaGetLastError();
        if(error != cudaSuccess) {
		// something's gone wrong
		// print out the CUDA error as a string
		printf("CUDA Error: %s\n", cudaGetErrorString(error));
	}
}










// Simply create particles and try the kernels
int main(int argc, char* argv[])
{
	const char* const filename       = FParameters::getStr(argc,argv,"-f", "../Data/test20k.fma");
	const unsigned int TreeHeight    = FParameters::getValue(argc, argv, "-h", 5);
        const unsigned int SubTreeHeight = FParameters::getValue(argc, argv, "-sh", 2);

	const unsigned int ORDER = CHEB_ORDER;
	const FReal epsilon = FReal(CHEB_EPSILON);

	std::cout << "CHEB_ORDER = " << ORDER << ", CHEB_EPSILON = " << epsilon << std::endl;

	//	// set threads
	//	omp_set_num_threads(NbThreads); 
	//	std::cout << "Using " << omp_get_max_threads() << " threads." << std::endl;

	//std::cout << "GCC:  sizeof(GChebCell) = " << sizeof(GChebCell) << std::endl;


	// init timer
	FTic time;

	// typedefs for STARPU
	typedef GParticle ParticleClass;
	typedef FVector<ParticleClass> ContainerClass;
	typedef FChebLeaf<ParticleClass,ContainerClass> LeafClass;
	typedef FChebMatrixKernelR MatrixKernelClass;
	typedef StarpuChebCell CellClass;
	typedef FOctree<ParticleClass,CellClass,ContainerClass,LeafClass> OctreeClass;
	//typedef FChebKernel<ParticleClass,CellClass,RealContainerClass,MatrixKernelClass,ORDER> KernelClass;
        typedef FChebFlopsSymKernel<ParticleClass,CellClass,ContainerClass,MatrixKernelClass,ORDER> FlopsKernelClass;
        typedef FChebSymKernel<ParticleClass,CellClass, ContainerClass,MatrixKernelClass,ORDER> KernelClass;

	// What we do //////////////////////////////////////////////////////
	std::cout << ">> Testing the Chebyshev interpolation base FMM algorithm.\n";
	
	// open particle file
	FFmaScanfLoader<ParticleClass> loader(filename);
	if(!loader.isOpen()) throw std::runtime_error("Particle file couldn't be opened!");
	
	// init oct-tree
	OctreeClass tree(TreeHeight, SubTreeHeight, loader.getBoxWidth(), loader.getCenterOfBox());

	// -----------------------------------------------------
	std::cout << "Creating and inserting " << loader.getNumberOfParticles() << " particles in a octree of height " << TreeHeight
						<< " ..." << std::endl;
	time.tic();
	loader.fillTree(tree);
	std::cout << "Done  " << "(" << time.tacAndElapsed() << ")." << std::endl;
	// -----------------------------------------------------


	
	// -----------------------------------------------------
	std::cout << "\nChebyshev FMM ... " << std::endl;
	time.tic();
	KernelClass kernels(TreeHeight, loader.getCenterOfBox(), loader.getBoxWidth(), epsilon);

	// -----------------------------------------------------
	// init constant memory on gpu
	char* constant_buffer = NULL;
	{
		const SymmetryHandler<ORDER> *const SymHandler = kernels.getPtrToSymHandler();
		int lowranks[16];
		int pindices[343];

		// global to local indirection
		unsigned int counter = 0;
		int g2l_idx[343];
		for (int idx=0; idx<343; ++idx)
			if (SymHandler->K[idx]!=NULL) g2l_idx[idx] = counter++;
			else                          g2l_idx[idx] = -1;

		// permutation vectors and permutated indices
		for (unsigned int idx=0; idx<343; ++idx) {
			const unsigned int gidx = SymHandler->pindices[idx];
			if (gidx != 0) {
				const unsigned int lidx = g2l_idx[gidx];
				pindices[idx] = lidx;
				lowranks[lidx] = SymHandler->LowRank[idx];
			} else pindices[idx] = -1;
		}

		// incremental low ranks
		for (unsigned int r=1; r<16; ++r)
			lowranks[r] = lowranks[r] + lowranks[r-1];

		//////////////////////////////////////////////
		// size of buffer to copu to gpu
		const unsigned int constant_size =
			343*NNODES            * sizeof(int) + // pvectors
			343                   * sizeof(int) + // pindices
			16                    * sizeof(int) + // lowranks
			2*lowranks[15]*NNODES * sizeof(FReal) + // K
			1                     * sizeof(FReal); // root box width
		//std::cout << "constant_size " << constant_size << std::endl;
		
		// memory //////////////////////////////////////
		constant_buffer = new char[constant_size];
		unsigned int msize;
		unsigned int offset = 0;
		// copy pvectors ///////////
		for (unsigned int idx=0; idx<343; ++idx) {
			msize = NNODES * sizeof(int);
			memcpy(constant_buffer + offset, SymHandler->pvectors[idx], msize);
			offset += msize;
		}
		// copy pindices ///////////
		msize = 343 * sizeof(int);
		memcpy(constant_buffer + offset, pindices, msize);
		offset += msize;
		// copy lowranks ///////////
		msize = 16 * sizeof(int);
		memcpy(constant_buffer + offset, lowranks, msize);
		offset += msize;
		// copy K //////////////////
		for (unsigned int idx=0; idx<343; ++idx)
			if (SymHandler->K[idx]!=NULL) {
				msize = 2*SymHandler->LowRank[idx]*NNODES * sizeof(FReal);
				memcpy(constant_buffer + offset, SymHandler->K[idx], msize);
				offset += msize;
			}
		// copy root box size //////
		msize = 1 * sizeof(FReal);
		const FReal width = loader.getBoxWidth();
		memcpy(constant_buffer + offset, &width, msize);
		offset += msize;
		//std::cout << "offset " << offset << std::endl;
	}

        std::cout << "init fmm & kernel " << time.tacAndElapsed() << "sec." << std::endl;

	// -----------------------------------------------------

        time.tic();
        m2l_copy_const_data_cuda_func(constant_buffer);
        std::cout << "Copy precompute to cuda " << time.tacAndElapsed() << "sec." << std::endl;

	// -----------------------------------------------------
	
        FlopsKernelClass flopskernels(TreeHeight, loader.getCenterOfBox(), loader.getBoxWidth(), epsilon);
        int nbCellsAtBottom = 0;
        {
            OctreeClass::Iterator iterator(&tree);
            iterator.gotoBottomLeft();
            do{
                nbCellsAtBottom++;
            }while(iterator.moveRight());
        }

        CellClass*const cells = new CellClass[nbCellsAtBottom];
        {
            const int idxLevel = TreeHeight-1;
            const CellClass* neighbors[343];
            OctreeClass::Iterator octreeIterator(&tree);
            octreeIterator.gotoBottomLeft();
            int counterCell = 0;
            do{
                //to have non zeros on multipole
                kernels.P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
                //copy
                cells[counterCell++] = *octreeIterator.getCurrentCell();
                //compute the flop
                const int counter = tree.getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);
                if(counter) flopskernels.M2L( octreeIterator.getCurrentCell() , neighbors, counter, idxLevel);

            }while(octreeIterator.moveRight());
        }
        const long long flops = flopskernels.flopsM2L;
        std::cout << "There are " << nbCellsAtBottom << " cells for a total of " << flops << " flops\n";

        // -----------------------------------------------------

        CellClass* gcells;

        time.tic();
        cudaMalloc(&gcells, sizeof(CellClass) * nbCellsAtBottom);
        cudaMemcpy(gcells, cells, sizeof(CellClass) * nbCellsAtBottom, cudaMemcpyHostToDevice);
        time.tac();
        std::cout <<"Allocate and copy to gpu " << time.elapsed() << std::endl;

        // -----------------------------------------------------
        time.tic();
        {
            const int idxLevel = TreeHeight-1;
            const CellClass* neighbors[343];
            for(int idx = 0 ; idx < nbCellsAtBottom ; ++idx){
                const int counter = tree.getInteractionNeighbors(neighbors, cells[idx].getCoordinate(), idxLevel);
                if(counter) kernels.M2L( &cells[idx] , neighbors, counter, idxLevel);
            }
        }
        time.tac();
        std::cout <<"CPU has taken " << time.elapsed() << "s => " << (FReal(flops)/time.elapsed()) << std::endl;
	// -----------------------------------------------------

        time.tic();
        {
            const int idxLevel = TreeHeight-1;
            eventTimerType t;
            initEventTimer(&t);
            startEventTimer(&t);
            perf_M2L_cuda_func( gcells,	nbCellsAtBottom, cells[0].getMortonIndex(), cells[nbCellsAtBottom-1].getMortonIndex(), idxLevel);
            stopEventTimer(&t);
            std::cout << "Timer " << getEventTimer(&t) << std::endl;
        }
        time.tac();
        std::cout <<"GPU has taken " << time.elapsed() << "s => " << (FReal(flops)/time.elapsed()) << std::endl;

        // -----------------------------------------------------

        CellClass*const resultcells = new CellClass[nbCellsAtBottom];
        cudaMemcpy(resultcells, gcells, sizeof(CellClass) * nbCellsAtBottom, cudaMemcpyDeviceToHost);

        FMath::FAccurater accurate;

        for(int idx = 0 ; idx < nbCellsAtBottom ; ++idx){
            for(int idxLocal = 0 ; idxLocal < CellClass::DataSize ; ++idxLocal){
                accurate.add( cells[idx].getLocal()[idxLocal], resultcells[idx].getLocal()[idxLocal] );
            }
        }

        std::cout << "Accuracy is inf = " << accurate.getInfNorm() << " l2norm = " << accurate.getL2Norm() << std::endl;

        // -----------------------------------------------------

        cudaFree(gcells);
        delete[] cells;
        delete[] resultcells;

	return 0;
}



