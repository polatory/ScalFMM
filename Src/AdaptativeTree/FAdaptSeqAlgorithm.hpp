// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, B��renger Bramas, Matthias Messner
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
#ifndef FADAPTSEQALGORITHM_HPP
#define FADAPTSEQALGORITHM_HPP

#include "../Utils/FGlobal.hpp"
#include "../Utils/FAssert.hpp"
#include "../Utils/FLog.hpp"
#include "../Utils/FTrace.hpp"
#include "../Utils/FTic.hpp"

#include "../Containers/FOctree.hpp"
#include "../Containers/FVector.hpp"

#include "FCoreCommon.hpp"

/**
 * @author Olivier Coulaud (olivier.coulaud@inria.fr)
 * @class FAdaptSeqAlgorithm
 * @brief
 * Please read the license
 *
 * This class is a basic FMM algorithm
 * It just iterates on a tree and call the kernels with good arguments.
 *
 * Of course this class does not deallocate pointer given in arguments.
 */
template<class OctreeClass, class CellClass, class ContainerClass, class KernelClass, class LeafClass>
class FAdaptSeqAlgorithm :  public FAbstractAlgorithm {

	OctreeClass* const tree;       //< The octree to work on
	KernelClass* const kernels;    //< The kernels

	const int OctreeHeight;

public:	
	/** The constructor need the octree and the kernels used for computation
	 * @param inTree the octree to work on
	 * @param inKernels the kernels to call
	 * An assert is launched if one of the arguments is null
	 */
	FAdaptSeqAlgorithm(OctreeClass* const inTree, KernelClass* const inKernels)
: tree(inTree) , kernels(inKernels), OctreeHeight(tree->getHeight()) {

		FAssertLF(tree, "tree cannot be null");
		FAssertLF(kernels, "kernels cannot be null");

		FLOG(FLog::Controller << "FFmmAlgorithm\n");
	}

	/** Default destructor */
	virtual ~FAdaptSeqAlgorithm(){
	}

	/**
	 * To execute the fmm algorithm
	 * Call this function to run the complete algorithm
	 */
	void execute(const unsigned operationsToProceed = FFmmNearAndFarFields){
		FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

		if(operationsToProceed & FFmmP2M) bottomPass();

		if(operationsToProceed & FFmmM2M) upwardPass();

		if(operationsToProceed & FFmmM2L) transferPass();

		if(operationsToProceed & FFmmL2L) downardPass();

		if( (operationsToProceed & FFmmP2P) || (operationsToProceed & FFmmL2P) ) directPass();
	}

private:
	/////////////////////////////////////////////////////////////////////////////
	// P2M
	/////////////////////////////////////////////////////////////////////////////

	/** P2M */
	void bottomPass(){
		FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
		FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
		FLOG(FTic counterTime);
		FLOG(FTic computationCounter);

		typename OctreeClass::Iterator octreeIterator(tree);

		// Iterate on leafs
		octreeIterator.gotoBottomLeft();
		do{
			// We need the current cell that represent the leaf
			// and the list of particles
			FLOG(computationCounter.tic());
			if( !octreeIterator.getCurrentCell()->iisSminMCriteria()){
				kernels->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
			}
			FLOG(computationCounter.tac());
		} while(octreeIterator.moveRight());

		FLOG( FLog::Controller << "\tFinished (@Bottom Pass (P2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
		FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
	}

	/////////////////////////////////////////////////////////////////////////////
	// Upward
	/////////////////////////////////////////////////////////////////////////////

	/** M2M */
	void upwardPass(){
		FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
		FLOG( FLog::Controller.write("\tStart Upward Pass\n").write(FLog::Flush); );
		FLOG(FTic counterTime);
		FLOG(FTic computationCounter);

		// Start from leal level - 1
		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.gotoBottomLeft();
		octreeIterator.moveUp();

		typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

		// for each levels
		for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
			FLOG(FTic counterTimeLevel);

			// for each cells
			do{
				// We need the current cell and the child
				// child is an array (of 8 child) that may be null
				FLOG(computationCounter.tic());
				//
				if(octreeIterator.getCurrentCell()->isAdaptative() &&  !octreeIterator.getCurrentCell()->isSminMCriteria() ){
					//		if( !octreeIterator.isSminMCriteria()){};
					// Need only adaptative child !!!!!
					const auto  * v =octreeIterator.getCurrentCell().getAdaptativeChild() ;
//					if (cell.sizeofAdaptativeChild()> 0 ){
//						for (int i=0; i < v->getSize() ; ++i){
//							output << v->operator [](i).cell->getGlobalId() << " ";
//						}
//					}

					//						kernels->M2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
				}
				FLOG(computationCounter.tac());
			} while(octreeIterator.moveRight());

			avoidGotoLeftIterator.moveUp();
			octreeIterator = avoidGotoLeftIterator;

			FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
		}


		FLOG( FLog::Controller << "\tFinished (@Upward Pass (M2M) = "  << counterTime.tacAndElapsed() << "s)\n" );
		FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
	}

	/////////////////////////////////////////////////////////////////////////////
	// Transfer
	/////////////////////////////////////////////////////////////////////////////

	/** M2L */
	void transferPass(){
		FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );

		FLOG( FLog::Controller.write("\tStart Downward Pass (M2L)\n").write(FLog::Flush); );
		FLOG(FTic counterTime);
		FLOG(FTic computationCounter);

		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.moveDown();

		typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

		const CellClass* neighbors[343];

		// for each levels
		for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
			FLOG(FTic counterTimeLevel);

			// for each cells
			do{
				const int counter = tree->getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);
				FLOG(computationCounter.tic());
				if(counter) kernels->M2L( octreeIterator.getCurrentCell() , neighbors, counter, idxLevel);
				FLOG(computationCounter.tac());
			} while(octreeIterator.moveRight());

			FLOG(computationCounter.tic());
			kernels->finishedLevelM2L(idxLevel);
			FLOG(computationCounter.tac());

			avoidGotoLeftIterator.moveDown();
			octreeIterator = avoidGotoLeftIterator;

			FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
		}
		FLOG( FLog::Controller << "\tFinished (@Downward Pass (M2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
		FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );
	}

	/////////////////////////////////////////////////////////////////////////////
	// Downward
	/////////////////////////////////////////////////////////////////////////////

	/** L2L */
	void downardPass(){
		FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
		FLOG( FLog::Controller.write("\tStart Downward Pass (L2L)\n").write(FLog::Flush); );
		FLOG(FTic counterTime);
		FLOG(FTic computationCounter );

		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.moveDown();

		typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

		const int heightMinusOne = OctreeHeight - 1;
		// for each levels exepted leaf level
		for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
			FLOG(FTic counterTimeLevel);

			// for each cells
			do{
				FLOG(computationCounter.tic());
				kernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
				FLOG(computationCounter.tac());
			} while(octreeIterator.moveRight());

			avoidGotoLeftIterator.moveDown();
			octreeIterator = avoidGotoLeftIterator;

			FLOG( FLog::Controller << "\t\t>> Level " << idxLevel << " = "  << counterTimeLevel.tacAndElapsed() << "s\n" );
		}

		FLOG( FLog::Controller << "\tFinished (@Downward Pass (L2L) = "  << counterTime.tacAndElapsed() << "s)\n" );
		FLOG( FLog::Controller << "\t\t Computation : " << computationCounter.cumulated() << " s\n" );


	}

	/////////////////////////////////////////////////////////////////////////////
	// Direct
	/////////////////////////////////////////////////////////////////////////////

	/** P2P */
	void directPass(){
		FTRACE( FTrace::FFunction functionTrace(__FUNCTION__, "Fmm" , __FILE__ , __LINE__) );
		FLOG( FLog::Controller.write("\tStart Direct Pass\n").write(FLog::Flush); );
		FLOG(FTic counterTime);
		FLOG(FTic computationCounterL2P);
		FLOG(FTic computationCounterP2P);

		const int heightMinusOne = OctreeHeight - 1;

		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.gotoBottomLeft();
		// There is a maximum of 26 neighbors
		ContainerClass* neighbors[27];
		// for each leafs
		do{
			FLOG(computationCounterL2P.tic());
			kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
			FLOG(computationCounterL2P.tac());
			// need the current particles and neighbors particles
			const int counter = tree->getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),heightMinusOne);
			FLOG(computationCounterP2P.tic());
			kernels->P2P(octreeIterator.getCurrentGlobalCoordinate(),octreeIterator.getCurrentListTargets(),
					octreeIterator.getCurrentListSrc(), neighbors, counter);
			FLOG(computationCounterP2P.tac());
		} while(octreeIterator.moveRight());


		FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
		FLOG( FLog::Controller << "\t\t Computation L2P : " << computationCounterL2P.cumulated() << " s\n" );
		FLOG( FLog::Controller << "\t\t Computation P2P : " << computationCounterP2P.cumulated() << " s\n" );

	}

};


#endif //FADAPTSEQALGORITHM_HPP

