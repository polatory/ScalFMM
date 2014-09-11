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
#ifndef FADAPTSEQALGORITHM_HPP
#define FADAPTSEQALGORITHM_HPP

#include <array>
#include "Utils/FGlobal.hpp"
#include "Utils/FAssert.hpp"
#include "Utils/FLog.hpp"

#include "Utils/FTic.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Core/FCoreCommon.hpp"

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

	OctreeClass* const tree;        //< The octree to work on
	KernelClass* const kernels;   //< The kernels

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
		FLOG( FLog::Controller.write("\tStart Bottom Pass\n").write(FLog::Flush) );
		FLOG(FTic counterTime);
		FLOG(FTic computationCounter);
		std::cout << "-->bottomPass"<<std::endl;

		typename OctreeClass::Iterator octreeIterator(tree);

		// Iterate on leafs
		octreeIterator.gotoBottomLeft();
		do{
			// We need the current cell that represent the leaf
			// and the list of particles
			FLOG(computationCounter.tic());
			if( !octreeIterator.getCurrentCell()->isSminMCriteria()){
				auto   * currentCell  = octreeIterator.getCurrentCell() ;

				//
				std::cout << "  P2M operator on cell " <<currentCell->getGlobalId() <<std::endl;
				// Classical P2M operator
				kernels->P2M( currentCell, octreeIterator.getCurrentListSrc());
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
		FLOG( FLog::Controller.write("\tStart Upward Pass\n").write(FLog::Flush); );
		FLOG(FTic counterTime);
		FLOG(FTic computationCounter);
		std::cout << "-->upwardPass"<<std::endl;
		// Start from leal level - 1
		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.gotoBottomLeft();
		octreeIterator.moveUp();

		typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

		// for each levels
		for(int idxLevel = OctreeHeight - 2 ; idxLevel > 1 ; --idxLevel ){
			FLOG(FTic counterTimeLevel);
			std::cout << "  Level "<<idxLevel<<std::endl;
			//
			std::array<int,8> m2m_child_level{ {0,0,0,0,0,0,0,0}}  ;
			//			std::array<KernelClass*,8> m2m_child_cell{ {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr}},
			CellClass*          m2m_child_cell[8] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};

			// for each cells
			do{
				// We need the current cell and the child
				// child is an array (of 8 child) that may be null
				FLOG(computationCounter.tic());
				//
				if(octreeIterator.getCurrentCell()->isadaptive() &&  !octreeIterator.getCurrentCell()->isSminMCriteria() ){
					m2m_child_level.fill(-1) ;
					// adaptive cell and enough work
					// Need only adaptive child !!!!!
					auto   * currentCell  = octreeIterator.getCurrentCell() ;
//				auto     currentKCell = currentCell->getKernelCell();
					int nbM2M =0;
					//
					for(int i=0 ; i<8 ; ++i){
						m2m_child_cell[i] =nullptr ;
					}
					if (currentCell->sizeofadaptiveChild()> 0 ){
						const auto  & cellChild      = currentCell->getAdaptiveChild() ;
						for (int i=0; i < cellChild.getSize() ; ++i){
							if (! cellChild[i].cell->isSminMCriteria() ){
								m2m_child_level[nbM2M] = cellChild[i].getLevel();
								std::cout << "     M2M operator to do between " <<currentCell->getGlobalId() <<" and "
										<< cellChild[i].cell->getGlobalId() << std::endl;
								nbM2M++ ;
							}
							else {
								std::cout << "     P2M operator to do between " <<currentCell->getGlobalId() <<" and "
										<< cellChild[i].cell->getGlobalId() << std::endl;
								kernels->P2MAdapt( currentCell, idxLevel );

							} //end if

						} // end loop
						kernels->M2MAdapt(currentCell, idxLevel, nbM2M,m2m_child_level.data(),m2m_child_cell);
					} // end if cell adaptive
					// Build Multipolaire vectors
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
		FLOG( FLog::Controller.write("\tStart Downward Pass (M2L)\n").write(FLog::Flush); );
		FLOG(FTic counterTime);
		FLOG(FTic computationCounter);
		std::cout << "-->transferPass"<<std::endl;

		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.moveDown();

		typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

		const CellClass* neighbors[343];

		// for each levels
		for(int idxLevel = 2 ; idxLevel < OctreeHeight ; ++idxLevel ){
			FLOG(FTic counterTimeLevel);
			std::cout << "  Level "<<idxLevel<<std::endl;

			// for each cells
			do{
				int counter = tree->getInteractionNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(), idxLevel);
				FLOG(computationCounter.tic());
				//				if(counter) kernels->M2L( octreeIterator.getCurrentCell() , neighbors, counter, idxLevel);
				if(counter) { //

					int nb = 0 ;

					auto  localExpCell = octreeIterator.getCurrentCell() ;
					if (localExpCell->sizeofadaptiveChild()  == 1){ // We only have 1
						localExpCell = (localExpCell->getAdaptiveChild(0))->cell ;
						//level ??
					}
					// The adaptive interaction list
					for (int i = 0 ; i < 343 ; ++i){
						if(neighbors[i]){
							if( ! neighbors[i]->isadaptive() ) { // Here we need Adaptive child
								nb+= neighbors[i]->sizeofadaptiveChild() ;
								neighbors[i] = nullptr ; counter-- ;	//nb+= neighbors[i]->sizeofadaptiveChild() ;
							}
						}
					}
					std::cout << "    M2L operator on cell " <<octreeIterator.getCurrentCell()->getGlobalId() <<" and set local expansion  on cell "<<localExpCell->getGlobalId()  <<
							"   "<<counter <<" cells. "		<< nb<<std::endl;
					if(nb>0 ){
						std::cout << "              WARNING need Adaptive Child in interaction list"<< std::endl;
					}
					if(!localExpCell->isSminMCriteria() ){
						for (int i = 0 ; i < 343 ; ++i){
							if(neighbors[i]){
								if( ! neighbors[i]->isSminMCriteria() ) { // Here we need Adaptive child
									std::cout << "          M2P operator on cell " <<neighbors[i]->getGlobalId()   <<" and  Leaf inside cell "<< localExpCell->getGlobalId()  <<std::endl;
								}
								else {
									std::cout << "          P2P operator on leaf inside cell " <<neighbors[i]->getGlobalId()    <<" and  leaf inside cell "<<localExpCell->getGlobalId()  <<std::endl;
								}
							}
						}
					}
					else {
						for (int i = 0 ; i < 343 ; ++i){
							if(neighbors[i]){
								if( ! neighbors[i]->isSminMCriteria() ) { // Here we need Adaptive child
									std::cout << "          M2L operator on cell " <<localExpCell->getGlobalId()   <<" and  cell "<<neighbors[i]->getGlobalId()  <<std::endl;
								}
								else {
									std::cout << "          P2L operator on cell " <<localExpCell->getGlobalId()   <<" and  leaf "<<neighbors[i]->getGlobalId()  <<std::endl;
								}
							}
						}
					}
				}
				FLOG(computationCounter.tac());
			} while(octreeIterator.moveRight());

			FLOG(computationCounter.tic());
			//			kernels->finishedLevelM2L(idxLevel);
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
		FLOG( FLog::Controller.write("\tStart Downward Pass (L2L)\n").write(FLog::Flush); );
		FLOG(FTic counterTime);
		FLOG(FTic computationCounter );
		std::cout << "-->downardPass"<<std::endl;

		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.moveDown();

		typename OctreeClass::Iterator avoidGotoLeftIterator(octreeIterator);

//		const int heightMinusOne = OctreeHeight - 1;
		// for each levels exepted leaf level
		//		for(int idxLevel = 2 ; idxLevel < heightMinusOne ; ++idxLevel ){
		for(int idxLevel = 3 ; idxLevel < OctreeHeight ; ++idxLevel ){
			FLOG(FTic counterTimeLevel);
			std::cout << "  Level "<<idxLevel<<std::endl;

			// for each cells
			do{
				FLOG(computationCounter.tic());
				if (octreeIterator.getCurrentCell()->isadaptive() ){
					const auto & parentCell = octreeIterator.getCurrentCell()->getadaptiveFather();
					if(parentCell.cell) {
						if (! octreeIterator.getCurrentCell()->isSminMCriteria() ){
							std::cout << "     L2L operator to do between " <<octreeIterator.getCurrentCell()->getGlobalId() <<" and "
									<< parentCell.cell->getGlobalId() << std::endl;
						}
						else {

							std::cout << "     L2P operator to do between " <<parentCell.cell->getGlobalId()  <<" and leaves: ";
							for (int i = 0; i <  octreeIterator.getCurrentCell()->getLeavesSize() ; ++i){
								std::cout << "  " << octreeIterator.getCurrentCell()->getLeaf(i)->getIndex();
							}
							std::cout << "  Inside cell: " << octreeIterator.getCurrentCell()->getGlobalId()<< std::endl ;
						}

						//				kernels->L2L( octreeIterator.getCurrentCell() , octreeIterator.getCurrentChild(), idxLevel);
					}
				}
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
		FLOG( FLog::Controller.write("\tStart Direct Pass\n").write(FLog::Flush); );
		FLOG(FTic counterTime);
		FLOG(FTic computationCounterL2P);
		FLOG(FTic computationCounterP2P);
		std::cout << "-->directPass"<<std::endl;

		const int heightMinusOne = OctreeHeight - 1;

		typename OctreeClass::Iterator octreeIterator(tree);
		octreeIterator.gotoBottomLeft();
		// There is a maximum of 26 neighbors
		ContainerClass* neighbors[27];
		// for each leafs
		do{
			FLOG(computationCounterL2P.tic());
			//			kernels->L2P(octreeIterator.getCurrentCell(), octreeIterator.getCurrentListTargets());
			if( !octreeIterator.getCurrentCell()->isSminMCriteria()){
				//
				std::cout << "  L2P operator on cell " <<octreeIterator.getCurrentCell()->getGlobalId() <<std::endl;
				//		kernels->P2M( octreeIterator.getCurrentCell() , octreeIterator.getCurrentListSrc());
			}
			FLOG(computationCounterL2P.tac());
			// need the current particles and neighbors particles
			const int counter = tree->getLeafsNeighbors(neighbors, octreeIterator.getCurrentGlobalCoordinate(),heightMinusOne);
			FLOG(computationCounterP2P.tic());
			//			kernels->P2P(octreeIterator.getCurrentGlobalCoordinate(),octreeIterator.getCurrentListTargets(),
			//					octreeIterator.getCurrentListSrc(), neighbors, counter);
			std::cout << "  P2P operator on leaf " <<octreeIterator.getCurrentCell()->getGlobalId() <<" and "<< counter <<" leafs. ";
			// Lister les cellules
			//			for(int ii = 0 ; ii < 27 ; ++ii){
			//				if(neighbors[ii]){
			//					std::cout << " " << neighbors[ii]->getIndex() ;
			//				}
			//			}
			std::cout << std::endl;
			FLOG(computationCounterP2P.tac());
		} while(octreeIterator.moveRight());


		FLOG( FLog::Controller << "\tFinished (@Direct Pass (L2P + P2P) = "  << counterTime.tacAndElapsed() << "s)\n" );
		FLOG( FLog::Controller << "\t\t Computation L2P : " << computationCounterL2P.cumulated() << " s\n" );
		FLOG( FLog::Controller << "\t\t Computation P2P : " << computationCounterP2P.cumulated() << " s\n" );

	}

};


#endif //FADAPTSEQALGORITHM_HPP


