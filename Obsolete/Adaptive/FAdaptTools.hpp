// ===================================================================================
// Copyright ScalFmm 2016 INRIA, Olivier Coulaud, Bérenger Bramas,
// Matthias Messner olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the
// FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
// An extension to the license is given to allow static linking of scalfmm
// inside a proprietary application (no matter its license).
// See the main license file for more details.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FADAPTTOOLS_HPP
#define FADAPTTOOLS_HPP
// Keep in private GIT
// @SCALFMM_PRIVATE
#include <fstream>

#include "Containers/FOctree.hpp"
//#include "Components/FSimpleLeaf.hpp"
#include "Components/FBasicParticleContainer.hpp"

#include "Adaptive/FAdaptiveCell.hpp"

//!  \fn   adaptiveTreeBuilSminC(OctreeClass  & tree) {


//! \brief Build father and child lists

//! \param tree Octree
//! \param sminM Criteria to use P2M rather M2M
//! \param sminL ??
//!
template  <class OctreeClass >
void   adaptiveTreeBuilSminC(OctreeClass & tree,const int sminM, const int sminL) {
	//
	typename OctreeClass::Iterator octreeIterator(&tree) ;
	int NbLevels = tree.getHeight();

    FSize  nbPart ;
	//
	octreeIterator.gotoBottomLeft();
	// Set s on the cells at leave level
	//
	//  Build SMIN criteria
	//
	do{  // leaves level = NbLevels - 1
		nbPart = octreeIterator.getCurrentListTargets()->getNbParticles() ;
		octreeIterator.getCurrentCell()->addPart(nbPart);
		octreeIterator.getCurrentCell()->setCelladaptive();
		if(nbPart <= sminM){
			octreeIterator.getCurrentCell()->setSminMCriteria();
		}
	}while(octreeIterator.moveRight());
	//
	octreeIterator.moveUp();
	octreeIterator.gotoLeft();
	//
	for(int idxLevel = NbLevels - 2 ; idxLevel >= 1 ; --idxLevel){
		int nbChildForMyCell;
		//
		do{
			// Check number of
			nbChildForMyCell=0 ;
			auto ** child      = octreeIterator.getCurrentChild();
			auto  &  myCell = *(octreeIterator.getCurrentCell());
			nbPart = 0 ; nbChildForMyCell =0;
			//						std::cout << "NB: ";
			for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
				if(child[idxChild]) {
					++nbChildForMyCell;
					nbPart += child[idxChild]->getnbPart();
					//								std::cout << "  "<< child[idxChild]->getnbPart();
				}
			}
			//						std::cout << std::endl;
			myCell.addPart(nbPart);
			if(nbChildForMyCell>1) {
				myCell.setCelladaptive(); ;
			}
			if(nbPart <= sminM){
				myCell.setSminMCriteria();
			}
		} while(octreeIterator.moveRight());
		//
		//  Go to next level
		octreeIterator.moveUp();
		octreeIterator.gotoLeft();
	}

}
//!  \fn   adaptiveTreeBuildLists(OctreeClass & tree) {


//! \brief Build father and child lists

//! \param tree Octree
//!
template  <class OctreeClass >
void   adaptiveTreeBuildLists(OctreeClass & tree) {
	//
	typename OctreeClass::Iterator octreeIterator(&tree) ;
	//
	int NbLevels = tree.getHeight();
	//
	octreeIterator.gotoBottomLeft();
	//
	// Add leaf pointer on cells at the leaves level
	do {
		octreeIterator.getCurrentCell()->addLeafptr(octreeIterator.getCurrentLeaf ()   ) ;
	}while(octreeIterator.moveRight());
	octreeIterator.gotoLeft();
	//
	for(int idxLevel = NbLevels ; idxLevel > 1 ; --idxLevel){
		//		std::cout << "idxLevel "<< idxLevel <<std::endl;
		do {
			//		FAdaptCell  &  myCell = *(octreeIterator.getCurrentCell());
			auto  &  myCell = *(octreeIterator.getCurrentCell());
			//			std::cout  << "     My ID "<< myCell.getGlobalId()
			//							<< "  isAdapt ? " << std::boolalpha<<myCell.isadaptive()<<  std::endl ;
			//
			typename OctreeClass::Iterator findParent(octreeIterator);
			// std::cout <<  "      Same cell  parent iterator "<< findParent.getCurrentCell()->getGlobalId() <<  std::endl ;
			while(findParent.moveUp() ) {
				// std::cout << "           Level " <<octreeIterator.level() <<std::endl;
				// std::cout << "              Check  cell (id, level, isAdap)  " <<  findParent.getCurrentCell()->getGlobalId()
				//				<< "    " <<findParent.level()  <<   "   "<< std::boolalpha<< findParent.getCurrentCell()->isadaptive()  <<  std::endl ;
				if ( findParent.getCurrentCell()->isadaptive() ){
					//
					//							std::cout << "                         I found " <<std::endl;
					myCell.addadaptiveParent(findParent.getCurrentCell(), findParent.level()) ;
					//					std::cout <<  "       FA      "  << findParent.getCurrentCell()->getGlobalId() << "  Level " << findParent.level()<< std::endl;

					if(myCell.isadaptive()){
						findParent.getCurrentCell()->addadaptiveChild(&myCell,idxLevel) ;
						//						std::cout <<  "       CA       "  << myCell <<std::endl;
					}
					break ;
				}
			}

		}while(octreeIterator.moveRight());
		//
		// Add links to leaves if needed
		//
		if (idxLevel != NbLevels) {
			octreeIterator.gotoLeft();
			do {
				if(octreeIterator.getCurrentCell()->isSminMCriteria() ){
					auto** child = octreeIterator.getCurrentChild();
					for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
						if(child[idxChild] ) {
							for(int i = 0 ; i < child[idxChild]->getLeavesSize()  ; ++i){
								octreeIterator.getCurrentCell()->addLeafptr(child[idxChild]->getLeaf(i)   ) ;
							}
						}
					}
				}
			}while(octreeIterator.moveRight());

		}
		//
		//  Go to next level
		octreeIterator.moveUp();
		octreeIterator.gotoLeft();
	}
};
/////////////////////////////////////////////////////////////////////////////////////////
//!  \fn   setGlobalID(FOctree< FAdaptCell,  FBasicParticleContainer<FReal,0,FReal> ,  FSimpleLeaf<FReal,  FBasicParticleContainer<FReal,0,FReal>> >  & tree) {
/////////////////////////////////////////////////////////////////////////////////////////


//! \brief Set a unique global Id to all cells into the octree

//! \param tree Octree
//! \return Return the total number of cells in the octree (including the leaves level.)
//template  <class Container, class Leaf >
template  <class OctreeClass >
long int  setGlobalID(OctreeClass  & tree) {
	//
	typename OctreeClass::Iterator octreeIterator(&tree) ;
	//
	//  Set Global id
	//
	int NbLevels = tree.getHeight();
	//
	octreeIterator.gotoTop() ;
	long int idCell =1 ;
	for(int idxLevel = 1 ; idxLevel < NbLevels ;  ++idxLevel){
		do{
			++idCell ;
			octreeIterator.getCurrentCell()->setGlobalId(idCell);
		} while(octreeIterator.moveRight());
		octreeIterator.moveDown() ;
		octreeIterator.gotoLeft();
	}
	octreeIterator.gotoLeft();
	do{
		long int id = octreeIterator.getCurrentCell()->getGlobalId();
		octreeIterator.getCurrentLeaf()->setIndex(id) ;
	} while(octreeIterator.moveRight());

	return idCell;
} ;
//! \brief export octree in Tulip format

//! \param nbCell number of cells (including the leaves level)
//!
template  <class OctreeClass >
void  TulipExport( std::ofstream& tlp, long int nbCell, OctreeClass & tree) {
	//
	tlp << "(tlp \"2.0\" "<< std::endl
			<< "(date \"04-11-2014\") "<< std::endl
			<< "(comments \"Octree \")"<< std::endl
			<<"(nodes ";

	int NbLevels = tree.getHeight();
	typename OctreeClass::Iterator octreeIterator(&tree) ;

	//
	for(int i = 1 ; i <= nbCell ;  ++i){
		tlp << "  " << i ;
	}
	tlp << " )" << std::endl
			<< ";(edge <edge_id> <source_id> <target_id>)" <<std ::endl;
	int edge = 0 ;
	// edge with the root
	octreeIterator.gotoTop() ;
	bool isadaptiveRoot = false;
	// Build edges between root and its child.
	do{
		++edge ;
		tlp <<"(edge " << edge << " 1 " <<octreeIterator.getCurrentCell()->getGlobalId() << " )" <<std::endl;
	} while(octreeIterator.moveRight());
	if (edge > 1) {isadaptiveRoot=true;}
	// End  for level 0
	//
	// Level 1 return to write
	octreeIterator.gotoLeft();
	for(int idxLevel = 1 ; idxLevel < NbLevels-1 ;  ++idxLevel){
		do{
			long int  myId= octreeIterator.getCurrentCell()->getGlobalId();

			auto** child = octreeIterator.getCurrentChild();
			//						std::cout << "NB: ";
			for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
				if(child[idxChild]) {
					++edge  ;
					tlp <<"(edge " << edge << "  " << myId << "  "  << child[idxChild]->getGlobalId() << " )" <<std::endl;
				}
			}
		} while(octreeIterator.moveRight());
		octreeIterator.moveDown() ;
		octreeIterator.gotoLeft();
	}
	//
	// Set adaptive property
	//
	std::string rootColor("  \"(255,0,0,255)\" ") ;      // Bleu
	std::string adaptiveColor("  \"(0,0,255,255)\" ") ;      // Bleu
	std::string defaultEdgeColor("  \"(255,191,0,255)\" ") ;
	std::string defaultNodeColor("  \"(0,0,0,255)\" ") ;    // Black
	std::string leafColor("  \"( 0,255,0,255)\" ") ;
	tlp << "(property  0 color \"viewColor\" " <<std::endl
			<< "(default " << defaultNodeColor <<defaultEdgeColor << "  ) " << std::endl ;
	octreeIterator.gotoTop() ;
	if (isadaptiveRoot) {
		tlp << "  (node 1 "  <<rootColor<<"  )" << std::endl ;
	}
	for(int idxLevel = 1 ; idxLevel < NbLevels-1 ;  ++idxLevel){
		do{
			if(octreeIterator.getCurrentCell()->isadaptive()){
				tlp << "  (node " <<  octreeIterator.getCurrentCell()->getGlobalId()   <<adaptiveColor<<"  )" << std::endl ;;
			}
		} while(octreeIterator.moveRight());
		octreeIterator.moveDown() ;
		octreeIterator.gotoLeft();
	}
	//	Set leaf color
	std::cout  << "  Iterator Level    " << octreeIterator.level()<<  "  is leaves level: " << octreeIterator.isAtLeafLevel()	<<std::endl;

	do{
		tlp << "  (node " <<  octreeIterator.getCurrentCell()->getGlobalId()   <<leafColor<<"  )" << std::endl ;;
	} while(octreeIterator.moveRight());
	tlp << ")" << std::endl;

	//
	// Set P2M property
	//
	std::cout <<std::endl<< "   TULIP OUTPUT PROPERTY P2M "<<std::endl<<std::endl;
	tlp << "(property  0 bool \"P2M\" " <<std::endl
			<< "(default \"false\" \"false\"  ) " << std::endl ;
	octreeIterator.gotoBottomLeft();		octreeIterator.moveUp() ;
	bool endProp = true;
	int elt =0;
	for(int idxLevel =  NbLevels-2 ; idxLevel >=1 ;  --idxLevel){
		endProp = true ;
		std::cout << " TLP " << idxLevel << "   " << octreeIterator.level() << "  "  << elt << std::endl;
		do{
			auto ** child      = octreeIterator.getCurrentChild();
			for(int idxChild = 0 ; idxChild < 8 ; ++idxChild){
				if(child[idxChild]) {
					if(child[idxChild]->isSminMCriteria()){
						endProp = false ; ++elt;
						tlp << "  (node " <<  child[idxChild]->getGlobalId()   <<"  \"" <<std::boolalpha <<child[idxChild]->isSminMCriteria()<<"\"  )" << std::endl ;;
					}
				}
			}
		} while(octreeIterator.moveRight());
		if(endProp){
			break ; // all upper levels are with smin > sminM
		}
		octreeIterator.gotoLeft();
		octreeIterator.moveUp() ;
	}
	if(elt != 0)
	{tlp << ")" << std::endl;}


	// END TLP FILE
	tlp << ")" << std::endl;
};
#endif
