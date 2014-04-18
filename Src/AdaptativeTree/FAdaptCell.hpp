// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas
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
#ifndef FADAPTCELL_HPP
#define FADAPTCELL_HPP

#include <cstddef>
#include <vector>
//
#include "Components/FBasicCell.hpp"
#include "Containers/FVector.hpp"

/**
 * @author Olivier Coulaud (Olivier.Coulaud@inria.fr)
 * @class FAdaptCell*
 * @brief This class defines adaptative cell.
 *
 *  ToDO
 *
 */
//class FAdaptCell;

template <class CellClass,  class LeafClass>
class FAdaptCell : public FBasicCell  {
public:
	struct  FExtACell {
		FAdaptCell<CellClass, LeafClass> * cell ;
		int               level  ;  //< Level in the octree of cell
		FExtACell() { cell=nullptr ; level =-1 ; } ;
	};
protected:
	// Global Index of the cell in the octree (This id is unique)
	long int gID;
	//! Number of particles inside the cell
	int nbP ;
	//! iAdaptive cell True, false
	bool adaptative ;
	bool sminMCriteria;
	//
	//  Lists
	//
	FExtACell                  adaptativeParent ; //< Adaptative Parent of the cell
	FVector<FExtACell> adaptativeChild ;   //< list of adaptative child of the cell
	FVector<LeafClass*> leaves ;                   //< list of leaf child of the cell
	//
	//
	CellClass   * trueFMMCell ;                        //<a pointer on the cell that contains Multipole and local values
public:
	FAdaptCell(): gID(-1), nbP(0), adaptative(false),trueFMMCell(nullptr) {
	}
	/** Default destructor */
	virtual ~FAdaptCell(){
	}
	//! add nbPart in cell
	void addPart(const int n){
		this->nbP += n ;
	}
	//! Return the number of particles inside the cell
	const int getnbPart() const{
		return this->nbP ;
	}
	//!return true if the sminM criteria is satisfied
	bool isSminMCriteria (){
		return this->sminMCriteria;
	}
	//!return true if the sminM criteria is satisfied
	bool isSminMCriteria () const{
		return this->sminMCriteria;
	}
	//! Set if the cell is adaptative or not
		void setSminMCriteria(const bool bb=true){
			this->sminMCriteria=bb ;
		}

	//! Return the global Id of the cell in the octree
	const long int  getGlobalId(){
		return this->gID ;
	}
	const long int  getGlobalId( ) const{
		return this->gID ;
	}
	//!  Set he global Id of the cell in the octree to id
	void setGlobalId(const long int & id){
		this->gID = id;  ;
	}
	//! Set if the cell is adaptative or not
	void setCellAdaptative(const bool bb=true){
		this->adaptative = bb;
	}
	//!return true if the cell is adaptative
	const bool isAdaptative() const{
		return this->adaptative;
	}
	//
	//! Add the adaptive parent of the cell
	void addAdaptativeParent( const FExtACell &ac) {
		this->adaptativeParent.cell   = ac.cell ;
		this->adaptativeParent.level = ac.level ;
		std::cout <<  "          addAdaptativeParent      "  << *(this) <<std::endl;
	}
	//
	//! Add the adaptive child of the cell
	void addAdaptativeParent( FAdaptCell<CellClass, LeafClass> * cell ,const int idxLevel)  {
		this->adaptativeParent.cell = cell ; 	this->adaptativeParent.level = idxLevel ;

	}
	FExtACell getAdaptativeFather() {
		return this->adaptativeParent ;
	}
	FExtACell getAdaptativeFather() const {
		return this->adaptativeParent ;
	}
	//
	//! Add the adaptive child of the cell
	void addAdaptativeChild( FAdaptCell<CellClass, LeafClass> * cell ,const int idxLevel)  {
		FExtACell AC ;
		AC.cell = cell ; 	AC.level = idxLevel ;
		this->adaptativeChild.push(AC);
	}
	//
	//! Add the adaptive child of the cell
	void addAdaptativeChild( const  FExtACell &AC){
		this->adaptativeChild.push(AC) ;
	}
	int sizeofAdaptativeChild() const{
		return this->adaptativeChild.getSize();
	}
	FVector<FExtACell>* getAdaptativeChild(){
		return &this->adaptativeChild;
	}
	const  FVector<FExtACell>* getAdaptativeChild() const{
		return &this->adaptativeChild;
	}
	//
	//! Add the adaptive child of the cell
	void addLeafptr( LeafClass * leaf)  {
		this->leaves.push(leaf) ;
	}
	//! Return the number of leaves
	int getLeavesSize() const{
		return this->leaves.getSize();
	}
	//! Return the number of leaves
	LeafClass* getLeaf(const int i ) {
		return this->leaves[i];
	}
	//! Return the number of leaves
	LeafClass* getLeaf(const int i ) const {
		return this->leaves[i];
	}
	/** Make it like the beginning */
	void resetToInitialState(){
		//this->dataDown = 0;
		//this->dataUp = 0;
		this->nbP = 0 ;
	}

	/////////////////////////////////////////////////

	/** Save the current cell in a buffer */
	template <class BufferWriterClass>
	void save(BufferWriterClass& buffer) const{
		FBasicCell::save(buffer);
		//buffer << dataDown << dataUp << nbP;
	}

	/** Restore the current cell from a buffer */
	template <class BufferReaderClass>
	void restore(BufferReaderClass& buffer){
		FBasicCell::restore(buffer);
	//	buffer >> dataDown >> dataUp >> nbP;
	}



	/////////////////////////////////////////////////

	/** Serialize only up data in a buffer */
	template <class BufferWriterClass>
	void serializeUp(BufferWriterClass& buffer) const {
		//buffer << this->dataUp;
	}
	/** Deserialize only up data in a buffer */
	template <class BufferReaderClass>
	void deserializeUp(BufferReaderClass& buffer){
		//buffer >> this->dataUp;
	}

	/** Serialize only down data in a buffer */
	template <class BufferWriterClass>
	void serializeDown(BufferWriterClass& buffer) const {
		//buffer << this->dataDown;
	}
	/** Deserialize only up data in a buffer */
	template <class BufferReaderClass>
	void deserializeDown(BufferReaderClass& buffer){
		//buffer >> this->dataDown;
	}
	/**
	 * Operator stream FAdaptCell to std::ostream
	 *
	 * @param[in,out] output where to write the adaptative cell
	 * @param[in] inPosition the cell to write out
	 * @return the output for multiple << operators
	 */
	template <class StreamClass>
	friend StreamClass& operator<<(StreamClass& output, const FAdaptCell<CellClass,LeafClass>&  cell){
		output << "(  Cell Id " << cell.getGlobalId()  << " Adpatative  " << std::boolalpha << cell.isAdaptative()
						<< "  sminM " << cell.isSminMCriteria()<< " "<< cell.getnbPart()  ;
		if(cell.getLeavesSize() >0){
			output << " LF={" ;
			for (int i=0; i	 <cell.getLeavesSize()  ; ++i){
				output << cell.getLeaf(i) << " ";
			}
			output << "}" ;
		}
		output <<" CA={ ";
		const FVector<FExtACell>  * v =cell.getAdaptativeChild() ;
		if (cell.sizeofAdaptativeChild()> 0 ){
			for (int i=0; i < v->getSize() ; ++i){
				output << v->operator [](i).cell->getGlobalId() << " ";
			}
		}
		output << "} " ;
		if(cell.getAdaptativeFather().cell){
			output << " FA={" << (cell.getAdaptativeFather()).cell->getGlobalId() << "} " ;
		}
		else
		{
			output <<  "   FA={} " ;
		}

		output << " )"	<<std::endl;
		return output;  // for multiple << operators.
	}

};


#endif //FADAPTCELL_HPP


