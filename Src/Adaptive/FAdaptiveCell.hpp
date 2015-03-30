#ifndef FAdaptiveCELL
#define FAdaptiveCELL
// Keep in private GIT
// @SCALFMM_PRIVATE
#include "Components/FBasicCell.hpp"

/**
 * This class is a wrapper to work with adaptive kernel
 * It contains a pointer to the real computation cell but use it only if
 * the cell is adaptive AND has been use for development.
 *
 * A cell is adaptive if:
 *  - it has more than one child
 *  -  it used for development
 *  -  it is use to store the leaves
 * Else it stores a pointer to the lower adaptive cell.
 */
template <class RealCell, class ContainerClass>
class FAdaptiveCell : public FBasicCell {
	/** The cell used for the computation */
	RealCell* realCell;
	/** To keep track of the cell state */
	bool IamAdaptive;

	/** If not adaptive then we need to know the lower adaptive cell */
	FAdaptiveCell<RealCell, ContainerClass>* subAdaptiveCell;
	/** The lower adaptive cell level */
	int subAdaptiveLevel;

	/** The leaves that have been skiped for the P2M/M2M... */
	FVector<ContainerClass*> subLeaves;
	//
	// Global Index of the cell in the octree (This id is unique)
    FSize gID;

public:
	/** Set has not Adaptive by default */
	FAdaptiveCell() : realCell(nullptr), IamAdaptive(false), subAdaptiveCell(nullptr), subAdaptiveLevel(0),gID(0){
	}

	~FAdaptiveCell(){
		if (realCell) delete realCell;
	}

	void resetToInitialState(){
		subLeaves.clear();
		if(realCell){
			realCell->resetToInitialState();
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	/// Get the real cell
	////////////////////////////////////////////////////////////////////////////////

	/** To say if it is used for development or not */
	void setHaveDevelopment(const bool inHaveDevelopment) {
		if(inHaveDevelopment && !realCell){
			// alloc and init the real cell
			realCell = new RealCell;
			realCell->setMortonIndex(this->getMortonIndex());
			realCell->setCoordinate(this->getCoordinate());
			// clean other information
			subAdaptiveCell  = nullptr;
			subAdaptiveLevel = 0;
			subLeaves.clear();
		}
		else if(!inHaveDevelopment && realCell){
			// clean real cell if needed
			delete realCell;
			realCell = nullptr;
		}
	}

    void resetSubLeaves() {
        subLeaves.clear();
    }

	bool hasDevelopment() const{
		return realCell != nullptr;
	}

	RealCell* getRealCell(){
		return realCell;
	}

	const RealCell* getRealCell() const {
		return realCell;
	}

	////////////////////////////////////////////////////////////////////////////////
	/// Set Adaptive
	////////////////////////////////////////////////////////////////////////////////

	bool isAdaptive() const {
		return IamAdaptive;
	}

	void setAdaptive(const bool inIsAdaptive) {
		IamAdaptive = inIsAdaptive;
	}

	////////////////////////////////////////////////////////////////////////////////
	/// Manage the sub leaves
	////////////////////////////////////////////////////////////////////////////////

	void addSubLeaf(const ContainerClass* aLeaf){
		subLeaves.push(const_cast<ContainerClass*>(aLeaf));
	}

	void addSubLeaf(ContainerClass* aLeaf){
		subLeaves.push(aLeaf);
	}

	void addSubLeaves(const ContainerClass*const* aLeavesToInsert, const int nbLeavesToInsert){
		subLeaves.memocopy(const_cast<ContainerClass*const*>(aLeavesToInsert),nbLeavesToInsert);
	}

	int getNbSubLeaves() const {
        return int(subLeaves.getSize());
	}

	ContainerClass* const * getSubLeaves() {
		return subLeaves.data();
	}

	const ContainerClass * const * getSubLeaves() const{
		return subLeaves.data();
	}

	ContainerClass* getSubLeaf(const int leafIdx) const{
		return subLeaves[leafIdx];
	}

	////////////////////////////////////////////////////////////////////////////////
	/// Manage the sub cell
	////////////////////////////////////////////////////////////////////////////////

	void setSubAdaptiveCell(FAdaptiveCell<RealCell,ContainerClass>* inSubAdaptiveCell, const int inSubAdaptiveLevel){
		subAdaptiveCell  = inSubAdaptiveCell;
		subAdaptiveLevel = inSubAdaptiveLevel;
	}

	void setSubAdaptiveCell(const FAdaptiveCell<RealCell,ContainerClass>* inSubAdaptiveCell, const int inSubAdaptiveLevel){
		subAdaptiveCell  = const_cast<FAdaptiveCell<RealCell,ContainerClass>*>(inSubAdaptiveCell);
		subAdaptiveLevel = inSubAdaptiveLevel;
	}

	FAdaptiveCell<RealCell,ContainerClass>* getSubAdaptiveCell() {
		return subAdaptiveCell;
	}

	FAdaptiveCell<RealCell,ContainerClass>* getSubAdaptiveCell() const {
		return subAdaptiveCell;
	}

	int getSubAdaptiveLevel() const {
		return subAdaptiveLevel;
	}
	//#ifdef DEBUG_ADAPTIVE
	////////////////////////////////////////////////////////////////////////////////
	/// Manage a global IG DEBUG PURPOSE
	////////////////////////////////////////////////////////////////////////////////
	//! Return the global Id of the cell in the octree
    const FSize  getGlobalId(){
		return this->gID ;
	}
    const FSize  getGlobalId( ) const{
		return this->gID ;
	}
	//!  Set he global Id of the cell in the octree to id
    void setGlobalId(const FSize & id){
		this->gID = id;  ;
	}
	//#endif
	/**
	 * Operator stream FAdaptCell to std::ostream
	 *
	 * @param[in,out] output where to write the adaptive cell
	 * @param[in] inPosition the cell to write out
	 * @return the output for multiple << operators
	 */
	template <class StreamClass>
	friend StreamClass& operator<<(StreamClass& output, const FAdaptiveCell<RealCell,  ContainerClass>&  cell){
		output << "(  Cell Id " << cell.getGlobalId()  << " Adaptive  " <<  std::boolalpha << cell.isAdaptive()
								<< "  sminM " << " cell.isSminMCriteria() "<< " "<< cell.getNbSubLeaves()  ;
		if(cell.getNbSubLeaves() >0){
			output << " LF = { " ;
			for (int i=0; i	 <cell.getNbSubLeaves()  ; ++i){
				output <<cell.getSubLeaf(i)->getNbParticles() << " ";
			}
			output << "}" ;
		}
		//		output <<" CA={ ";
		//		const FVector<FExtACell>  * v =cell.getadaptiveChild() ;
		//		if (cell.sizeofadaptiveChild()> 0 ){
		//			for (int i=0; i < v->getSize() ; ++i){
		//				output << v->operator [](i).cell->getGlobalId() << " ";
		//			}
		//		}
		//		output << "} " ;
		//		if(cell.getadaptiveFather().cell){
		//			output << " FA={" << (cell.getadaptiveFather()).cell->getGlobalId() << "} " ;
		//		}
		//		else
		//		{
		//			output <<  "   FA={} " ;
		//		}
		if(cell.hasDevelopment()){
			output <<*(cell.getRealCell()) ;
		}
		//	output << std::endl << "    Multipoles " << *(cell.getRealCell()) << std::endl;
		output << " )"	<<std::endl;
		return output;  // for multiple << operators.
	}
};

#endif
