#ifndef FADAPTATIVECELL
#define FADAPTATIVECELL

#include "../Components/FBasicCell.hpp"
#include "../Containers/FVector.hpp"

/**
 * This class is a wrapper to work with adaptative kernel
 * It contains a pointer to the real computation cell but use it only if
 * the cell is adaptative AND has been use for development.
 *
 * A cell is adaptative if:
 * × it has more than one child
 * × it used for development
 * × it is usde to store the leaves
 * Else it stores a pointer to the lower adaptive cell.
 */
template <class RealCell, class ContainerClass>
class FAdaptativeCell : public FBasicCell {
    /** The cell used for the computation */
    RealCell* realCell;
    /** To keep track of the cell state */
    bool IamAdaptative;

    /** If not adaptative then we need to know the lower adaptative cell */
    FAdaptativeCell<RealCell, ContainerClass>* subAdaptativeCell;
    /** The lower adaptative cell level */
    int subAdaptativeLevel;

    /** The leaves that have been skiped for the P2M/M2M... */
    FVector<ContainerClass*> subLeaves;

public:
    /** Set has not adaptative by default */
    FAdaptativeCell() : realCell(nullptr), IamAdaptative(false), subAdaptativeCell(nullptr), subAdaptativeLevel(0){
    }

    ~FAdaptativeCell(){
        delete realCell;
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
            subAdaptativeCell  = nullptr;
            subAdaptativeLevel = 0;
            subLeaves.clear();
        }
        else if(!inHaveDevelopment && realCell){
            // clean real cell if needed
            delete realCell;
            realCell = nullptr;
        }
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
    /// Set adaptative
    ////////////////////////////////////////////////////////////////////////////////

    bool isAdaptative() const {
        return IamAdaptative;
    }

    void setAdaptative(const bool inIsAdaptative) {
        IamAdaptative = inIsAdaptative;
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
        return subLeaves.getSize();
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

    void setSubAdaptativeCell(FAdaptativeCell<RealCell,ContainerClass>* inSubAdaptativeCell, const int inSubAdaptativeLevel){
        subAdaptativeCell  = inSubAdaptativeCell;
        subAdaptativeLevel = inSubAdaptativeLevel;
    }

    void setSubAdaptativeCell(const FAdaptativeCell<RealCell,ContainerClass>* inSubAdaptativeCell, const int inSubAdaptativeLevel){
        subAdaptativeCell  = const_cast<FAdaptativeCell<RealCell,ContainerClass>*>(inSubAdaptativeCell);
        subAdaptativeLevel = inSubAdaptativeLevel;
    }

    FAdaptativeCell<RealCell,ContainerClass>* getSubAdaptativeCell() {
        return subAdaptativeCell;
    }

    FAdaptativeCell<RealCell,ContainerClass>* getSubAdaptativeCell() const {
        return subAdaptativeCell;
    }

    int getSubAdaptativeLevel() const {
        return subAdaptativeLevel;
    }
};

#endif
