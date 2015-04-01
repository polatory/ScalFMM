// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _FCOSTCELL_HPP_
#define _FCOSTCELL_HPP_

/**
 * \brief Cell with a cost memory for balance computations.
 * \author Quentin Khan
 *
 * This class extends BaseClass to add simple computation cost memory.
 *
 * \tparam BaseClass The base cell class to extend. The constructors are
 * inherited automatically.
 * \tparam The type to use in order to store the cost. Defaults to FSize.
 */
template<typename BaseClass, typename CostType = FSize>
class FCostCell : public BaseClass {
    /// The cost of the cell. It is declared mutable because the existing
    /// algorithms use const cells.
    mutable CostType _cost = 0;
public:
    using BaseClass::BaseClass;

    /// Debug member, used to check whether the cell was already visited.
    bool _visited = false;

    /**
     * \brief Gets the cost of the cell.
     * \return The cost of the cell
     */
    CostType getCost() const {
        return _cost;
    }

    /**
     * \brief Sets the cost of the cell.
     */
    void setCost(CostType newCost) {
        _cost = newCost;
    }

    /** 
     * \brief Add a cost to the cell.
     * \return The cost of the cell
     * \warning Can be used on const cells !
     */
    CostType addCost(CostType cost) const {
        _cost += cost;
        return _cost;
    }
};

#endif
