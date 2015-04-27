// ==== CMAKE ====
// Keep in private GIT
// @SCALFMM_PRIVATE

#ifndef _FCOSTCELL_HPP_
#define _FCOSTCELL_HPP_

#include <type_traits>


/** 
 * \brief Empty trait class.
 * \author Quentin Khan
 *
 * This class is used to check whether a cell class has FCostCell in its
 * inheritance tree.
 */
class FCostCellTypeTrait {};


/**
 * \brief Cell with a cost memory for balance computations.
 * \author Quentin Khan
 *
 * This class extends BaseClass to add simple computation cost memory.
 *
 * \tparam BaseClass The base cell class to extend. The constructors are
 * inherited automatically.
 * \tparam CostType The type to use in order to store the cost. Defaults to FSize.
 */
template<typename BaseClass, typename CostType = FSize>
class FCostCell : public BaseClass, virtual public FCostCellTypeTrait {
    static_assert(std::is_arithmetic<CostType>::value,
                  "The cell cost type must be an arithmetic type.");

    /// The far-field cost of the cell.
    /** Declared mutable because the existing algorithms use const cells.*/
    mutable CostType _cost = 0;

    /// The near-field cost of the cell.
    /** Declared mutable because the existing algorithms use const cells.*/
    mutable CostType _leafCost = 0;

public:
    /// Type definition that can be retrieved by other classes
    using costtype = CostType;

    using BaseClass::BaseClass;

    /// Debug member, used to check whether the cell was already visited.
    bool _visited = false;

    /**
     * \brief Gets the far-field cost of the cell.
     * \return The far-field cost of the cell
     */
    CostType getCost() const {
        return _cost;
    }

    /**
     * \brief Gets the near-field cost of the cell.
     * \return The near-field cost of the cell
     */
    CostType getNearCost() const {
        return _leafCost;
    }

    /**
     * \brief Sets the cost of the cell.
     */
    void setCost(CostType newCost) {
        _cost = newCost;
    }

    /**
     * \brief Sets the near-field cost of the cell.
     */
    void setNearCost(CostType newCost) {
        _leafCost = newCost;
    }

    /** 
     * \brief Add a far-field cost to the cell.
     * \return The cost of the cell
     * \warning Can be used on const cells !
     */
    CostType addCost(CostType cost) const {
        _cost += cost;
        return _cost;
    }

    /** 
     * \brief Add a near-field cost to the cell.
     * \return The cost of the cell
     * \warning Can be used on const cells !
     */
    CostType addNearCost(CostType cost) const {
        _leafCost += cost;
        return _leafCost;
    }
};

#endif
