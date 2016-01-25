#ifndef SCALFMM_NODE_HPP_
#define SCALFMM_NODE_HPP_

#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_set>

#include "FBox.hpp"
#include "Utils/FPoint.hpp"
#include "Utils/FConstFuncs.hpp"

namespace scalfmm {
    namespace tests {
        struct test_Node;
        struct test_NodeIterator;
        struct test_InOrderNodeIterator;
        struct test_PreOrderNodeIterator;
        struct test_PostOrderNodeIterator;
    }
}

struct NodeEmptyData {};

/** \brief Tree node implementation
 *
 */
template<class _Tree, class _ParticleContainer, class NodeData>
class FNode {
public:
    /// Tree type this class belongs to
    using tree_t = _Tree;
    /// Type used to represent real numbers
    using Real = typename tree_t::Real;
    /// Space dimension count
    constexpr static const std::size_t Dim = _Tree::Dim;
    /// Child count if the node is not a leaf
    constexpr static std::size_t child_count = Fpow(2,Dim);

    /// Node position type
    using position_t = FPoint<Real, Dim>;
    /// Node bounding box type
    using box_t =  FBox<position_t>;
    /// Interaction lists type
    using interaction_list_t = std::unordered_set<FNode*>;
    /// Children array type
    using child_node_array_t = std::array<FNode*, child_count>;
    /// Particle container type
    using particle_container_t = _ParticleContainer;
    /// Particle type
    using particle_t = typename particle_container_t::value_type;
    /// Node data structure
    using data_t = NodeData;


private:

    friend struct scalfmm::tests::test_Node;
    friend struct scalfmm::tests::test_NodeIterator;
    friend struct scalfmm::tests::test_InOrderNodeIterator;
    friend struct scalfmm::tests::test_PreOrderNodeIterator;
    friend struct scalfmm::tests::test_PostOrderNodeIterator;
    friend tree_t;

    /// Children array, filled with nullptr is Node is a leaf
    child_node_array_t _children{{}};
    /// Node parent, nullptr if the node is a root
    FNode* _parent = nullptr;
    /// Node spatial bounding box
    box_t _box;
    /// Node depth in its tree
    std::size_t _depth = 0;
    /// Node index in parent child array
    std::size_t _m_idx = 0;
    /// Particle container
    std::unique_ptr<particle_container_t> _p_container{new particle_container_t};
    /// Node data
    data_t _data = data_t();

    /// Tree the Node belongs to
    tree_t* _tree;

    /// Indicates whether node is a leaf
    bool _is_leaf = true;

public:
    /// Near-field leaves interaction list
    interaction_list_t U;
    /// Mid-field node interaction list
    interaction_list_t V;
    /// Near-field node interaction list
    interaction_list_t W;
    /// Mid-field leaf interaction list
    interaction_list_t X;

    /** Constructor called from parent
     *
     * \param parent The parent node
     * \param child_index The index of this node in the parent children array
     */
    FNode(FNode& parent, const std::size_t& child_index) :
        _parent(&parent),
        _box   (parent.getBox().center(), parent.getBox().corner(child_index) ),
        _depth (parent.getDepth()+1),
        _m_idx((parent.getIndex() << Dim) + child_index),
        _tree  (parent._tree )
    {
        if (child_index >= child_count) {
            throw std::invalid_argument(std::string("Wrong child index in node contructor: got ")
                                        + std::to_string(child_index)
                                        + std::string(", expected at most ")
                                        + std::to_string(Dim)
                                        + ".");
        }
    }

    /** Root constructor called from tree */
    FNode(tree_t* tree) :
        _box(tree->box()),
        _tree(tree)
    {
        tree->leaves().insert(this);
    }

    /** Default constructor */
    FNode() = delete;

    /** Copy constructor */
    FNode(const FNode& other) = delete;

    /** Copy operator */
    FNode& operator=(const FNode& other) = delete;

    /** Move constructor */
    FNode(FNode&& other) = delete;

    /** Move operator */
    FNode& operator=(FNode&& other) = delete;

    /** Destructor */
    ~FNode() {
        for(auto&& child : _children) {
            delete child;
        }
    }


    /// Data accessor
    data_t& getData() noexcept {
        return _data;
    }
    /// Data const accessor
    const data_t& getData() const noexcept {
        return _data;
    }

    /// Children container accessor
    child_node_array_t& getChildren() noexcept {
        return _children;
    }
    /// Children const container accessor
    const child_node_array_t& getChildren() const noexcept {
        return _children;
    }

    /** Child container accessor
     *
     * \param index Child index
     */
    FNode* getChild(const std::size_t& index) noexcept {
        return getChildren()[index];
    }
    /** Child container const accessor
     *
     * \param index Child index
     */
    const FNode* getChild(const std::size_t& index) const noexcept {
        return getChildren()[index];
    }

    /// Parent accessor
    FNode* getParent() noexcept {
        return _parent;
    }
    /// Parent const accessor
    const FNode* getParent() const noexcept {
        return _parent;
    }

    /// Depth accessor
    const std::size_t& getDepth() const noexcept {
        return _depth;
    }

    /// Morton index accessor
    const std::size_t& getIndex() const noexcept {
        return _m_idx;
    }

    /// Tree accessor
    tree_t& getTree() noexcept {
        return *_tree;
    }
    /// Tree const accessor
    const tree_t& getTree() const noexcept {
        return *_tree;
    }

    /// Box const accessor
    const box_t& getBox() const noexcept {
        return _box;
    }

    /// Particle container accessor
    particle_container_t* getParticleContainer() noexcept {
        return _p_container.get();
    }
    /// Particle container accessor
    const particle_container_t* getParticleContainer() const noexcept {
        return _p_container.get();
    }

    /// Particle count for the container
    std::size_t getParticleCount() const noexcept {
        if(getParticleContainer()) {
            return getParticleContainer()->size();
        }
        return 0;
    }

    /** Returns true if this node and the 'other' node are adjacent
     *
     * The nodes are assumed to belong to the same tree.
     *
     * To check whether nodes are adjacent, on each axis, the distance between
     * the nodes' center is compared to the sum of their half diameter. For at
     * least one of the axes, the two must be equal. For the others, the
     * distance must be less than or equal to the sum. This ensures that a node
     * is not adjacent to one of its descendants.
     *
     * \param other The node to test adjacency with.
     *
     * \return  true if this FNode and the 'other' FNode are adjacent.
     */
    bool is_adjacent(const FNode& other) const noexcept {
        // Sum of the half side lengh of the two nodes boxes.
        // Boxes are cubes, we only need one side.
        Real centers_distance = getBox().center()[0] - getBox().c1[0]
            + other.getBox().center()[0] - other.getBox().c1[0];
        // Used to check that the other box isn't overlapping with this box
        bool one_axis_is_at_exact_distance = false;

        position_t my_center = getBox().center();
        position_t other_center = other.getBox().center();

        for(std::size_t i = 0; i < Dim; ++i) {
            Real distance = fabs(my_center[i] - other_center[i]);
            if( Ffeq(distance, centers_distance) ) {
                one_axis_is_at_exact_distance = true;
            } else if(distance > centers_distance) {
                return false;
            }
        }

        return one_axis_is_at_exact_distance;
    }

    /** Tests whether this node and the 'other' node are adjacent
     *
     * \return true if the nodes are adjacent.
     */
    bool is_adjacent(const FNode* other) const noexcept {
        if(nullptr == other) {
            return false;
        } else if (this == other){
            return true;
        }

        return is_adjacent(*other);
    }

    /** Tests whether this node is a leaf
     *
     * \return true if this node is a leaf.
     */
    bool is_leaf() const noexcept {
        return _is_leaf;
    }

    /** Inserts a particle in the node
     *
     * Inserts a particle in the node if it is a leaf. If if isn't, the
     * particle is forwarded to the relevant child.
     */
    void insert(const particle_t& p) {
        if(! is_leaf()) {
            getChildren()[box_t::space_filling_curve_t::index(p.position(), getBox().center())]->insert(p);
        } else {
            getParticleContainer()->push_back(p);
            if(getParticleContainer()->size() > getTree().leaf_max_particle_count()) {
                split();
            }
        }
    }

    /** Extracts a particle from a leaf
     *
     * If the nodes is not a leaf, does nothing.
     */
    particle_t extract(const std::size_t& idx) {
        particle_t p = *(getParticleContainer()->begin() + idx);
        if(getParticleContainer()) {
            getParticleContainer()->erase(getParticleContainer()->begin() + idx);
        }
        return p;
    }

    /** Splits or fuses nodes after tree modifications
     *
     * Iterates over all the sub-tree particles and splits the leafs that
     * hold too many particles or fuses those that hold too few.
     *
     * This may be necessary when the tree parameters have been changed
     * (leaf max particle count modified) or when particles have moved and
     * been reinserted in the tree to be in the right boxes.
     */
    void reshape() {
        if(is_leaf()) {
            if(getParticleCount() > getTree().leaf_max_particle_count()) {
                split();
            }
        } else {
            std::size_t p_count = 0;
            for(auto child : getChildren()) {
                child->reshape();
                p_count += child->getParticleCount();
            }

            if(p_count <= getTree().leaf_max_particle_count()) {
                fuse();
            }
        }
    }

    /** Applies a function to the node and it descendants in order */
    void for_each_in_order(std::function<void(FNode*)> lambda) {
        std::size_t idx = 0;
        if(! is_leaf()) {
            for(; idx < getChildren().size()/2; ++idx) {
                getChildren()[idx]->for_each_in_order(lambda);
            }
        }
        lambda(this);
        if(! is_leaf()) {
            for(; idx < getChildren().size(); ++idx) {
                getChildren()[idx]->for_each_in_order(lambda);
            }
        }
    }

    /** Applies a function to the node and it descendants post order (children first) */
    void for_each_post_order(std::function<void(FNode*)> lambda) {
        if(! is_leaf()) {
            for(std::size_t idx = 0; idx < getChildren().size(); ++idx) {
                getChildren()[idx]->for_each_post_order(lambda);
            }
        }
        lambda(this);
    }

    /** Applies a function to the node and it descendants pre order (parent first) */
    void for_each_pre_order(std::function<void(FNode*)> lambda) {
        lambda(this);
        if(! is_leaf()) {
            for(std::size_t idx = 0; idx < getChildren().size(); ++idx) {
                getChildren()[idx]->for_each_pre_order(lambda);
            }
        }
    }

    /** Equality test operator */
    bool operator==(const FNode& other) const {
        return other.getParent() == getParent()
            && other.getDepth() == getDepth()
            && other.getIndex() == getIndex();
    }

private:
    /** Tree setter */
    void setTree(tree_t* t) {
        _tree = t;
        if(! is_leaf()) {
            for(FNode*& child : getChildren()) {
                child->setTree(t);
            }
        }
    }


    /** Creates or allocates the particle container */
    void create_particle_container() {
        _p_container.reset(new particle_container_t);
    }

    /** Deletes the particle container to save space when it is not needed */
    void delete_particle_container() {
        _p_container.reset(nullptr);
    }

    /** Allocates this node's children */
    void create_children() {
        std::size_t idx = 0;

        getTree().leaves().erase(this);
        for(FNode*& child : getChildren()) {
            child = new FNode(*this, idx);
            getTree().leaves().insert(child);
            ++idx;
        }
        _is_leaf = false;
    }

    /* Deletes this node's children */
    void delete_children() {
        for(FNode*& child : getChildren()) {
            getTree().leaves().erase(child);
            delete child;
            child = nullptr;
        }
        getTree().leaves().insert(this);
        _is_leaf = true;
    }

    /** Adds children to a leaf node
     *
     * Adds children to a leaf node and redistributes its particles among
     * the newly created leaves.
     */
    void split() {
        if(getDepth()+1 > getTree().max_height()) {
            // TODO: log that there were too many particles
            return;
        }

        create_children();

        // Update interaction lists
        this->U.erase(this);
        for(FNode* child : getChildren()) {
            // Children are adjacent to each other
            child->U.insert(getChildren().begin(), getChildren().end());
            // Find where to put U list items in child
            for(FNode* u_item : this->U) {
                if(child->is_adjacent(u_item)) {
                    // Adjacent to child
                    child->U.insert(u_item);
                    u_item->U.insert(child);
                } else if(u_item->getDepth() < child->getDepth()) {
                    // Adjacent to parent (this) but not to child
                    child->X.insert(u_item);
                    u_item->W.insert(child);
                } else { // u_item->getDepth() >= child->getDepth()
                    // Find ancestor of u_item that is adjacent to child
                    while(u_item->getDepth() > child->getDepth()) {
                        if(child->is_adjacent(u_item->getParent())) {
                            // Parent is adjacent -> W list
                            child->W.insert(u_item);
                            u_item->X.insert(child);
                            break;
                        } else {
                            u_item = u_item->getParent();
                        }
                    }
                    if(u_item->getDepth() == child->getDepth()) {
                        // No adjacent ancestor -> add a neighbour
                        child->V.insert(u_item);
                        u_item->V.insert(child);
                    }
                }
            }

            // Find where to put W list items in child
            for(FNode* w_item : this->W) {
                // Find first ancestor of w_item that is adjacent to child
                // not needed, done in U list treatment, only check parent
                if(child->getDepth() < w_item->getDepth()) {
                    if(child->is_adjacent(w_item->getParent())) {
                        child->W.insert(w_item);
                        w_item->X.insert(child);
                    }
                } else if(child->getDepth() == w_item->getDepth()) {
                    // No adjacent ancestor -> add a neighbour
                    child->V.insert(w_item);
                    w_item->V.insert(child);
                }
            }
        }
        // Remove this from other lists
        for(FNode* w_item : this->W) {
            w_item->X.erase(this);
        }
        for(FNode* u_item : this->U) {
            u_item->U.erase(this);
        }
        // Clear leaf-only lists
        this->U.clear();
        this->W.clear();

        move_particles_to_children();
    }

    void move_particles_to_children() {
        for(auto&& p : *getParticleContainer()) {
            insert(p);
        }
        delete_particle_container();
    }


    /** Fuses the children nodes until this node becomes a leaf */
    void fuse() {
        if(is_leaf()) {
            return; // In a leaf, there's nothing to do
        }

        for(FNode* child : getChildren()) {
            child->fuse(); // Fuse children into leaves
        }

        create_particle_container();
        _is_leaf = true;

        // Remove children from U lists
        for(FNode* child_1 : getChildren()) {
            for(FNode* child_2 : getChildren()) {
                child_1->U.erase(child_2);
            }
        }
        // Use the children interaction lists to update this one
        for(FNode* child : getChildren()) {
            // Child U list items get into this U list
            for(FNode* u_item : child->U) {
                // Remove child from u_item U list
                u_item->U.erase(child);
                // Add this and u_item in each other's U lists
                this->U.insert(u_item);
                u_item->U.insert(this);
            }
            child->U.clear();

            // Child X list items get into the this U list
            for(FNode* x_item : child->X) {
                // Remove child from x_item W list
                x_item->W.erase(child);
                // Add this and x_item in each other's U lists
                this->U.insert(x_item);
                x_item->U.insert(this);
            }

            // Child W items get into this W list
            for(FNode* w_item : child->W) {
                // Remove child from w_item X list
                w_item->X.erase(child);
                // Add this and w_item in each other's W and X lists
                // when w_item is not adjacent to this
                if(! is_adjacent(w_item)) {
                    this->W.insert(w_item);
                    w_item->X.insert(this);
                }
            }

            // Child V list items get into this W list
            for(FNode* v_item : child->V) {
                // Remove child from the v_item V list
                v_item->V.erase(child);
                // Add this and v_item in each other's W and X lists
                // when v_item is not adjacent to this
                if(! is_adjacent(v_item)) {
                    this->W.insert(v_item);
                    v_item->X.insert(this);
                }
            }

            for(auto&& p : *(child->getParticleContainer())) {
                insert(p);
            }
        }
        // This belongs to its own U list
        this->U.insert(this);
        // Clear the children after everything is updated
        delete_children();

    }

};


#endif
