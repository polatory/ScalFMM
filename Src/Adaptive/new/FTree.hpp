#ifndef SCALFMM_TREE_HPP_
#define SCALFMM_TREE_HPP_

#include <unordered_set>

#include "FBox.hpp"
#include "FNode.hpp"
#include "FNodeIteratorBox.hpp"
#include "FInOrderNodeIterator.hpp"
#include "FPrePostOrderNodeIterator.hpp"

/** \brief Adaptive FMM tree
 *
 * \author Quentin Khan <quentin.khan@inria.fr>
 *
 * Implements the adaptive tree to be used with an adaptive FMM algorithm.
 *
 * \tparam _ParticleContainer The class use to store particles. Expected to
 * expose a value_type type definition that exposes the particle type.
 * \tparam _NodeData The cell data class that is used by the FMM kernel.
 */
template<
    class _ParticleContainer,
    class _NodeData >
class FTree {
public:
    /// Internal node structure
    using node_t = FNode<FTree, _ParticleContainer, _NodeData> ;
    /// Particle container type
    using particle_container_t = typename _ParticleContainer::value_type;
    /// Particle type
    using particle_t = typename _ParticleContainer::value_type;
    /// Floating point numbers type
    using Real = typename particle_t::Real;
    /// Space dimension count
    constexpr static std::size_t Dim = particle_t::Dim;
    /// Particle position type
    using position_t = typename particle_t::position_t;
    /// Box type use to slice space
    using box_t = FBox<position_t>;
    /// List type to store leaves
    using leaf_list_t = std::unordered_set<node_t*>;
    /// In order tree traversal mock container type
    using in_order_walk_t   = FNodeIteratorBox<node_t, FInOrderNodeIterator>;
    /// Pre-order tree traversal mock container type
    using pre_order_walk_t  = FNodeIteratorBox<node_t, FPreOrderNodeIterator>;
    /// Post-order tree traversal mock container type
    using post_order_walk_t = FNodeIteratorBox<node_t, FPostOrderNodeIterator>;

private:
    /// Tree maximum heigh
    std::size_t _max_height = 50;
    /// Maximum particle per leaf density
    std::size_t _leaf_max_particle_count = 50;
    /// Tree space bounding box
    box_t _box;
    /// Tree root node
    node_t* _root = nullptr;
    /// Tree leaf list
    leaf_list_t _leaves;


public:
    /** \brief Swaps two trees */
    void swap(FTree& tree) {
        using std::swap;
        swap(_max_height, tree._max_height);
        swap(_box, tree._box);
        swap(_root, tree._root);
        swap(_leaves, tree._leaves);
        swap(_leaf_max_particle_count, tree._leaf_max_particle_count);

        if(tree.root() && root()) {
            tree.root()->setTree(&(root()->getTree()));
        }
        if(root()) {
            root()->setTree(this);
        }
    }

    /** \brief Builds a tree
     *
     * \param box Tree bounding box
     *
     */
    FTree(box_t box_) : _box(box_) {
        _root = new node_t(this);
    }

    FTree(FTree&) = delete;

    /** \brief Move constructor */
    FTree(FTree&& tree) {
        swap(tree);
    }

    FTree& operator=(FTree&) = delete;

    /** \brief Move assignment */
    FTree& operator=(FTree&& tree) {
        swap(tree);
    }

    /** \brief Destructor */
    ~FTree() {
        delete _root;
    }

    /** \brief Maximum height accessor
     * \return Tree maximum height
     */
    std::size_t max_height() const {
        return _max_height;
    }

    /** \brief Maximum height setter
     * \param height New maximum height
     */
    void max_height(const std::size_t& height) {
        _max_height = height;
    }

    /** \brief Maximum leaf particle density accessor
     * \return Tree maximum leaf particle density
     */
    std::size_t leaf_max_particle_count() const {
        return _leaf_max_particle_count;
    }
    /** \brief Maximum leaf particle density accessor
     * \param count New maximum density
     */
    void leaf_max_particle_count(std::size_t count) {
        _leaf_max_particle_count = count;
    }

    /** \brief Tree root accessor
     * \return Root pointer
     */
    node_t* root() {
        return _root;
    }

    /** \brief Tree root const accessor
     * \return Root pointer
     */
    const node_t* root() const {
        return _root;
    }

    /** \brief Bounding box accessor
     * \return Tree bounding box
     */
    const box_t& box() const {
        return _box;
    }

    /** \brief Leaf list accessor
    * \return Tree leaf list
    */
    leaf_list_t& leaves() {
        return _leaves;
    }

    /** \brief Leaf list const accessor
    * \return Tree leaf list
    */
    const leaf_list_t& leaves() const {
        return _leaves;
    }

    /** \brief In order walk accessor */
    in_order_walk_t in_order_walk() {
        return in_order_walk_t(root());
    }

    /** \brief Pre-order walk accessor */
    pre_order_walk_t pre_order_walk() {
        return pre_order_walk_t(root());
    }

    /** \brief Post-order walk accessor */
    post_order_walk_t post_order_walk() {
        return post_order_walk_t(root());
    }

    /** \brief Proxy call for FNode#insert applied to #root*/
    void insert(const particle_t& p) {
        root()->insert(p);
    }


    /** \brief Proxy call for FNode#extract appliced to #root*/
    void extract(const particle_t& p) {
        root()->extract(p);
    }

    /** \brief Updates the underlying graph after moving some particles
     *
     * After an FMM run, particles may have moved out of their leaf's bounding
     * box. This extracts the particles from the leaf and inserts them again at
     * the right place.
     *
     */
    void reshape() {
        root()->reshape();
    }

};

#endif
