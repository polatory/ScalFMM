#ifndef SCALFMM_TREE_HPP_
#define SCALFMM_TREE_HPP_

#include <unordered_set>

#include "FBox.hpp"
#include "FNode.hpp"
#include "FNodeIteratorBox.hpp"
#include "FInOrderNodeIterator.hpp"
#include "FPrePostOrderNodeIterator.hpp"

template<
    class _ParticleContainer,
    class _NodeData >
class FTree {
public:
    using node_t = FNode<FTree, _ParticleContainer, _NodeData> ;
    using particle_container_t = typename _ParticleContainer::value_type;
    using particle_t = typename _ParticleContainer::value_type;

    using Real = typename particle_t::Real;
    constexpr static std::size_t Dim = particle_t::Dim;

    using position_t = typename particle_t::position_t;
    using box_t = FBox<position_t>;

    using leaf_list_t = std::unordered_set<node_t*>;

    using in_order_walk_t   = FNodeIteratorBox<node_t, FInOrderNodeIterator>;
    using pre_order_walk_t  = FNodeIteratorBox<node_t, FPreOrderNodeIterator>;
    using post_order_walk_t = FNodeIteratorBox<node_t, FPostOrderNodeIterator>;

private:
    std::size_t _max_height = 50;
    std::size_t _leaf_max_particle_count = 50;

    box_t _box;
    node_t* _root = nullptr;
    leaf_list_t _leaves;


public:
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

    FTree(box_t box_) : _box(box_) {
        _root = new node_t(this);
    }

    FTree(FTree&) = delete;

    FTree(FTree&& tree) {
        swap(tree);
    }

    FTree& operator=(FTree&) = delete;

    FTree& operator=(FTree&& tree) {
        swap(tree);
    }

    ~FTree() {
        delete _root;
    }

    std::size_t max_height() const {
        return _max_height;
    }

    void max_height(const std::size_t& height) {
        _max_height = height;
    }

    std::size_t leaf_max_particle_count() const {
        return _leaf_max_particle_count;
    }
    void leaf_max_particle_count(std::size_t count) {
        _leaf_max_particle_count = count;
    }

    node_t* root() {
        return _root;
    }

    const node_t* root() const {
        return _root;
    }

    const box_t& box() const {
        return _box;
    }

    leaf_list_t& leaves() {
        return _leaves;
    }

    const leaf_list_t& leaves() const {
        return _leaves;
    }

    in_order_walk_t in_order_walk() {
        return in_order_walk_t(root());
    }

    pre_order_walk_t pre_order_walk() {
        return pre_order_walk_t(root());
    }

    post_order_walk_t post_order_walk() {
        return post_order_walk_t(root());
    }


    void insert(particle_t p) {
        root()->insert(p);
    }

    void extract(const particle_t& p) {
        root()->extract(p);
    }

    void reshape() {
        root()->reshape();
    }

    bool operator==(const FTree& other) {
        return this == &other;
    }

    bool operator!=(const FTree& other) {
        return ! (*this == other);
    }
};

#endif
