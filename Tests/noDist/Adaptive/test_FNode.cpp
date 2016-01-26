#include <cassert>
#include <random>

#include <unordered_set>

#include "test_Abstract.hpp"
#include "MockParticle.hpp"

#include "Adaptive/new/FNode.hpp"


namespace scalfmm {
    namespace tests {

        struct test_Node : test_Abstract<test_Node> {

            template<typename _Real, std::size_t _Dim>
            struct DummyTree {
                constexpr static const std::size_t Dim = _Dim;
                using Real = _Real;
                using particle_t = MockParticle<Real, Dim>;
                using node_t = FNode<DummyTree, std::vector<particle_t >, NodeEmptyData >;
                using leaf_list_t = std::unordered_set<node_t*>;
                using position_t = typename node_t::position_t;
                using box_t = typename node_t::box_t;

                box_t _box;
                node_t* _root;
                leaf_list_t _leaves;
                std::size_t _leaf_max_particle_count = 10;

                DummyTree(Real side_length) {
                    _box.set(_box.c1, position_t()+side_length);
                    _root = new node_t(this);
                }

                ~DummyTree() {
                    delete _root;
                }

                box_t box() const {return _box;}
                node_t* root() const {return _root;}
                leaf_list_t& leaves() {return _leaves;}
                std::size_t max_height() const {return 20;}
                std::size_t leaf_max_particle_count() const {return _leaf_max_particle_count;}
                friend bool operator==(const DummyTree& t1, const DummyTree& t2) {
                    return &t1 == &t2;
                }
            };



            using Real = double;
            static constexpr std::size_t Dim = 3;
            using tree_t = DummyTree<Real, Dim>;
            using node_t = typename tree_t::node_t;
            using box_t = typename tree_t::box_t;
            using position_t = typename tree_t::position_t;
            using particle_t = typename tree_t::particle_t;


            static constexpr Real box_side_length = 100;
            static constexpr std::size_t particle_count = 1000;

            tree_t* tree;

            std::random_device rd;
            std::mt19937 random_generator;
            std::uniform_int_distribution<std::size_t> random_index;
            std::uniform_real_distribution<Real> random_real;

            test_Node() :
                random_generator(rd()),
                random_index(0, Fpow(2,Dim)-1),
                random_real(0, box_side_length)
                {}

            void set_up() {
                tree = new tree_t(box_side_length);
            }

            void tear_down() {
                delete tree;
            }

            void run_all() {
                RUN(test_constructor_from_tree);
                RUN(test_constructor_from_parent);
                RUN(test_getters);
                RUN(test_is_adjacent);
                RUN(test_in_order_traversal);
                RUN(test_post_order_traversal);
                RUN(test_pre_order_traversal);
                RUN(test_split_simple);
                RUN(test_split_max_height);
                RUN(test_split_interaction_lists_update);
                RUN(test_fuse_interaction_lists_update);
                RUN(test_insert);
                RUN(test_insert_random);
                RUN(test_reshape_fuse);
                RUN(test_reshape_split);
            }


            void test_constructor_from_tree() {
                assert(tree->root()->getDepth()  == 0);
                assert(tree->root()->getParent() == nullptr);
                assert(tree->root()->getBox()    == tree->box());
                assert(tree->root()->getIndex()  == 0);
                for(auto child : tree->root()->getChildren()) {
                    assert(child == nullptr);
                }

                assert(tree->leaves().size() == 1);
                assert(tree->leaves().count(tree->root()) == 1);

            }

            void test_constructor_from_parent() {
                try {
                    node_t* node = new node_t(*(tree->root()), Fpow(2,Dim)+1);
                    assert(false);
                    node->getParent();
                } catch(std::invalid_argument& e) { }


                const node_t* node = new node_t(*(tree->root()), 5);

                box_t b = box_t( tree->box().center(), tree->box().corner(5) );

                assert(node->getDepth() == 1);
                assert(node->getParent() == tree->root());
                assert(node->getBox() == b);
                assert(node->getIndex() == 5);
                for(auto child : node->getChildren()) {
                    assert(child == nullptr);
                }

                delete node;
            }

            void test_getters() {
                tree->root()->create_children();
                tree->root()->delete_particle_container();
                node_t* node = tree->root();
                node_t* leaf = tree->root()->getChild(0);

                // Test getters on a leaf
                //assert(leaf->data() == NodeEmptyData{});
                assert(leaf->getChildren() == typename node_t::child_node_array_t{{}});
                assert(leaf->getChild(0) == nullptr);
                assert(leaf->getParent() == node);
                assert(leaf->getDepth() == 1);
                assert(leaf->getIndex() == 0);
                assert(leaf->getTree() == *tree);
                assert(leaf->getBox() == box_t(node->getBox().center(), node->getBox().corner(0)));
                assert(leaf->getParticleContainer() != nullptr);
                leaf->getParticleContainer()->push_back({{0,0,0}});
                assert(leaf->getParticleCount() == 1);


                // Test getters on an inner node
                //assert(leaf->data() == NodeEmptyData());
                assert(node->getChild(0) == leaf);
                assert(node->getChild(0) == leaf);
                assert(node->getParent() == nullptr);
                assert(node->getDepth() == 0);
                assert(node->getIndex() == 0);
                assert(node->getTree() == *tree);
                assert(node->getBox() == tree->box() );
                assert(node->getParticleContainer() == nullptr);
                assert(node->getParticleCount() == 0);

                assert(*(node->getChild(0)) == *leaf);

            }

            void test_is_adjacent() {
                // A node is adjacent to itself
                node_t* root = tree->root();
                assert( root->is_adjacent(root) );
                assert(! root->is_adjacent(nullptr) );

                tree->root()->create_children();
                // A node is not adjacent to its children
                assert( ! root->is_adjacent(root->getChild(0)) );
                assert( ! root->is_adjacent(root->getChild(1)) );
                assert( ! root->is_adjacent(root->getChild(2)) );
                assert( ! root->is_adjacent(root->getChild(3)) );
                assert( ! root->is_adjacent(root->getChild(4)) );
                assert( ! root->is_adjacent(root->getChild(5)) );
                assert( ! root->is_adjacent(root->getChild(6)) );
                assert( ! root->is_adjacent(root->getChild(7)) );

                // Sibling nodes are adjacent
                auto check_siblings = [](node_t* node, node_t* parent) {
                    for(std::size_t i = 0; i < node_t::child_count; ++i) {
                        assert( node->is_adjacent(parent->getChild(0)) );
                    }
                };
                for(std::size_t i = 0; i < node_t::child_count; ++i) {
                    node_t* child = tree->root()->getChild(i);
                    check_siblings(child, tree->root());

                    child->create_children();
                    for(std::size_t j = 0; j < node_t::child_count; ++j) {
                        assert( ! child->is_adjacent(child->getChild(j)) );
                        check_siblings(child->getChild(j), child);
                    }
                }

                // Some random checks
                auto get_node = [&root](std::vector<std::size_t> indexes) {
                    node_t* node = root;
                    for(std::size_t p : indexes)
                        node = node->getChild(p);
                    return node;
                };

                assert(get_node({0,7})->is_adjacent(get_node({1,0})));
                //TODO: add tests btw nodes that have different parents
            }


            void test_in_order_traversal() {
                std::vector<std::pair<std::size_t,std::size_t>> indexes;
                std::function<void(node_t*)> lambda = [&](node_t* node) {
                    indexes.emplace_back(node->getDepth(), node->getIndex());
                };

                tree->root()->create_children();
                tree->root()->getChild(1)->create_children();

                tree->root()->for_each_in_order(lambda);

                assert((indexes[ 0] == std::pair<std::size_t, std::size_t>(1, 0)) );
                assert((indexes[ 1] == std::pair<std::size_t, std::size_t>(2, 8)) );
                assert((indexes[ 2] == std::pair<std::size_t, std::size_t>(2, 9)) );
                assert((indexes[ 3] == std::pair<std::size_t, std::size_t>(2,10)) );
                assert((indexes[ 4] == std::pair<std::size_t, std::size_t>(2,11)) );
                assert((indexes[ 5] == std::pair<std::size_t, std::size_t>(1, 1)) );
                assert((indexes[ 6] == std::pair<std::size_t, std::size_t>(2,12)) );
                assert((indexes[ 7] == std::pair<std::size_t, std::size_t>(2,13)) );
                assert((indexes[ 8] == std::pair<std::size_t, std::size_t>(2,14)) );
                assert((indexes[ 9] == std::pair<std::size_t, std::size_t>(2,15)) );
                assert((indexes[10] == std::pair<std::size_t, std::size_t>(1, 2)) );
                assert((indexes[11] == std::pair<std::size_t, std::size_t>(1, 3)) );
                assert((indexes[12] == std::pair<std::size_t, std::size_t>(0, 0)) );
                assert((indexes[13] == std::pair<std::size_t, std::size_t>(1, 4)) );
                assert((indexes[14] == std::pair<std::size_t, std::size_t>(1, 5)) );
                assert((indexes[15] == std::pair<std::size_t, std::size_t>(1, 6)) );
                assert((indexes[16] == std::pair<std::size_t, std::size_t>(1, 7)) );
            }

            void test_post_order_traversal() {
                std::vector<std::pair<std::size_t,std::size_t>> indexes;
                std::function<void(node_t*)> lambda = [&](node_t* node) {
                    indexes.emplace_back(node->getDepth(), node->getIndex());
                };

                tree->root()->create_children();
                tree->root()->getChild(1)->create_children();

                tree->root()->for_each_post_order(lambda);

                assert((indexes[ 0] == std::pair<std::size_t, std::size_t>(1, 0)) );
                assert((indexes[ 1] == std::pair<std::size_t, std::size_t>(2, 8)) );
                assert((indexes[ 2] == std::pair<std::size_t, std::size_t>(2, 9)) );
                assert((indexes[ 3] == std::pair<std::size_t, std::size_t>(2,10)) );
                assert((indexes[ 4] == std::pair<std::size_t, std::size_t>(2,11)) );
                assert((indexes[ 5] == std::pair<std::size_t, std::size_t>(2,12)) );
                assert((indexes[ 6] == std::pair<std::size_t, std::size_t>(2,13)) );
                assert((indexes[ 7] == std::pair<std::size_t, std::size_t>(2,14)) );
                assert((indexes[ 8] == std::pair<std::size_t, std::size_t>(2,15)) );
                assert((indexes[ 9] == std::pair<std::size_t, std::size_t>(1, 1)) );
                assert((indexes[10] == std::pair<std::size_t, std::size_t>(1, 2)) );
                assert((indexes[11] == std::pair<std::size_t, std::size_t>(1, 3)) );
                assert((indexes[12] == std::pair<std::size_t, std::size_t>(1, 4)) );
                assert((indexes[13] == std::pair<std::size_t, std::size_t>(1, 5)) );
                assert((indexes[14] == std::pair<std::size_t, std::size_t>(1, 6)) );
                assert((indexes[15] == std::pair<std::size_t, std::size_t>(1, 7)) );
                assert((indexes[16] == std::pair<std::size_t, std::size_t>(0, 0)) );

            }

            void test_pre_order_traversal() {
                std::vector<std::pair<std::size_t,std::size_t>> indexes;
                std::function<void(node_t*)> lambda = [&](node_t* node) {
                    indexes.emplace_back(node->getDepth(), node->getIndex());
                };

                tree->root()->create_children();
                tree->root()->getChild(1)->create_children();

                tree->root()->for_each_pre_order(lambda);

                assert((indexes[ 0] == std::pair<std::size_t, std::size_t>(0, 0)) );
                assert((indexes[ 1] == std::pair<std::size_t, std::size_t>(1, 0)) );
                assert((indexes[ 2] == std::pair<std::size_t, std::size_t>(1, 1)) );
                assert((indexes[ 3] == std::pair<std::size_t, std::size_t>(2, 8)) );
                assert((indexes[ 4] == std::pair<std::size_t, std::size_t>(2, 9)) );
                assert((indexes[ 5] == std::pair<std::size_t, std::size_t>(2,10)) );
                assert((indexes[ 6] == std::pair<std::size_t, std::size_t>(2,11)) );
                assert((indexes[ 7] == std::pair<std::size_t, std::size_t>(2,12)) );
                assert((indexes[ 8] == std::pair<std::size_t, std::size_t>(2,13)) );
                assert((indexes[ 9] == std::pair<std::size_t, std::size_t>(2,14)) );
                assert((indexes[10] == std::pair<std::size_t, std::size_t>(2,15)) );
                assert((indexes[11] == std::pair<std::size_t, std::size_t>(1, 2)) );
                assert((indexes[12] == std::pair<std::size_t, std::size_t>(1, 3)) );
                assert((indexes[13] == std::pair<std::size_t, std::size_t>(1, 4)) );
                assert((indexes[14] == std::pair<std::size_t, std::size_t>(1, 5)) );
                assert((indexes[15] == std::pair<std::size_t, std::size_t>(1, 6)) );
                assert((indexes[16] == std::pair<std::size_t, std::size_t>(1, 7)) );

            }


            struct interaction_list_struct {
                using interaction_list_t = typename node_t::interaction_list_t;
                interaction_list_t U;
                interaction_list_t V;
                interaction_list_t W;
                interaction_list_t X;

                interaction_list_struct(tree_t* tree, node_t* node) {
                    tree->root()->for_each_in_order([&](node_t* n) {
                            if(n == node && n->is_leaf()) {
                                this->U.insert(n);
                                return;
                            }

                            if(node->is_leaf()) {
                                if(n->is_leaf() && n->is_adjacent(node)) {
                                    this->U.insert(n);
                                }
                                if(n->getParent() && node->is_adjacent(n->getParent())
                                   && ! node->is_adjacent(n)
                                   && node->getDepth() < n->getDepth() ) {
                                    this->W.insert(n);
                                }
                            }

                            if(n->getDepth() == node->getDepth()) {
                                if(node->getParent()
                                   && n->getParent()
                                   && ! node->is_adjacent(n)
                                   && n->getParent() != node->getParent()
                                   && node->getParent()->is_adjacent(n->getParent())) {
                                    this->V.insert(n);
                                }
                            }

                            if(n->is_leaf()) {
                                if(node->getParent()
                                   && n->is_adjacent(node->getParent())
                                   && ! n->is_adjacent(node)
                                   && n->getDepth() < node->getDepth()) {
                                    this->X.insert(n);
                                }
                            }
                        });
                }
            };

            void test_split_simple() {
                // Split the root
                tree->root()->split();
                // The root is no longer a leaf
                assert(! tree->root()->is_leaf());
                // The root's children have been created
                assert(tree->root()->getChild(0) != nullptr);
                assert(tree->root()->getChild(1) != nullptr);
                assert(tree->root()->getChild(2) != nullptr);
                assert(tree->root()->getChild(3) != nullptr);
                assert(tree->root()->getChild(4) != nullptr);
                assert(tree->root()->getChild(5) != nullptr);
                assert(tree->root()->getChild(6) != nullptr);
                assert(tree->root()->getChild(7) != nullptr);
                // Tree leaf list update : each child must be in the list once
                assert(tree->leaves().count(tree->root()) == 0);
                assert(tree->leaves().count(tree->root()->getChild(0)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(1)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(2)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(3)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(4)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(5)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(6)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(7)) == 1);
                // Split a child
                node_t* node = tree->root()->getChild(1);
                node->split();
                // Node is no longer a leaf
                assert(! node->is_leaf());
                // Check that the node's children have been added to the leaf list
                assert(tree->leaves().count(node) == 0);
                assert(tree->leaves().count(node->getChild(0)) == 1);
                assert(tree->leaves().count(node->getChild(1)) == 1);
                assert(tree->leaves().count(node->getChild(2)) == 1);
                assert(tree->leaves().count(node->getChild(3)) == 1);
                assert(tree->leaves().count(node->getChild(4)) == 1);
                assert(tree->leaves().count(node->getChild(5)) == 1);
                assert(tree->leaves().count(node->getChild(6)) == 1);
                assert(tree->leaves().count(node->getChild(7)) == 1);
                // other leaves are still in the list
                assert(tree->leaves().count(tree->root()->getChild(0)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(2)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(3)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(4)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(5)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(6)) == 1);
                assert(tree->leaves().count(tree->root()->getChild(7)) == 1);
            }

            void test_split_max_height() {
                node_t* node = tree->root();
                std::size_t i;
                // Check that the split limit works
                for(i = 0; i < tree->max_height() && node != nullptr; ++i) {
                    node->split();
                    node = node->getChild(1);
                }

                assert(node != nullptr);
                assert(node->getDepth() == tree->max_height());

                node->split();
                for(auto c : node->getChildren()) {
                    assert(c == nullptr);
                }

            }

            std::size_t random_child_index() {
                return random_index(random_generator);
            }

            void setup_random_tree(const std::size_t& split_count) {
                // Randomly split N leaves
                for(std::size_t i = 0; i < split_count; ++i) {
                    node_t* node = tree->root();
                    // Follow random path until we hit a leaf
                    while( ! node->is_leaf()) {
                        node = node->getChild(random_child_index());
                    }
                    node->split();
                }
            }



            void test_split_interaction_lists_update() {
                setup_random_tree(20);

                tree->root()->for_each_in_order([&](node_t* node) {
                        interaction_list_struct lists(tree, node);
                        assert(node->U == lists.U);
                        assert(node->V == lists.V);
                        assert(node->W == lists.W);
                        assert(node->X == lists.X);
                    });
            }

            void test_fuse_interaction_lists_update() {
                setup_random_tree(10);

                while(! tree->root()->is_leaf()) {
                    node_t* cursor = tree->root();
                    while(! cursor->is_leaf()) {
                        cursor = cursor->getChild(random_child_index());
                    }

                    cursor->getParent()->fuse();

                    tree->root()->for_each_in_order([&](node_t* node) {
                            interaction_list_struct lists(tree, node);
                            assert(node->U == lists.U);
                            assert(node->V == lists.V);
                            assert(node->W == lists.W);
                            assert(node->X == lists.X);
                        });
                }
            }

            position_t random_position() {
                position_t p;
                for(auto&& i : p)
                    i = random_real(random_generator);
                return p;
            }

            void test_insert() {
                tree->_leaf_max_particle_count = 2;

                node_t* root = tree->root();

                root->insert({{24, 24, 24},1}); // goes in first child
                root->insert({{76, 76, 76},2}); // goes in last child
                root->insert({{26, 30, 74},3}); // goes in second child

                // _leaf_max_particle_count is exceeded, the node should have split

                assert(! root->is_leaf());
                assert(root->getChild(0)->is_leaf());
                assert(root->getChild(1)->is_leaf());
                assert(root->getChild(2)->is_leaf());
                assert(root->getChild(3)->is_leaf());
                assert(root->getChild(4)->is_leaf());
                assert(root->getChild(5)->is_leaf());
                assert(root->getChild(6)->is_leaf());
                assert(root->getChild(7)->is_leaf());

                std::size_t particle_count = 0;
                for(auto child : root->getChildren()) {
                    for(auto&& particle : *(child->getParticleContainer())) {
                        assert(child->getBox().contains(particle));
                        ++particle_count;
                    }
                }

                assert(particle_count == 3);

            }


            void test_insert_random() {
                std::vector<particle_t> particle_list;
                particle_list.reserve(particle_count);

                for(std::size_t i = 0; i < particle_count; ++i) {
                    particle_list.push_back(particle_t{random_position(), i});
                    tree->root()->insert(particle_list.back());
                }

                std::size_t particle_counter = 0;

                tree->root()->for_each_in_order([&](node_t* node) {
                        if(! node->is_leaf()) {
                            assert(node->getParticleContainer() == nullptr);
                        } else {
                            for(auto&& p : *(node->getParticleContainer())) {
                                ++particle_counter;
                                assert(node->getBox().contains(p.position()));
                            }
                        }
                    });

                assert(particle_count == particle_counter);
            }


            void test_reshape_fuse() {
                // Set max particle count
                tree->_leaf_max_particle_count = 8;
                // Insert exactly max particle count particles
                tree->root()->insert({{12.5,12.5,12.5}});
                tree->root()->insert({{12.5,12.5,37.5}});
                tree->root()->insert({{12.5,37.5,12.5}});
                tree->root()->insert({{12.5,37.5,37.5}});
                tree->root()->insert({{37.5,12.5,12.5}});
                tree->root()->insert({{37.5,12.5,37.5}});
                tree->root()->insert({{37.5,37.5,12.5}});
                tree->root()->insert({{37.5,37.5,37.5}});
                // Check that root has not yet been divided
                assert(tree->root()->is_leaf());
                // Add one more particle
                tree->root()->insert({{12.5,12.5,62.5}});
                // Check that root has been divided
                assert(! tree->root()->is_leaf());
                // Add particles until 2*max_particle_count
                tree->root()->insert({{12.5,12.5,87.5}});
                tree->root()->insert({{12.5,37.5,62.5}});
                tree->root()->insert({{12.5,37.5,87.5}});
                tree->root()->insert({{37.5,12.5,62.5}});
                tree->root()->insert({{37.5,12.5,87.5}});
                tree->root()->insert({{37.5,37.5,62.5}});
                tree->root()->insert({{37.5,37.5,87.5}});

                // Remove max_particle_count particles
                std::size_t i = 0;
                for(node_t* leaf : tree->leaves()) {
                    while(leaf->getParticleCount() > tree->_leaf_max_particle_count / 2) {
                        leaf->extract(0);
                        if(++i >= tree->_leaf_max_particle_count) {
                            break;
                        }
                    }
                }

                tree->root()->reshape();

                assert(tree->root()->is_leaf());
            }

            void test_reshape_split() {
                // Set max particle count
                tree->_leaf_max_particle_count = 8;
                // Insert exactly max particle count particles
                tree->root()->getParticleContainer()->push_back({{12.5,12.5,12.5}});
                tree->root()->getParticleContainer()->push_back({{12.5,12.5,37.5}});
                tree->root()->getParticleContainer()->push_back({{12.5,37.5,12.5}});
                tree->root()->getParticleContainer()->push_back({{12.5,37.5,37.5}});
                tree->root()->getParticleContainer()->push_back({{37.5,12.5,12.5}});
                tree->root()->getParticleContainer()->push_back({{37.5,12.5,37.5}});
                tree->root()->getParticleContainer()->push_back({{37.5,37.5,12.5}});
                tree->root()->getParticleContainer()->push_back({{37.5,37.5,37.5}});

                assert(tree->root()->is_leaf());
                tree->root()->reshape();
                assert(tree->root()->is_leaf());

                // Insert max particle count particles again
                tree->root()->getParticleContainer()->push_back({{12.5,12.5,62.5}});
                tree->root()->getParticleContainer()->push_back({{12.5,12.5,87.5}});
                tree->root()->getParticleContainer()->push_back({{12.5,37.5,62.5}});
                tree->root()->getParticleContainer()->push_back({{12.5,37.5,87.5}});
                tree->root()->getParticleContainer()->push_back({{37.5,12.5,62.5}});
                tree->root()->getParticleContainer()->push_back({{37.5,12.5,87.5}});
                tree->root()->getParticleContainer()->push_back({{37.5,37.5,62.5}});
                tree->root()->getParticleContainer()->push_back({{37.5,37.5,87.5}});
                // Check that root has not yet been divided
                assert(tree->root()->is_leaf());
                // Reshape tree
                tree->root()->reshape();
                // Root shouldn't be a leaf anymore
                assert(! tree->root()->is_leaf());
            }

        };
    }
}


int main() {
    scalfmm::tests::test_Node().run_all();
}
