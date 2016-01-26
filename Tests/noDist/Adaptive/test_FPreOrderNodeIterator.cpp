#include <cassert>
#include <iterator>

#include "test_Abstract.hpp"
#include "MockParticle.hpp"

#include "Adaptive/new/FTree.hpp"
#include "Adaptive/new/FNode.hpp"
#include "Adaptive/new/FPrePostOrderNodeIterator.hpp"

#include "Utils/make_unique.hpp"

namespace scalfmm {
    namespace tests {
        struct test_PreOrderNodeIterator : public test_Abstract<test_PreOrderNodeIterator> {

            using Real = double;
            constexpr static std::size_t Dim = 3;

            using tree_t = FTree<std::vector<MockParticle<Real, Dim> >, NodeEmptyData >;
            using node_t = typename tree_t::node_t;
            using box_t  = typename tree_t::box_t;
            using iterator_t = FPreOrderNodeIterator<node_t>;


            tree_t* tree = nullptr;
            node_t* root = nullptr;
            box_t box  {{0,0,0}, 100};


        public:

            /** Sets up a tree in which only the root node and its second child
             *  have children.
             *
             * 0       ___________________o____ ...
             *         |     |     |    |
             * 1       o  ___o____ o    o ...
             *            ||||||||
             * 2          oooooooo
             */
            void set_up() {
                assert(nullptr == tree);
                tree = new tree_t(box);
                root = tree->root();
                root->split();
                root->getChild(1)->split();
            }

            void tear_down() {
                delete tree;
                tree = nullptr;
            }

            void run_all() {
                RUN(test_constructor_default);
                RUN(test_constructor_custom);
                RUN(test_constructor_copy);
                RUN(test_constructor_move);
                RUN(test_operator_copy);
                RUN(test_operator_move);

                RUN(test_prefix_increment);
                RUN(test_prefix_decrement);
                RUN(test_incr_inverse_decr);
                RUN(test_single_node_tree);
            }

            /** Tested for NodeIterator, testing that it is correctly accessible */
            void test_constructor_default() {
                iterator_t a;

                assert(a._root == nullptr);
                assert(a._cur == nullptr);
                assert(a._is_end == false);
                assert(a._is_rend == false);
            }

            void test_constructor_custom() {
                iterator_t a(root);

                assert(a._root == root);
                assert(a._cur  == root);
                assert(a._is_end == false);
                assert(a._is_rend == false);

                assert(a._cur_state == decltype(a._cur_state) ({-1}));
            }


            /** Tested for NodeIterator, testing that it is correctly accessible */
            void test_constructor_copy () {
                iterator_t a(root);
                iterator_t b(a);
            }

            /** Tested for NodeIterator, testing that it is correctly accessible */
            void test_constructor_move () {
                auto a = std::make_unique<iterator_t>(root);
                iterator_t b(std::move(*a));
            }

            /** Tested for NodeIterator, testing that it is correctly accessible */
            void test_operator_copy () {
                iterator_t a(root);
                iterator_t b;
                b = a;
            }

            /** Tested for NodeIterator, testing that it is correctly accessible */
            void test_operator_move () {
                auto a = std::make_unique<iterator_t>(root);
                iterator_t b;
                b = std::move(*a);
            }

            /** Tests a single node tree */
            void test_single_node_tree() {
                tree_t stree(box);
                node_t* sroot = stree.root();

                iterator_t it(sroot);
                assert(it.operator->() == sroot);
                assert(it._is_end == false);
                assert(it._is_rend == false);
                ++it;
                assert(it._is_end == true);
                assert(it._is_rend == false);
                --it;
                assert(it.operator->() == sroot);
                --it;
                assert(it._is_end == false);
                assert(it._is_rend == true);
                ++it;
                assert(it.operator->() == sroot);
                assert(it._is_end == false);
                assert(it._is_rend == false);
            }

            /** Walks the tree, checking each step */
            void test_prefix_increment() {
                iterator_t it(root);

                assert(root == &(*it));
                ++it;
                assert(root->getChild(0) == &(*it));
                ++it;
                assert(root->getChild(1) == &(*it));
                ++it;
                assert(root->getChild(1)->getChild(0) == &(*it));
                ++it;
                assert(root->getChild(1)->getChild(1) == &(*it));
                ++it;
                assert(root->getChild(1)->getChild(2) == &(*it));
                ++it;
                assert(root->getChild(1)->getChild(3) == &(*it));
                ++it;
                assert(root->getChild(1)->getChild(4) == &(*it));
                ++it;
                assert(root->getChild(1)->getChild(5) == &(*it));
                ++it;
                assert(root->getChild(1)->getChild(6) == &(*it));
                ++it;
                assert(root->getChild(1)->getChild(7) == &(*it));
                ++it;
                assert(root->getChild(2) == &(*it));
                ++it;
                assert(root->getChild(3) == &(*it));
                ++it;
                assert(root->getChild(4) == &(*it));
                ++it;
                assert(root->getChild(5) == &(*it));
                ++it;
                assert(root->getChild(6) == &(*it));
                ++it;
                assert(root->getChild(7) == &(*it));
                iterator_t it2(++it);
                ++it;
                assert(it2 == it);

                assert(it2._is_end == true);
                assert(it2._is_rend == false);
            }

            /** Moves iterator past the end and plays back the walk */
            void test_prefix_decrement() {
                iterator_t it(root);

                for(int i = 0; i < 17; ++i) {
                    ++it;
                }
                assert(it._is_end == true);
                --it;
                assert(it._is_end == false);

                assert(root->getChild(7) == &(*it));
                --it;
                assert(root->getChild(6) == &(*it));
                --it;
                assert(root->getChild(5) == &(*it));
                --it;
                assert(root->getChild(4) == &(*it));
                --it;
                assert(root->getChild(3) == &(*it));
                --it;
                assert(root->getChild(2) == &(*it));
                --it;
                assert(root->getChild(1)->getChild(7) == &(*it));
                --it;
                assert(root->getChild(1)->getChild(6) == &(*it));
                --it;
                assert(root->getChild(1)->getChild(5) == &(*it));
                --it;
                assert(root->getChild(1)->getChild(4) == &(*it));
                --it;
                assert(root->getChild(1)->getChild(3) == &(*it));
                --it;
                assert(root->getChild(1)->getChild(2) == &(*it));
                --it;
                assert(root->getChild(1)->getChild(1) == &(*it));
                --it;
                assert(root->getChild(1)->getChild(0) == &(*it));
                --it;
                assert(root->getChild(1) == &(*it));
                --it;
                assert(root->getChild(0) == &(*it));
                --it;
                assert(root == &(*it));
                --it;
                assert(it._is_rend == true);
            }

            void print(iterator_t it) {
                int i = 0;
                while(! it._is_end) {
                    std::cout << i++ << " " << (it++)._cur << std::endl;
                }
            }

            void test_incr_inverse_decr() {
                iterator_t it(root);

                while(! (it)._is_end) {
                    auto cur = it._cur;
                    ++it;
                    auto cur2 = it._cur;
                    --it;
                    assert(cur == it._cur);
                    ++it;
                    assert(cur2 == it._cur);
                }

                while(! (it)._is_rend) {
                    auto cur = it._cur;
                    --it;
                    auto cur2 = it._cur;
                    ++it;
                    assert(cur == it._cur);
                    --it;
                    assert(cur2 == it._cur);
                }

            }


        };
    }
}

int main() {
        scalfmm::tests::test_PreOrderNodeIterator().run_all();
}
