#include <cassert>

#include "test_Abstract.hpp"
#include "MockParticle.hpp"

#include "Adaptive/new/FTree.hpp"
#include "Adaptive/new/FNode.hpp"
#include "Adaptive/new/FNodeIterator.hpp"

#include "Utils/make_unique.hpp"

namespace scalfmm {
    namespace tests {

        struct test_NodeIterator : public test_Abstract<test_NodeIterator> {

            using Real = double;
            constexpr static std::size_t Dim = 3;

            using tree_t = FTree<std::vector<MockParticle<Real, Dim> >, NodeEmptyData >;
            using node_t = typename tree_t::node_t;
            using box_t  = typename tree_t::box_t;

            struct test_iterator : public FNodeIterator<test_iterator, node_t> {
                using base_t = FNodeIterator<test_iterator, node_t>;

                int pp = 0;
                int mm = 0;

                using base_t::base_t;
                using base_t::operator++;
                using base_t::operator--;

                template<base_t::Direction> void move_to_boundary() {}
                test_iterator& operator++(){++pp; return *this;}
                test_iterator& operator--(){++mm; return *this;}
            };

            tree_t* tree = nullptr;
            node_t* root = nullptr;
            box_t box  {{0,0,0}, 100};


        public:

            void set_up() {
                assert(nullptr == tree);
                tree = new tree_t(box);
                root = tree->root(); // Allocated by tree constructor
            }

            void tear_down() {
                delete tree;
                tree = nullptr;
                root = nullptr; // Deleted by tree destructor
            }

            void run_all() {
                RUN(test_constructor_default);
                RUN(test_constructor_custom);
                RUN(test_constructor_copy);
                RUN(test_constructor_move);
                RUN(test_operator_copy);
                RUN(test_operator_move);
                RUN(test_operator_dereference);
                RUN(test_increment_postfix);
                RUN(test_decrement_postfix);
            }

            void test_constructor_default() {
                test_iterator a;
                assert(a._root == nullptr);
                assert(a._cur == nullptr);
                assert(a._is_end == false);
                assert(a._is_rend == false);
            }

            void test_constructor_custom() {
                test_iterator a(root);

                assert(a._root == root);
                assert(a._is_end == false);
                assert(a._is_rend == false);
            }

            void test_constructor_copy() {
                test_iterator a(root);
                a._is_end = true;
                a._is_rend = true;
                test_iterator b(a);

                assert(a._root      == b._root     );
                assert(a._cur       == b._cur      );
                assert(a._is_end    == b._is_end   );
                assert(a._is_rend   == b._is_rend  );
                assert(a._cur_state == b._cur_state);
            }

            void test_constructor_move() {
                auto a = std::make_unique<test_iterator>(root);
                a->_is_end = true;
                a->_is_rend = true;

                auto _root = a->_root;
                auto _cur  = a->_cur;
                auto _cur_state  = a->_cur_state;
                auto _is_end = a->_is_end;
                auto _is_rend = a->_is_rend;

                test_iterator b(std::move(*a));

                assert(_root      == b._root     );
                assert(_cur       == b._cur      );
                assert(_is_end    == b._is_end   );
                assert(_is_rend   == b._is_rend  );
                assert(_cur_state == b._cur_state);
            }

            void test_operator_copy() {
                test_iterator a(root);
                a._is_end = true;
                a._is_rend = true;
                test_iterator b;
                b = a;

                assert(a._root      == b._root     );
                assert(a._cur       == b._cur      );
                assert(a._is_end    == b._is_end   );
                assert(a._is_rend   == b._is_rend  );
                assert(a._cur_state == b._cur_state);
            }

            void test_operator_move() {
                auto a = std::make_unique<test_iterator>(root);
                a->_is_end = true;
                a->_is_rend = true;

                auto _root = a->_root;
                auto _cur  = a->_cur;
                auto _cur_state  = a->_cur_state;
                auto _is_end = a->_is_end;
                auto _is_rend = a->_is_rend;

                test_iterator b;
                b = std::move(*a);

                assert(_root      == b._root     );
                assert(_cur       == b._cur      );
                assert(_is_end    == b._is_end   );
                assert(_is_rend   == b._is_rend  );
                assert(_cur_state == b._cur_state);
            }

            void test_operator_dereference() {
                test_iterator a(root);
                a._cur = root;
                assert(&(*a) == root);

                const test_iterator b(a);
                assert(&(*b) == root);

                auto c = std::make_unique<test_iterator>(a);
                assert(&((*c)->getTree()) == tree);

                auto d = std::make_unique<const test_iterator>(a);
                assert(&((*d)->getTree()) == tree);
            }

            void test_increment_postfix() {
                test_iterator a(root);
                int plus = (a++).pp;
                assert(plus == 0);
                assert(a.pp == 1);
            }

            void test_decrement_postfix() {
                test_iterator a(root);
                int minus = (a--).mm;
                assert(minus == 0);
                assert(a.mm == 1);
            }


        };
    }
}


int main() {
    scalfmm::tests::test_NodeIterator().run_all();
}
