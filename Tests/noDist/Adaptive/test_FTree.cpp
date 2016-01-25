#include <cassert>

#include "test_Abstract.hpp"
#include "MockParticle.hpp"

#include "Adaptive/new/FNode.hpp"
#include "Adaptive/new/FTree.hpp"

#define BOX_SIDE_LEN 100


namespace scalfmm {
    namespace tests {

        class test_Tree : public test_Abstract<test_Tree> {

            using Real = double;
            constexpr static std::size_t Dim = 3;
            using particle_t           = MockParticle<Real, Dim>;
            using particle_container_t = std::vector<particle_t>;
            using tree_t               = FTree<particle_container_t, NodeEmptyData>;
            using node_t               = typename tree_t::node_t;
            using box_t                = typename tree_t::box_t;
            using position_t           = typename box_t::position_t;

            box_t box {position_t(), 100};
            tree_t* tree = nullptr;

        public:
            void run_all() {
                RUN(test_constructor);
                RUN(test_swap);
                RUN(test_constructor_move);
            }

            void set_up() {
                assert(tree == nullptr);

                this->tree = new tree_t(box);
            }

            void tear_down() {
                delete this->tree;
                this->tree = nullptr;
            }

            void test_constructor() {
                assert(tree->root() != nullptr);
            }

            void test_swap() {
                box_t t_box = {position_t(), 50};
                tree_t t(t_box);

                auto root = tree->root();
                auto t_root = t.root();

                auto leaves = tree->leaves();
                auto t_leaves = t.leaves();

                auto max_height = tree->max_height();
                auto t_max_height = t.max_height();

                auto leaf_max_particle_count = tree->leaf_max_particle_count();
                auto t_leaf_max_particle_count = t.leaf_max_particle_count();

                t.swap(*tree);

                assert(tree->box()        == t_box);
                assert(tree->root()       == t_root);
                assert(tree->leaves()     == t_leaves);
                assert(tree->max_height() == t_max_height);
                assert(tree->leaf_max_particle_count() == t_leaf_max_particle_count);
                assert(tree->root()->getTree() == *tree);

                assert(t.box()        == box);
                assert(t.root()       == root);
                assert(t.leaves()     == leaves);
                assert(t.max_height() == max_height);
                assert(t.leaf_max_particle_count() == leaf_max_particle_count);
                assert(t.root()->getTree() == t);
            }

            void test_constructor_move() {
                auto root = tree->root();
                auto leaves = tree->leaves();
                auto max_height = tree->max_height();
                auto leaf_max_particle_count = tree->leaf_max_particle_count();

                tree_t t(std::move(*tree));

                assert(t.box()        == box       );
                assert(t.root()       == root      );
                assert(t.leaves()     == leaves    );
                assert(t.max_height() == max_height);
                assert(t.leaf_max_particle_count() == leaf_max_particle_count);

                tree_t s {tree_t(box)};
            }

        };
    }
}


int main() {
    scalfmm::tests::test_Tree().run_all();
}
