#include <cassert>
#include <iostream>

#include "test_Abstract.hpp"

#include "Utils/FConstFuncs.hpp"
#include "Adaptive/new/FZCurve.hpp"

namespace scalfmm {
    namespace tests {


        struct test_ZCurve : public test_Abstract<test_ZCurve> {
            using Real = double;
            constexpr static const std::size_t Dim = 3;
            constexpr static const std::size_t pos_count = FZCurve<Dim>::pos_count;
            using boolpos_t = FZCurve<Dim>::position_t;
            using position_t = FPoint<Real, Dim>;

            void run_all() {
                RUN(test_curve_1);
                RUN(test_curve_2);
            }


            // Tests that the position and index functions are the inverse of each other
            void test_curve_1 () {
                FZCurve<Dim> indexer;

                std::array<boolpos_t, pos_count> pos;

                //std::cout << std::endl;
                for(std::size_t i = 0; i < pos_count; ++i) {
                    pos[i] = indexer.position(i);
                    std::size_t j = indexer.index(pos[i]);
                    //std::cout << i << pos[i] << j << std::endl;
                    assert(i == j);
                }
            }

            // Tests that the position and index 2 overload functions are the inverse of each other
            void test_curve_2 () {
                FZCurve<Dim> indexer;
                position_t center = {50,50,50};
                std::array<position_t, pos_count> pos;
                for(std::size_t i = 0; i < pos_count; ++i) {
                    pos[i] = position_t(indexer.position(i)) * 75;
                    std::size_t j = indexer.index(pos[i], center);
                    assert(i == j);
                }
            }

        };

    }
}

int main() {
    scalfmm::tests::test_ZCurve().run_all();
}
