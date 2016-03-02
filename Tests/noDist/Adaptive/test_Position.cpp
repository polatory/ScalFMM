#include <cassert>
#include <iostream>

#include "test_Abstract.hpp"

#include "Utils/FPoint.hpp"
#include "Utils/FConstFuncs.hpp"

namespace scalfmm {
    namespace tests {

        using Real = double;
        constexpr const size_t Dim = 3;
        using position = FPoint<Real, Dim>;


        struct test_Position : public test_Abstract<test_Position> {

            void run_all() {
                RUN(test_constructor_default   );
                RUN(test_constructor_arguments   );
                RUN(test_constructor_copy   );
                RUN(test_constructor_move   );
                RUN(test_operator_assignment_copy   );
                RUN(test_operator_assignment_move   );
                RUN(test_operator_equal_int);
                RUN(test_operator_equal_float);
                RUN(test_operator_plus );
                RUN(test_operator_mult );
                RUN(test_iterators     );
            }

            void test_constructor_default() {
                position a;
                for(Real& i : a) {
                    assert(Ffeq(i, 0));
                }
            }

            void test_constructor_arguments() {
                position b(0.1, 0.2, 0.3);
                assert(Ffeq(b[0], 0.1));
                assert(Ffeq(b[1], 0.2));
                assert(Ffeq(b[2], 0.3));
            }

            void test_constructor_copy() {
                position b(0.1, 0.2, 0.3);
                position c(b);
                assert(Ffeq(b[0], c[0]));
                assert(Ffeq(b[1], c[1]));
                assert(Ffeq(b[2], c[2]));

                c[0] = b[0]+1;
                assert(! Ffeq(b[0], c[0]));
            }

            void test_constructor_move() {
                position a(position(0.1, 0.2, 0.3));
                assert(Ffeq(a[0], 0.1));
                assert(Ffeq(a[1], 0.2));
                assert(Ffeq(a[2], 0.3));
            }

            void test_operator_assignment_copy() {
                position a, b(0.1, 0.5, 28.33);
                a = b;
                assert(Ffeq(b[0], a[0]));
                assert(Ffeq(b[1], a[1]));
                assert(Ffeq(b[2], a[2]));

                a[0] = b[0] + 1;
                assert(! Ffeq(b[0], a[0]));
            }

            void test_operator_assignment_move() {
                position a, b(0.1, 0.5, 28.33);
                a = std::move(b);

                assert(Ffeq(b[0], a[0]));
                assert(Ffeq(b[1], a[1]));
                assert(Ffeq(b[2], a[2]));
            }

            void test_operator_equal_int() {
                FPoint<int, 3> a (1, 2, 3);
                FPoint<int, 3> b (1, 2, 3);
                FPoint<int, 3> c (1, 1, 0);
                FPoint<int, 3> d (0, 1, 3);

                assert(a == b);
                assert(b == a);
                assert(a != c);
                assert(c != a);
                assert(a != d);
                assert(d != a);
            }

            void test_operator_equal_float() {
                position a (0.1, 0.2, 0.3);
                position b (0.1, 0.2, 0.3);
                position c (0.1, 0.1, 0);
                position d (0, 0.1, 0.3);

                assert(a == b);
                assert(b == a);
                assert(a != c);
                assert(c != a);
                assert(a != d);
                assert(d != a);
            }

            void test_operator_plus() {
                position a(0.1, 0.2, 0.3);
                position b(0.5, 0.6, 0.7);

                position c = a + b;
                position d;
                d = a + b;
                position e = d - b;

                assert(c == a + b);
                assert(d == a + b);

                assert(a == d - b);
                assert(a == e);

                position f = a + position(1,2,3);
                assert(a != f);
                f -= position(1,2,3);
                assert(a == f);

                position g(a);
                g = position(1,0,0);
                assert(a != g);
            }

            void test_operator_mult() {
                const int mul = 3;
                position a (1,2,3);
                position b (a[0]*mul, a[1]*mul, a[2]*mul);

                position c = a * mul;
                position d = mul * a;

                assert(b == c);
                assert(b == d);
                assert(mul * a == b);

                position e = c / mul;
                assert(a == e);

            }

            void test_iterators() {
                position a(0.1, 0.2, 0.3);
                int i = 0;
                for(auto z : a) {
                    (void)z;
                    i++;
                }
                assert(i == 3);

                i = 0;
                for(auto z = a.rbegin(); z != a.rend(); ++z) {
                    i++;
                }
                assert(i == 3);
            }
        };

    }
}


int main() {
    scalfmm::tests::test_Position().run_all();
}
