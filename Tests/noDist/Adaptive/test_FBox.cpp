#include <cassert>
#include <cmath>
#include <iostream>

#include "test_Abstract.hpp"

#include "Utils/FPoint.hpp"
#include "Utils/FConstFuncs.hpp"
#include "Adaptive/new/FBox.hpp"

namespace scalfmm {
    namespace tests {


        struct test_Box : public test_Abstract<test_Box> {

            constexpr static const std::size_t Dim = 3;
            using position_t = FPoint<double, Dim>;
            using box_t      = FBox<position_t>;

            void run_all() {
                RUN(test_constructor_1       );
                RUN(test_constructor_2       );
                RUN(test_constructor_3       );
                RUN(test_constructor_4       );
                RUN(test_constructor_5       );
                RUN(test_assign_copy_operator);
                RUN(test_assign_move_operator);
                RUN(test_center              );
                RUN(test_corner              );
                RUN(test_contains            );
            }


            /// Test default constructor
            void test_constructor_1() {
                box_t b;
                position_t p {0,0,0};
                assert(b.c1 == p);
                assert(b.c2 == p);
            }

            /// Test constructor with minimum and maximum corners
            void test_constructor_2() {
                position_t p1 = {0,0,0};
                position_t p2 = {100,100,100};

                box_t b(p1, p2);
                assert(b.c1 == p1);
                assert(b.c2 == p2);
            }

            /// Test constructor with initializer lists
            void test_constructor_3() {
                position_t p {0,0,0};
                box_t d({0,0,0},{0,0,0});

                assert(d.c1 == p);
                assert(d.c2 == p);
            }

            /// Test copy constructor
            void test_constructor_4() {
                position_t p1 = {0,0,0};
                position_t p2 = {100,100,100};
                box_t c(p1, p2);

                box_t e(c);
                assert(c.c1 == e.c1);
                assert(c.c2 == e.c2);

                // test deep copy
                c.set(p2, p1 + position_t{-1,0,0});

                assert(c.c1 != e.c1);
                assert(c.c2 == e.c2);
            }

            void test_constructor_5() {
                position_t p1 = {0,100,0};
                position_t p2 = {100,0,100};
                box_t c(p1, p2);

                for(std::size_t i = 0; i < Dim; ++i) {
                    assert(Ffeq(c.c1[i], std::fmin(p1[i], p2[i])));
                    assert(Ffeq(c.c2[i], std::fmax(p1[i], p2[i])));
                }
            }

            void test_assign_copy_operator() {
                position_t p1 = {0,0,0};
                position_t p2 = {100,100,100};
                box_t a(p1, p2);

                box_t b;
                b = a;
                assert(a.c1 == b.c1);
                assert(a.c2 == b.c2);
                assert(a.center() == b.center());
            }

            void test_assign_move_operator() {
                position_t p1 = {0,0,0};
                position_t p2 = {100,100,100};
                box_t a(p1, p2);

                box_t b;
                b = std::move(a);
                assert(a.c1 == b.c1);
                assert(a.c2 == b.c2);
            }


            void test_center() {
                position_t p1 = {0,10,50};
                position_t p2 = {100,100,100};
                position_t p3 = {50,55,75};

                box_t b{p1, p2};

                assert(b.center() == p3);

                // Test center after copy
                box_t c(b);
                assert(c.center() == p3);

                b.set({20,30,40}, {120, 130, 140});
                assert(b.center() == position_t(70, 80, 90));
            }

            void test_corner() {
                position_t p0 = {  0,  0,  0};
                position_t p1 = {  0,  0,100};
                position_t p2 = {  0,100,  0};
                position_t p3 = {  0,100,100};
                position_t p4 = {100,  0,  0};
                position_t p5 = {100,  0,100};
                position_t p6 = {100,100,  0};
                position_t p7 = {100,100,100};

                box_t b = {p0, p7};

                assert(b.corner(0) == p0);
                assert(b.corner(1) == p1);
                assert(b.corner(2) == p2);
                assert(b.corner(3) == p3);
                assert(b.corner(4) == p4);
                assert(b.corner(5) == p5);
                assert(b.corner(6) == p6);
                assert(b.corner(7) == p7);

                b.corner(0, {1,1,1});
                assert(b.corner(0) == b.c1);
                assert(b.corner(0) == position_t(1,1,1));
                b.corner(0, {0,0,0});

                b.corner(2, {20, 80, 30});
                assert(b.c1 == position_t(20, 0, 30));
                assert(b.c2 == position_t(100, 80, 100));

                b.corner(5, {0, 20, 25});
                assert(b.c1 == position_t(0, 20, 25));
                assert(b.c2 == position_t(20, 80, 30));
            }

            struct DummyObject {
                position_t p;
                position_t position() const {return p;}
            };

            void test_contains() {
                box_t b = {{0,0,0}, {100,100,100}};

                assert(b.contains(b.c1));
                assert(b.contains(b.c2));
                assert(b.contains({50,50,50}));
                assert(! b.contains({150,50,50}));
                assert(! b.contains({50,150,50}));
                assert(! b.contains({50,50,150}));
                assert(! b.contains({50,150,150}));
                assert(! b.contains({150,150,50}));
                assert(! b.contains({150,50,150}));
                assert(! b.contains({150,150,150}));

                assert(b.contains(DummyObject  {{ 50, 50, 50}}));
                assert(! b.contains(DummyObject{{150, 50, 50}}));
                assert(! b.contains(DummyObject{{ 50,150, 50}}));
                assert(! b.contains(DummyObject{{ 50, 50,150}}));
                assert(! b.contains(DummyObject{{ 50,150,150}}));
                assert(! b.contains(DummyObject{{150,150, 50}}));
                assert(! b.contains(DummyObject{{150, 50,150}}));
                assert(! b.contains(DummyObject{{150,150,150}}));


            }

        };

    }
}

int main() {
    scalfmm::tests::test_Box().run_all();
}
