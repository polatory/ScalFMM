#ifndef SCALFMM_TEST_ABSTRACT_HPP_
#define SCALFMM_TEST_ABSTRACT_HPP_

#include <type_traits>
#include <iostream>

#define RUN(test_function_name)                                         \
    std::cout << #test_function_name << std::endl;                      \
    run_test(& std::remove_reference<decltype((*this))>::type :: test_function_name);

template<class T>
struct test_Abstract {
    virtual void set_up() {}
    virtual void tear_down() {}
    virtual void run_all() = 0;

    virtual void run_test(void (T::*test_function)()) {
        this->set_up();
        (dynamic_cast<T*>(this)->*test_function)();
        this->tear_down();

    }
};

#endif
