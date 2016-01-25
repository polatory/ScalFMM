#ifndef _SCALFMM_TEST_MOCKPARTICLE_HPP_
#define _SCALFMM_TEST_MOCKPARTICLE_HPP_

#include "Utils/FPoint.hpp"

namespace scalfmm {
    namespace tests {

        template<typename _Real, std::size_t _Dim>
        struct MockParticle {

            constexpr static const std::size_t Dim = _Dim;
            using Real = _Real;
            using position_t =  FPoint<Real, Dim>;

            static std::size_t idx;

            MockParticle(position_t pos) :
                _pos(pos),
                index(++idx)
            {}

            MockParticle(position_t pos, std::size_t new_index) :
                _pos(pos),
                index(new_index)
            {
                idx = new_index+1;
            }


            position_t _pos;
            std::size_t index;

            position_t& position(const position_t& p) {return (_pos = p);}
            position_t& position() {return _pos;}
            const position_t& position() const {return _pos;}
        };

        template<typename Real, std::size_t Dim>
        std::size_t MockParticle<Real, Dim>::idx = 0;
    }
}


#endif
