#ifndef FCUDAEMPTYCELL_HPP
#define FCUDAEMPTYCELL_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../FStarPUDefaultAlign.hpp"

struct alignas(FStarPUDefaultAlign::StructAlign) FCudaEmptyCell {
    MortonIndex mortonIndex;
    int coordinates[3];
};

#endif // FCUDAEMPTYCELL_HPP

