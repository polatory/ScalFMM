#ifndef FHUTILS_HPP
#define FHUTILS_HPP

#include "Utils/FGlobal.hpp"

#include <cstring>

template <class Type>
void FSetToZeros(Type array[], const int length){
    memset(array, 0, length*sizeof(Type));
}


struct FBlockDescriptor {
    int row, col, nbRows, nbCols, level;
};

#endif // FHUTILS_HPP

