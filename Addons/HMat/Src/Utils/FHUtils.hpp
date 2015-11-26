#ifndef FHUTILS_HPP
#define FHUTILS_HPP

#include "Utils/FGlobal.hpp"

#include <cstring>

template <class Type>
void FSetToZeros(Type array[], const int length){
    memset(array, 0, length*sizeof(Type));
}

#endif // FHUTILS_HPP

