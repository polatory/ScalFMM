// ===================================================================================
// Copyright ScalFmm 2011 INRIA, Olivier Coulaud, Bérenger Bramas, Matthias Messner
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info".
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FTEMPLATE_HPP
#define FTEMPLATE_HPP

///////////////////////////////////////////////////////////////////////////////////////
/// This file provide useful method to work with template.
/// It provide solution in order to build several methods
/// and run them accordingly to a given condition.
/// We recommand to look at the testTemplateExample.cpp in order
/// to see the usage.
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
/// FForAll : Compile all and exec all
///////////////////////////////////////////////////////////////////////////////////////

#include <functional>

namespace FForAll{

template <class IterType, const IterType CurrentIter, const IterType iterTo, const IterType IterStep,
          class Func, bool IsNotOver, typename... Args>
struct Evaluator{
    static void Run(Args... args){
        Func::template For<CurrentIter>(args...);
        Evaluator<IterType, CurrentIter+IterStep, iterTo, IterStep, Func, (CurrentIter+IterStep < iterTo), Args...>::Run(args...);
    }
};

template <class IterType, const IterType CurrentIter, const IterType iterTo, const IterType IterStep,
          class Func, typename... Args>
struct Evaluator< IterType, CurrentIter, iterTo, IterStep, Func, false, Args...>{
    static void Run(Args... args){
    }
};

template <class IterType, const IterType IterFrom, const IterType iterTo, const IterType IterStep,
          class Func, typename... Args>
void For(Args... args){
    Evaluator<IterType, IterFrom, iterTo, IterStep, Func, (IterFrom<iterTo), Args...>::Run(args...);
}

}


///////////////////////////////////////////////////////////////////////////////////////
/// FRunIf : Compile all and exec only one (if the template variable is equal to
/// the first variable)
///////////////////////////////////////////////////////////////////////////////////////

#include <iostream>

namespace FRunIf{

template <class IterType, const IterType CurrentIter, const IterType iterTo, const IterType IterStep,
          class Func, bool IsNotOver, typename... Args>
struct Evaluator{
    static void Run(IterType value, Args... args){
        if(CurrentIter == value){
            Func::template Run<CurrentIter>(args...);
        }
        else{
            Evaluator<IterType, CurrentIter+IterStep, iterTo, IterStep, Func, (CurrentIter+IterStep < iterTo), Args...>::Run(value, args...);
        }
    }
};

template <class IterType, const IterType CurrentIter, const IterType iterTo, const IterType IterStep,
          class Func, typename... Args>
struct Evaluator< IterType, CurrentIter, iterTo, IterStep, Func, false, Args...>{
    static void Run(IterType value, Args... args){
        std::cout << __FUNCTION__ << " no matching value found\n";
    }
};

template <class IterType, const IterType IterFrom, const IterType iterTo, const IterType IterStep,
          class Func, typename... Args>
void Run(IterType value, Args... args){
    Evaluator<IterType, IterFrom, iterTo, IterStep, Func, (IterFrom<iterTo), Args...>::Run(value, args...);
}

}


///////////////////////////////////////////////////////////////////////////////////////
/// FRunIfFunctional : Compile all and exec only those whose respect a condition
///////////////////////////////////////////////////////////////////////////////////////

namespace FRunIfFunctional{

template <class IterType, const IterType CurrentIter, const IterType iterTo, const IterType IterStep,
          class Func, bool IsNotOver, typename... Args>
struct Evaluator{
    static void Run(std::function<bool(IterType)> test, Args... args){
        if(test(CurrentIter)){
            Func::template Run<CurrentIter>(args...);
        }
        Evaluator<IterType, CurrentIter+IterStep, iterTo, IterStep, Func, (CurrentIter+IterStep < iterTo), Args...>::Run(test, args...);
    }
};

template <class IterType, const IterType CurrentIter, const IterType iterTo, const IterType IterStep,
          class Func, typename... Args>
struct Evaluator< IterType, CurrentIter, iterTo, IterStep, Func, false, Args...>{
    static void Run(std::function<bool(IterType)> test, Args... args){
    }
};

template <class IterType, const IterType IterFrom, const IterType iterTo, const IterType IterStep,
          class Func, typename... Args>
void Run(std::function<bool(IterType)> test,  Args... args){
    Evaluator<IterType, IterFrom, iterTo, IterStep, Func, (IterFrom<iterTo), Args...>::Run(test, args...);
}

}

#endif // FTEMPLATE_HPP
