#ifndef FPARAMETERS_H
#define FPARAMETERS_H

#include <sstream>
#include <iostream>
#include <string.h>

/** This file proposes some methods
  * to work with user input parameters.
  */

namespace FParameters{
    /**
    * This function gives a parameter in a standart type
    * @parameter inArg parameter position has to be strictly less than argc/userParemetersCount
    * @return argv[inArg] in the template VariableType form
    * @warning VariableType need to work with istream >> operator
    * <code> const int argInt = userParemetersAt<int>(1,-1); </code>
    */
    template <class VariableType>
    const VariableType StrToOther(const char* const str, const VariableType& defaultValue = VariableType()){
            std::istringstream iss(str,std::istringstream::in);
            VariableType value;
            iss >> value;
            if( /*iss.tellg()*/ iss.eof() ) return value;
            return defaultValue;
    }

    /** To put a char into lower format
      *
      */
    char toLower(const char c){
        return char('A' <= c && c <= 'Z' ? (c - 'A') + 'a' : c);
    }

    /** To know if two char are equals
      *
      */
    bool areCharsEquals(const char c1, const char c2, const bool caseSensible = false){
        return (caseSensible && c1 == c2) || (!caseSensible && toLower(c1) == toLower(c2));
    }

    /** To know if two str are equals
      *
      */
    bool areStrEquals(const char* const inStr1, const char* const inStr2, const bool caseSensible = false){
        int idxStr = 0;
        while(inStr1[idxStr] != '\0' && inStr2[idxStr] != '\0'){
            if(!areCharsEquals(inStr1[idxStr] ,inStr2[idxStr],caseSensible)){
                return false;
            }
            ++idxStr;
        }
        return inStr1[idxStr] == inStr2[idxStr];
    }

    /** To find a parameters from user format char parameters
      *
      */
    int findParameter(const int argc, const char* const * const argv, const char* const inName, const bool caseSensible = false){
        for(int idxArg = 0; idxArg < argc ; ++idxArg){
            if(areStrEquals(inName, argv[idxArg], caseSensible)){
                return idxArg;
            }
        }
        return -1;
    }

    /** To get a value like :
      * getValue(argc,argv, "Toto", 0, false);
      * will return 55 if the command contains : -Toto 55
      * else 0
      */
    template <class VariableType>
    const VariableType getValue(const int argc, const char* const * const argv, const char* const inName, const VariableType& defaultValue = VariableType(), const bool caseSensible = false){
        const int position = findParameter(argc,argv,inName,caseSensible);
        if(position == -1 || position == argc - 1){
            return defaultValue;
        }
        return StrToOther(argv[position+1],defaultValue);
    }

    /** Get a str from argv
      */
    const char* getStr(const int argc, const char* const * const argv, const char* const inName, const char* const inDefault, const bool caseSensible = false){
        const int position = findParameter(argc,argv,inName,caseSensible);
        if(position == -1 || position == argc - 1){
            return inDefault;
        }
        return argv[position+1];
    }
}



#endif // FPARAMETERS_H
