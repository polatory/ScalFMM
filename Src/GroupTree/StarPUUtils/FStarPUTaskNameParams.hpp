#ifndef FSTARPUTASKNAMEPARAMS_HPP
#define FSTARPUTASKNAMEPARAMS_HPP

#include "../../Utils/FGlobal.hpp"

#include <list>
#include <cstring>
#include <cstdio>

/**
 * This class creates task name for starpu
 * it is used for simgrid (to pass task parameters)
 */
class FStarPUTaskNameParams{
protected:
    std::list<const char*> names;

public:
    FStarPUTaskNameParams(){
    }

    ~FStarPUTaskNameParams(){
        clear();
    }

    void clear(){
        while(names.size()){
            delete[] names.front();
            names.pop_front();
        }
    }

    template <typename ... Params>
    const char* print(const char format[], Params... args ){
        const size_t length = 512;
        char* name = new char[length+1];
        snprintf(name, length, format, args...);
        name[length] = '\0';
        names.push_back(name);
        return name;
    }

    const char* add(const char* strToCpy){
        const size_t length = strlen(strToCpy);
        char* cpy = new char[length+1];
        memcpy(cpy, strToCpy, length);
        cpy[length] = '\0';
        names.push_back(cpy);
        return cpy;
    }
};

#endif // FSTARPUTASKNAMEPARAMS_HPP
