// /!\ Please, you must read the license at the bottom of this page


/**
 * This file illustrate how to include MPI in your application.
 * 1 - Create an application that inherite from an virtual name (here ApplicationImplementation)
 * 2 - Use the block between == to say what you class will inherit from
 * 3 - Compile as needed
 */

#include "../Src/Utils/FGlobal.hpp"

//================================================================================================
#ifdef SCALFMM_USE_MPI
// Compile by mpic++ testApplication.cpp ../Src/Utils/FAssertable.cpp -o testApplication.exe
// run by mpirun -np 4 ./testApplication.exe
#include "../Src/Utils/FMpiApplication.hpp"
typedef FMpiApplication ApplicationImplementation;
#else
// Compile by g++ testApplication.cpp ../Src/Utils/FAssertable.cpp -o testApplication.exe
#include "../Src/Utils/FSingleApplication.hpp"
typedef FSingleApplication ApplicationImplementation;
#endif
//================================================================================================



#include <stdio.h>


/**
* FApp is an example of the FApplication
* It inherite from ApplicationImplementation
*/
class FApp : public ApplicationImplementation{
public:
    FApp(const int inArgc, char ** const inArgv )
        : ApplicationImplementation(inArgc,inArgv) {
    }

protected:
    void initMaster(){
        printf("I am %d on %d, I am master\n", processId(), processCount());

        const std::string argStr = userParemeterAt<std::string>(0);
        printf("[Master] arg str = %s\n", argStr.c_str());	// will print ./testApplication
        const int argInt = userParemeterAt<int>(0,-1);
        printf("[Master] arg int = %d\n", argInt);		// will print -1
    }
    void initSlave(){
        printf("I am %d on %d, I am slave\n", processId(), processCount());
    }

    void run(){
        printf("I am %d, I start to work\n",processId());
        for(long idx = 0 ; idx < 50000000 ; ++idx) {++idx;--idx;}
        processBarrier();
        printf("I am %d, I just finished\n",processId());
    }
};


// Usual Main
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how mpi MAY be included in the lib (no define!)\n";
    //////////////////////////////////////////////////////////////
    FApp app(argc,argv);
    return app.execute();
}


// [--LICENSE--]
