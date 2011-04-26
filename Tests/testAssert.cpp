// /!\ Please, you must read the license at the bottom of this page


// Compile by : g++ testAssert.cpp ../Src/Utils/FAssertable.cpp -o testAssert.exe

/**
* In this file we show how to use assert and error managing module
*/



#include <iostream>

#include "../Src/Utils/FSingleApplication.hpp"
#include "../Src/Utils/FAssertable.hpp"


// This class is a basic application that need to be assertable
class FApp : public FSingleApplication, protected FAssertable {
public:
	FApp(const int inArgc, char ** const inArgv )
		: FSingleApplication(inArgc,inArgv) {
	}

protected:
	void run(){
		int* pt = new int;
		assert(pt, "pt allocation failled!", __LINE__, __FILE__);
		delete pt;

		assert(false, "Error in doing some stuff", __LINE__, __FILE__);
	}

	// optional : doing something before assert calls FApp->exit
	void calledBeforeExit() const {
		std::cout << "assert is false we will quit! what can I do...?\n";
	}
};

// Usual Main
int main(int argc, char ** argv){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> the FAssert system.\n";
    //////////////////////////////////////////////////////////////
	FApp app(argc,argv);
	return app.execute();
}


// [--LICENSE--]
