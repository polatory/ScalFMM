// /!\ Please, you must read the license at the bottom of this page

#include <iostream>

#include "../Sources/Utils/FSingleApplication.hpp"
#include "../Sources/Utils/FAssertable.hpp"

// Compile by : g++ testAssert.cpp ../Sources/Utils/FAssertable.cpp -o testAssert.exe

/**
* In this file we show how to use assert and error managing module
*/

// This class is a basic application that need to be assertable
class FApp : public FSingleApplication, public FAssertable {
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
	FApp app(argc,argv);
	return app.execute();
}


// [--LICENSE--]
