// [--License--]

#include "../Src/Utils/FDebug.hpp"

// Compile by : g++ testDebug.cpp ../Src/Utils/FDebug.cpp -o testDebug.exe

/**
* In this file we show how to use the debug module.
* Warning, in FGlobal.hpp (included in FDebug.hpp) SCALFMM_USE_DEBUG might be undefined.
*/

int main(void){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> to understand the FDebug system.\n";
    //////////////////////////////////////////////////////////////
	// Print data simply
	FDEBUG( FDebug::Controller << "Hello Wordl\n");

	// Print a variable (formated print)
	int i = 50;
	FDEBUG( FDebug::Controller.writeVariableFromLine( "i", i, __LINE__, __FILE__););

	// Write a developer information
	FDEBUG( FDebug::Controller.writeFromLine("Strange things are there!", __LINE__, __FILE__); )

	// Flush
	FDEBUG( FDebug::Controller << FDebug::Flush );

	// Change stream type
	FDEBUG( FDebug::Controller.writeToFile("testDebug.out.temp"); )
	FDEBUG( FDebug::Controller << "Hello Wordl 2 the return\n");

	return 0;
}


// [--END--]
