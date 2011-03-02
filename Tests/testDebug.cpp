// /!\ Please, you must read the license at the bottom of this page

#include "../Sources/Utils/FDebug.hpp"

// Compile by : g++ testDebug.cpp ../Sources/Utils/FDebug.cpp -o testFDebug.exe

/**
* In this file we show how to use the debug module
* please refere to the source of testDebug.cpp directly to knwo more
*/

int main(void){
	// Print data simply
	FDEBUG( FDebug::Controller << "Hello Wordl\n");

	// Print a variable (formated print)
	int i = 50;
	FDEBUG( FDebug::Controller.writeVariableFromLine( "i", i, __LINE__, __FILE__););

	// Write a developer information
	FDEBUG( FDebug::Controller.writeFromLine("Strange things happend here!", __LINE__, __FILE__); )

	// Change stream type
	FDEBUG( FDebug::Controller.writeToFile("FDebug.out"); )
	FDEBUG( FDebug::Controller << "Hello Wordl 2 the return\n");

	return 0;
}


// [--LICENSE--]
