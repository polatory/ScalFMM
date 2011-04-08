#include <iostream>
#include "../Sources/Utils/FTic.hpp"

#include <stdlib.h>
#include <unistd.h>

/**
* Here we show an example of using FTic
* g++ testTic.cpp -o testTic.exe
*/

int main(){
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use FTic time counter.\n";
    //////////////////////////////////////////////////////////////

	FTic counter;	

	counter.tic();
	usleep(1500000);
	//Sleep(1500); //on windows
	counter.tac();

	std::cout << counter.elapsed() << " (s)\n";

	return 0;
}

