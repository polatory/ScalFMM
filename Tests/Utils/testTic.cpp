// ===================================================================================
// Copyright ScalFmm 2016 INRIA
//
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by Mozilla Public License Version 2.0 (MPL 2.0) and
// abiding by the rules of distribution of free software.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// Mozilla Public License Version 2.0 (MPL 2.0) for more details.
// https://www.mozilla.org/en-US/MPL/2.0/
// ===================================================================================
#include <iostream>
#include "../../Src/Utils/FTic.hpp"
#include "../../Src/Utils/FParameterNames.hpp"

#include <cstdlib>
#include <unistd.h>

/**
* Here we show an example of using FTic
*/

int main(int argc, char ** argv){
    FHelpDescribeAndExit(argc, argv,
                         "Only the code is interesting in order to understand the use of timers.");
    ///////////////////////What we do/////////////////////////////
    std::cout << ">> This executable is useless to execute.\n";
    std::cout << ">> It is only interesting to wath the code to understand\n";
    std::cout << ">> how to use FTic time counter.\n";
    //////////////////////////////////////////////////////////////
    {
	FTic counter;	
	counter.tic();
	usleep(1500000);
	//Sleep(1500); //on windows
	counter.tac();
	std::cout << counter.elapsed() << " (s)\n";
    }
    {
        FTic counter;
        usleep(1500000);
        //Sleep(1500); //on windows
        std::cout << counter.tacAndElapsed() << " (s)\n";
    }
    {
        FTic counter;
        usleep(1500000);
        //Sleep(1500); //on windows
        counter.tac();
        counter.tic();
        usleep(1500000);
        //Sleep(1500); //on windows
        std::cout << counter.tacAndElapsed() << " (s)\n";
        std::cout << counter.cumulated() << " (s)\n";
    }
    return 0;
}

